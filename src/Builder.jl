"""
Builder constructs BeamletOptics.jl optical objects from parsed ZMX elements.
"""

using BeamletOptics
using LinearAlgebra

"""
    build_surface(s::ZmxSurface, diameter::Real)

Instantiates the appropriate BMO surface geometry:
- `CircularFlatSurface` for plano surfaces (R = Inf)
- `SphericalSurface` for spherical surfaces
- `EvenAsphericalSurface` for even aspheric surfaces
"""
function build_surface(s::ZmxSurface, diameter::Real)
    if s.surface_type == :EVENASPH
        coeffs = isempty(s.parms) ? Float64[] : s.parms
        return EvenAsphericalSurface(s.radius, Float64(diameter), s.conic, coeffs)
    elseif isinf(s.radius)
        return CircularFlatSurface(Float64(diameter))
    else
        return SphericalSurface(s.radius, Float64(diameter))
    end
end

"""
    build_bmo_object(el::ZmxElement)

Instantiates a concrete 3D BeamletOptics object for a given `ZmxElement` and translates it
to its absolute position along the optical axis (+Y).
Returns `(object, is_detector)`.
"""
function build_bmo_object(el::ZmxSinglet)
    n = resolve_refractive_index(el.glass_name, el.nd, el.vd)
    
    is_front_asph = el.front_surface.surface_type == :EVENASPH
    is_back_asph = el.back_surface.surface_type == :EVENASPH
    
    local lens
    if !is_front_asph && !is_back_asph
        try
            lens = SphericalLens(el.front_surface.radius, el.back_surface.radius, el.center_thickness, el.diameter, n)
        catch e
            lens = ThinLens(el.front_surface.radius, el.back_surface.radius, el.diameter, n)
        end
    else
        try
            s_front = build_surface(el.front_surface, el.diameter)
            s_back = build_surface(el.back_surface, el.diameter)
            lens = Lens(s_front, s_back, el.center_thickness, n)
        catch e
            lens = ThinLens(el.front_surface.radius, el.back_surface.radius, el.diameter, n)
        end
    end
    
    translate3d!(lens, [0.0, el.axial_position, 0.0])
    return (lens, false)
end

function build_bmo_object(el::ZmxDoublet)
    n1 = resolve_refractive_index(el.glass1, el.nd1, el.vd1)
    n2 = resolve_refractive_index(el.glass2, el.nd2, el.vd2)
    
    local dl
    try
        dl = SphericalDoubletLens(
            el.surface1.radius, el.surface2.radius, el.surface3.radius,
            el.thickness1, el.thickness2, el.diameter, n1, n2
        )
    catch e
        # Fallback to chained thin lenses
        front = ThinLens(el.surface1.radius, el.surface2.radius, el.diameter, n1)
        back = ThinLens(el.surface2.radius, el.surface3.radius, el.diameter, n2)
        translate3d!(back, [0.0, el.thickness1, 0.0])
        dl = DoubletLens(front, back)
    end
    translate3d!(dl, [0.0, el.axial_position, 0.0])
    return (dl, false)
end

function build_bmo_object(el::ZmxTriplet)
    n1 = resolve_refractive_index(el.glass1, el.nd1, el.vd1)
    n2 = resolve_refractive_index(el.glass2, el.nd2, el.vd2)
    n3 = resolve_refractive_index(el.glass3, el.nd3, el.vd3)
    
    tl = SphericalTripletLens(
        el.surface1.radius, el.surface2.radius, el.surface3.radius, el.surface4.radius,
        el.thickness1, el.thickness2, el.thickness3, el.diameter,
        n1, n2, n3
    )
    translate3d!(tl, [0.0, el.axial_position, 0.0])
    return (tl, false)
end

function build_bmo_object(el::ZmxMirror)
    sub_thick = max(el.diameter * 0.1, 0.005)
    local mirror
    if isinf(el.surface.radius)
        mirror = RoundPlanoMirror(el.diameter, sub_thick)
    else
        mirror = SphericalMirror(el.surface.radius, sub_thick, el.diameter)
    end
    translate3d!(mirror, [0.0, el.axial_position, 0.0])
    return (mirror, false)
end

function build_bmo_object(el::ZmxDetector)
    det = Detector(el.diameter)
    translate3d!(det, [0.0, el.axial_position, 0.0])
    return (det, true)
end

function build_bmo_object(el::ZmxStop)
    return (nothing, false)
end

"""
    build_bmo_system(zmx_sys::ZmxSystem; elements::Vector{ZmxElement} = group_elements(zmx_sys))::BMOImportResult

Builds a complete `BeamletOptics.System` and returns a `BMOImportResult`.
For folded/catadioptric systems where a detector is placed in front of a mirror, an absorber
is placed behind the detector to intercept incoming forward rays (central obscuration).
"""
function build_bmo_system(zmx_sys::ZmxSystem; elements::Vector{ZmxElement} = group_elements(zmx_sys))::BMOImportResult
    optical_objects = []
    local detector = nothing
    local stop_elem = nothing

    # Detect if this is a folded/catadioptric system with a mirror reflecting backwards
    has_mirror = any(el -> el isa ZmxMirror, elements)
    mirror_pos = has_mirror ? maximum(el.axial_position for el in elements if el isa ZmxMirror) : 0.0

    for el in elements
        if el isa ZmxStop
            stop_elem = el
            continue
        end

        obj, is_det = build_bmo_object(el)
        if obj !== nothing
            if is_det
                detector = obj
                # If folded catadioptric system with detector in front of mirror:
                if has_mirror && el.axial_position < mirror_pos
                    # Place an absorbing disc behind the detector facing incoming rays (-Y)
                    # to intercept forward rays (camera obscuration)
                    gap = min(0.001, (mirror_pos - el.axial_position) * 0.02)
                    absorber = Detector(el.diameter, true)
                    translate3d!(absorber, [0.0, el.axial_position - gap, 0.0])
                    push!(optical_objects, absorber)

                    # Rotate detector by π around X so its active surface faces the mirror (+Y)
                    rotate3d!(detector, [1.0, 0.0, 0.0], π)
                end
            end
            push!(optical_objects, obj)
        end
    end

    system = System(optical_objects)
    return BMOImportResult(system, detector, stop_elem, elements, zmx_sys)
end

"""
    import_zmx(filepath::AbstractString; configuration::Int=1)::BMOImportResult

Main entry point: loads a Zemax .zmx file and builds the reconstructed `BeamletOptics.System`.
Applies multi-configuration overrides (`THIC`) for the specified `configuration` (default: 1).
"""
function import_zmx(filepath::AbstractString; configuration::Int=1)::BMOImportResult
    zmx_sys = parse_zmx(filepath; configuration=configuration)
    elements = group_elements(zmx_sys)
    return build_bmo_system(zmx_sys; elements = elements)
end

"""
    suggest_source(res::BMOImportResult; num_rays::Int=200, num_rings::Int=6)::AbstractBeamGroup

Generates a beam source for the optical system with the primary wavelength.
Automatically discriminates between:
- **Finite object conjugate** (`SURF 0 DISZ < ∞`): Creates a `PointSource` located at the object distance
  aimed at the entrance pupil / aperture stop with numerical aperture from `OBNA` or pupil geometry.
- **Infinite object conjugate** (`SURF 0 DISZ == ∞`): Creates an on-axis `CollimatedSource` in front of the system.
"""
function suggest_source(res::BMOImportResult; num_rays::Int=200, num_rings::Int=6)::BeamletOptics.AbstractBeamGroup
    # BMO requirement: num_rays >= 20 * num_rings
    num_rays = max(num_rays, 20 * num_rings)

    # 1. Determine axial position of front optical component and aperture stop
    y_front = Inf
    for el in res.elements
        if !(el isa ZmxStop)
            y_front = min(y_front, el.axial_position)
        end
    end
    if isinf(y_front)
        y_front = 0.0
    end
    y_pupil = res.stop !== nothing ? res.stop.axial_position : y_front

    # 2. Determine pupil/beam diameter
    D = 0.0
    if res.zmx_system.enpd > 0.0
        D = res.zmx_system.enpd
    elseif res.stop !== nothing && res.stop.diameter > 0.0
        D = res.stop.diameter
    elseif !isempty(res.elements) && hasproperty(res.elements[1], :diameter)
        D = res.elements[1].diameter * 0.8
    else
        D = 0.02 # fallback 20 mm
    end

    # 3. Primary wavelength
    wvl_idx = clamp(res.zmx_system.primary_wavelength_idx, 1, length(res.zmx_system.wavelengths))
    λ = res.zmx_system.wavelengths[wvl_idx]

    # 4. Check finite vs. infinite conjugate
    is_finite = isfinite(res.zmx_system.object_distance) && res.zmx_system.object_distance > 0.0

    if is_finite
        y_obj = y_front - res.zmx_system.object_distance
        θ = if res.zmx_system.obna > 0.0
            asin(clamp(res.zmx_system.obna, 0.0, 0.999))
        else
            dist_to_pupil = max(1e-6, abs(y_pupil - y_obj))
            atan((D / 2.0) / dist_to_pupil)
        end
        return PointSource([0.0, y_obj, 0.0], [0.0, 1.0, 0.0], θ, λ;
                           num_rings = num_rings, num_rays = num_rays)
    else
        y_start = y_front - 0.010 # 10 mm in front of first component
        return CollimatedSource(
            [0.0, y_start, 0.0],
            [0.0, 1.0, 0.0],
            D,
            λ;
            num_rings = num_rings,
            num_rays = num_rays
        )
    end
end

"""
    generate_zmx_rays(
        res::BMOImportResult;
        fields = :all,
        wavelengths = :all,
        num_rings::Int = 4,
        num_rays::Int = 80,
        combined::Bool = false
    )

Reproduces the full ray bundles defined in the Zemax file across specified fields (`fields`)
and wavelengths (`wavelengths`).

Automatically respects:
- **Finite object conjugates** (`SURF 0 DISZ < ∞`): Spawns `PointSource` objects from the object plane.
- **Infinite object conjugates** (`SURF 0 DISZ == ∞`): Spawns `CollimatedSource` objects under the field angles.
- **Field types**: `:angle` (degrees), `:object_height`, `:paraxial_image_height`, `:real_image_height` (meters).
- **Aperture definition**: `OBNA` or entrance pupil / aperture stop geometry.

Keyword arguments:
- `fields`: `:all`, an integer field index (e.g. `1`), or a vector of indices.
- `wavelengths`: `:all`, `:primary`, an integer wavelength index, or a vector of indices.
- `combined`: If `true`, groups all generated beams into a single source object. If `false`, returns `Vector{AbstractBeamGroup}`.
"""
function generate_zmx_rays(
    res::BMOImportResult;
    fields = :all,
    wavelengths = :all,
    num_rings::Int = 4,
    num_rays::Int = 80,
    combined::Bool = false
)
    num_rays = max(num_rays, 20 * num_rings)

    # 1. Front component and pupil axial locations
    y_front = Inf
    for el in res.elements
        if !(el isa ZmxStop)
            y_front = min(y_front, el.axial_position)
        end
    end
    if isinf(y_front)
        y_front = 0.0
    end
    y_pupil = res.stop !== nothing ? res.stop.axial_position : y_front

    # 2. Pupil / beam diameter
    D = 0.0
    if res.zmx_system.enpd > 0.0
        D = res.zmx_system.enpd
    elseif res.stop !== nothing && res.stop.diameter > 0.0
        D = res.stop.diameter
    elseif !isempty(res.elements) && hasproperty(res.elements[1], :diameter)
        D = res.elements[1].diameter * 0.8
    else
        D = 0.02
    end

    # 3. Selected wavelengths
    all_wvls = res.zmx_system.wavelengths
    selected_wvls = if wavelengths === :all
        all_wvls
    elseif wavelengths === :primary
        [all_wvls[clamp(res.zmx_system.primary_wavelength_idx, 1, length(all_wvls))]]
    elseif wavelengths isa Integer
        [all_wvls[clamp(wavelengths, 1, length(all_wvls))]]
    elseif wavelengths isa AbstractVector
        [all_wvls[clamp(i, 1, length(all_wvls))] for i in wavelengths]
    else
        [Float64(wavelengths)]
    end

    # 4. Selected fields
    all_fields = isempty(res.zmx_system.fields) ? [(0.0, 0.0, 1.0)] : res.zmx_system.fields
    selected_fields = if fields === :all
        all_fields
    elseif fields isa Integer
        [all_fields[clamp(fields, 1, length(all_fields))]]
    elseif fields isa AbstractVector
        [all_fields[clamp(i, 1, length(all_fields))] for i in fields]
    else
        [(0.0, 0.0, 1.0)]
    end

    # 5. Check finite vs. infinite conjugate
    is_finite = isfinite(res.zmx_system.object_distance) && res.zmx_system.object_distance > 0.0
    y_obj = is_finite ? (y_front - res.zmx_system.object_distance) : (y_front - 0.010)

    sources = BeamletOptics.AbstractBeamGroup[]
    all_beams = Beam[]

    for (fx, fy, fw) in selected_fields
        if is_finite
            # Finite conjugate: Field determines position on the object plane
            x_obj = if res.zmx_system.field_type == :angle
                tan(deg2rad(fx)) * res.zmx_system.object_distance
            else
                fx # Already scaled to meters in Parser
            end
            z_obj = if res.zmx_system.field_type == :angle
                tan(deg2rad(fy)) * res.zmx_system.object_distance
            else
                fy
            end
            P_obj = [x_obj, y_obj, z_obj]

            # Chief ray direction towards center of pupil
            vec_to_pupil = [0.0 - x_obj, y_pupil - y_obj, 0.0 - z_obj]
            L = norm(vec_to_pupil)
            dir = normalize(vec_to_pupil)

            # Cone half-angle
            θ = if res.zmx_system.obna > 0.0
                asin(clamp(res.zmx_system.obna, 0.0, 0.999))
            else
                atan((D / 2.0) / max(1e-6, L))
            end

            for λ in selected_wvls
                src = PointSource(P_obj, dir, θ, λ; num_rings=num_rings, num_rays=num_rays)
                push!(sources, src)
                append!(all_beams, src.beams)
            end
        else
            # Infinite conjugate: Field determines beam incident angle
            if res.zmx_system.field_type == :angle
                ang_x = deg2rad(fx)
                ang_z = deg2rad(fy)
                kx = tan(ang_x)
                kz = tan(ang_z)
                dir = normalize([kx, 1.0, kz])
            else
                dir = [0.0, 1.0, 0.0]
            end

            for λ in selected_wvls
                center = [0.0, y_obj, 0.0]
                src = CollimatedSource(center, dir, D, λ; num_rings=num_rings, num_rays=num_rays)
                push!(sources, src)
                append!(all_beams, src.beams)
            end
        end
    end

    if combined
        return CollimatedSource(all_beams, D, [0.0, y_obj, 0.0], [0.0, 1.0, 0.0])
    else
        return sources
    end
end

"""
    find_best_focus(detector::Detector) -> NamedTuple

Analytically determines the optimal focal plane displacement along the optical axis (Y)
that minimizes the geometric RMS spot radius on the detector.

Returns a named tuple with:
- `delta_y`: Axial displacement from the detector plane (in meters)
- `y_opt`: New axial position `y_nom + delta_y` (in meters)
- `rms_opt`: Geometric RMS spot radius at optimal focus (in meters)
- `rms_nom`: Geometric RMS spot radius at the current detector plane (in meters)
"""
function find_best_focus(detector::Detector)
    hits = detector.hits
    if hits === nothing || isempty(hits)
        return (delta_y=0.0, y_opt=0.0, rms_opt=0.0, rms_nom=0.0)
    end
    N = length(hits)
    pts = [BeamletOptics.hit_point(h) for h in hits]
    dirs = [h.ray.dir for h in hits]

    xs = [p[1] for p in pts]
    ys = [p[2] for p in pts]
    zs = [p[3] for p in pts]

    y_nom = ys[1]

    us = [d[1] for d in dirs]
    vs = [d[2] for d in dirs]
    ws = [d[3] for d in dirs]

    # Ray slopes relative to optical axis Y
    slope_x = us ./ vs
    slope_z = ws ./ vs

    cx0 = sum(xs) / N
    cz0 = sum(zs) / N
    csx = sum(slope_x) / N
    csz = sum(slope_z) / N

    tx = xs .- cx0
    tz = zs .- cz0
    tu = slope_x .- csx
    tw = slope_z .- csz

    A = sum(tu.^2 .+ tw.^2) / N
    B = sum(tx .* tu .+ tz .* tw) / N
    C = sum(tx.^2 .+ tz.^2) / N

    rms_nom = sqrt(max(0.0, C))
    if A <= 1e-18
        return (delta_y=0.0, y_opt=y_nom, rms_opt=rms_nom, rms_nom=rms_nom)
    end

    delta_y = -B / A
    y_opt = y_nom + delta_y
    rms_opt_sq = C - (B^2 / A)
    rms_opt = rms_opt_sq > 0 ? sqrt(rms_opt_sq) : 0.0

    return (delta_y=delta_y, y_opt=y_opt, rms_opt=rms_opt, rms_nom=rms_nom)
end

find_best_focus(res::BMOImportResult) = res.detector !== nothing ? find_best_focus(res.detector) : (delta_y=0.0, y_opt=0.0, rms_opt=0.0, rms_nom=0.0)

"""
    refocus!(res::BMOImportResult; delta_y=nothing, y_opt=nothing)

Shifts the detector in `res.system` along the optical axis (Y).
If neither `delta_y` nor `y_opt` is specified, it analytically calculates `find_best_focus(res.detector)`
to find the optimal focal plane automatically and shifts the detector to it.
Returns the applied `delta_y` in meters.
"""
function refocus!(res::BMOImportResult; delta_y=nothing, y_opt=nothing)
    res.detector === nothing && return 0.0
    shift = 0.0
    if delta_y !== nothing
        shift = Float64(delta_y)
    elseif y_opt !== nothing
        cur_y = position(res.detector)[2]
        shift = Float64(y_opt) - cur_y
    else
        opt = find_best_focus(res.detector)
        shift = opt.delta_y
    end
    translate3d!(res.detector, [0.0, shift, 0.0])
    return shift
end


