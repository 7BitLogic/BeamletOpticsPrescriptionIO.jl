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
        lens = SphericalLens(el.front_surface.radius, el.back_surface.radius, el.center_thickness, el.diameter, n)
    else
        s_front = build_surface(el.front_surface, el.diameter)
        s_back = build_surface(el.back_surface, el.diameter)
        lens = Lens(s_front, s_back, el.center_thickness, n)
    end
    
    translate3d!(lens, [0.0, el.axial_position, 0.0])
    return (lens, false)
end

function build_bmo_object(el::ZmxDoublet)
    n1 = resolve_refractive_index(el.glass1, el.nd1, el.vd1)
    n2 = resolve_refractive_index(el.glass2, el.nd2, el.vd2)
    
    dl = SphericalDoubletLens(
        el.surface1.radius, el.surface2.radius, el.surface3.radius,
        el.thickness1, el.thickness2, el.diameter, n1, n2
    )
    translate3d!(dl, [0.0, el.axial_position, 0.0])
    return (dl, false)
end

function build_bmo_object(el::ZmxTriplet)
    n1 = resolve_refractive_index(el.glass1)
    n2 = resolve_refractive_index(el.glass2)
    n3 = resolve_refractive_index(el.glass3)
    
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
"""
function build_bmo_system(zmx_sys::ZmxSystem; elements::Vector{ZmxElement} = group_elements(zmx_sys))::BMOImportResult
    optical_objects = []
    local detector = nothing

    for el in elements
        obj, is_det = build_bmo_object(el)
        if obj !== nothing
            push!(optical_objects, obj)
            if is_det
                detector = obj
            end
        end
    end

    system = System(optical_objects)
    return BMOImportResult(system, detector, elements, zmx_sys)
end

"""
    import_zmx(filepath::AbstractString)::BMOImportResult

Main entry point: loads a Zemax .zmx file and builds the reconstructed `BeamletOptics.System`.
"""
function import_zmx(filepath::AbstractString)::BMOImportResult
    zmx_sys = parse_zmx(filepath)
    elements = group_elements(zmx_sys)
    return build_bmo_system(zmx_sys; elements = elements)
end

