"""
ElementGrouper converts sequential Zemax surfaces into discrete optical elements
(Singlet lenses, Cemented doublets/triplets, Mirrors, Stops, and Image detectors).
"""

function is_air(glass_name::AbstractString)::Bool
    g = uppercase(strip(glass_name))
    return isempty(g) || g == "AIR" || g == "0"
end

function is_mirror_glass(glass_name::AbstractString)::Bool
    return uppercase(strip(glass_name)) == "MIRROR"
end

"""
    safe_lens_diameter(nominal_diam::Float64, radii::Vector{Float64})::Float64

Ensures diameter is physically consistent with surface radii of curvature.
In BMO, a spherical surface requires radius >= diameter / 2.
"""
function safe_lens_diameter(nominal_diam::Float64, radii::Vector{Float64})::Float64
    d = nominal_diam
    if d <= 0.0
        # Estimate from finite radii or fallback to 25.4mm
        finite_r = [abs(r) for r in radii if !isinf(r) && abs(r) > 1e-6]
        d = isempty(finite_r) ? 25.4e-3 : minimum(finite_r) * 0.8
    end
    
    # Check against spherical radius limits
    for r in radii
        if !isinf(r) && abs(r) > 1e-6
            max_d = 1.95 * abs(r)
            if d > max_d
                d = max_d
            end
        end
    end
    return d
end

"""
    group_elements(zmx_sys::ZmxSystem)::Vector{ZmxElement}

Parses the sequential surface stream into discrete `ZmxElement` components, computing
absolute axial vertex coordinates along the BMO optical axis (+Y).
"""
function group_elements(zmx_sys::ZmxSystem)::Vector{ZmxElement}
    surfaces = zmx_sys.surfaces
    n_surfs = length(surfaces)
    n_surfs == 0 && return ZmxElement[]

    start_idx = 1
    if surfaces[1].index == 0
        start_idx = 2
    end

    pos_y = zeros(Float64, n_surfs)
    if start_idx <= n_surfs
        pos_y[start_idx] = 0.0
        for k in start_idx:(n_surfs - 1)
            t = isinf(surfaces[k].thickness) ? 0.0 : surfaces[k].thickness
            pos_y[k + 1] = pos_y[k] + t
        end
    end

    elements = ZmxElement[]
    elem_idx = 1
    i = start_idx

    while i < n_surfs
        surf = surfaces[i]
        
        # 1. Mirror check
        if surf.is_mirror || is_mirror_glass(surf.glass_name)
            nominal_d = 2.0 * surf.semi_diameter
            diam = safe_lens_diameter(nominal_d, [surf.radius])
            push!(elements, ZmxMirror(
                name = isempty(surf.comment) ? "Mirror_$elem_idx" : surf.comment,
                surface = surf,
                diameter = diam,
                axial_position = pos_y[i]
            ))
            elem_idx += 1
            i += 1
            continue
        end

        # 2. Stop check (if stop is situated in an air gap)
        if surf.is_stop && is_air(surf.glass_name)
            diam = max(2.0 * surf.semi_diameter, 0.001)
            push!(elements, ZmxStop(
                name = isempty(surf.comment) ? "Stop" : surf.comment,
                surface = surf,
                diameter = diam,
                axial_position = pos_y[i]
            ))
            i += 1
            continue
        end

        # 3. Refractive glass elements (Singlet, Doublet, Triplet)
        if !is_air(surf.glass_name)
            k = i
            while k < n_surfs && !is_air(surfaces[k].glass_name) && !surfaces[k].is_mirror && !is_mirror_glass(surfaces[k].glass_name)
                k += 1
            end
            num_glass = k - i

            if num_glass == 1
                s_front = surfaces[i]
                s_back = surfaces[i + 1]
                cthick = s_front.thickness
                nominal_d = max(2.0 * s_front.semi_diameter, 2.0 * s_back.semi_diameter)
                diam = safe_lens_diameter(nominal_d, [s_front.radius, s_back.radius])
                name = isempty(s_front.comment) ? "Lens_$elem_idx" : s_front.comment
                
                push!(elements, ZmxSinglet(
                    name = name,
                    front_surface = s_front,
                    back_surface = s_back,
                    center_thickness = cthick,
                    diameter = diam,
                    glass_name = s_front.glass_name,
                    nd = s_front.nd,
                    vd = s_front.vd,
                    axial_position = pos_y[i]
                ))
                elem_idx += 1
                i = i + 1
            elseif num_glass == 2
                s1 = surfaces[i]
                s2 = surfaces[i + 1]
                s3 = surfaces[i + 2]
                nominal_d = max(2.0 * s1.semi_diameter, 2.0 * s2.semi_diameter, 2.0 * s3.semi_diameter)
                diam = safe_lens_diameter(nominal_d, [s1.radius, s2.radius, s3.radius])
                name = isempty(s1.comment) ? "Doublet_$elem_idx" : s1.comment
                
                push!(elements, ZmxDoublet(
                    name = name,
                    surface1 = s1,
                    surface2 = s2,
                    surface3 = s3,
                    thickness1 = s1.thickness,
                    thickness2 = s2.thickness,
                    diameter = diam,
                    glass1 = s1.glass_name,
                    glass2 = s2.glass_name,
                    nd1 = s1.nd,
                    vd1 = s1.vd,
                    nd2 = s2.nd,
                    vd2 = s2.vd,
                    axial_position = pos_y[i]
                ))
                elem_idx += 1
                i = i + 2
            elseif num_glass == 3
                s1 = surfaces[i]
                s2 = surfaces[i + 1]
                s3 = surfaces[i + 2]
                s4 = surfaces[i + 3]
                nominal_d = max(2.0 * s1.semi_diameter, 2.0 * s2.semi_diameter, 2.0 * s3.semi_diameter, 2.0 * s4.semi_diameter)
                diam = safe_lens_diameter(nominal_d, [s1.radius, s2.radius, s3.radius, s4.radius])
                name = isempty(s1.comment) ? "Triplet_$elem_idx" : s1.comment
                
                push!(elements, ZmxTriplet(
                    name = name,
                    surface1 = s1,
                    surface2 = s2,
                    surface3 = s3,
                    surface4 = s4,
                    thickness1 = s1.thickness,
                    thickness2 = s2.thickness,
                    thickness3 = s3.thickness,
                    diameter = diam,
                    glass1 = s1.glass_name,
                    glass2 = s2.glass_name,
                    glass3 = s3.glass_name,
                    axial_position = pos_y[i]
                ))
                elem_idx += 1
                i = i + 3
            else
                s1 = surfaces[i]
                s2 = surfaces[i + 1]
                s3 = surfaces[i + 2]
                nominal_d = max(2.0 * s1.semi_diameter, 2.0 * s2.semi_diameter, 2.0 * s3.semi_diameter)
                diam = safe_lens_diameter(nominal_d, [s1.radius, s2.radius, s3.radius])
                push!(elements, ZmxDoublet(
                    name = "Cemented_$elem_idx",
                    surface1 = s1,
                    surface2 = s2,
                    surface3 = s3,
                    thickness1 = s1.thickness,
                    thickness2 = s2.thickness,
                    diameter = diam,
                    glass1 = s1.glass_name,
                    glass2 = s2.glass_name,
                    nd1 = s1.nd,
                    vd1 = s1.vd,
                    nd2 = s2.nd,
                    vd2 = s2.vd,
                    axial_position = pos_y[i]
                ))
                elem_idx += 1
                i = i + 2
            end
            continue
        end

        i += 1
    end

    if n_surfs >= start_idx
        img_surf = surfaces[end]
        det_diam = max(2.0 * img_surf.semi_diameter, 0.005)
        push!(elements, ZmxDetector(
            name = "Detector",
            surface = img_surf,
            diameter = det_diam,
            axial_position = pos_y[end]
        ))
    end

    return elements
end

