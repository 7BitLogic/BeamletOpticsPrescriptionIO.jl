"""
ElementGrouper converts sequential Zemax surfaces into discrete optical elements
(Singlet lenses, Cemented doublets/triplets, Mirrors, Stops, and Image detectors).
"""

"""
    is_air(glass_name::AbstractString)::Bool

Returns true if the glass represents air or free space.
"""
function is_air(glass_name::AbstractString)::Bool
    g = uppercase(strip(glass_name))
    return isempty(g) || g == "AIR" || g == "0"
end

"""
    is_mirror_glass(glass_name::AbstractString)::Bool

Returns true if the glass specification indicates a reflective mirror.
"""
function is_mirror_glass(glass_name::AbstractString)::Bool
    return uppercase(strip(glass_name)) == "MIRROR"
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

    # Determine starting optical surface (skip surface 0 / OBJ)
    start_idx = 1
    if surfaces[1].index == 0
        start_idx = 2
    end

    # Calculate cumulative axial positions along +Y for all surfaces
    # pos_y[k] is the vertex coordinate of surface k
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
            diam = max(2.0 * surf.semi_diameter, 0.005)
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
            diam = max(2.0 * surf.semi_diameter, 0.005)
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
            # Count consecutive glass segments
            k = i
            while k < n_surfs && !is_air(surfaces[k].glass_name) && !surfaces[k].is_mirror && !is_mirror_glass(surfaces[k].glass_name)
                k += 1
            end
            num_glass = k - i

            if num_glass == 1
                # Singlet lens: surfaces[i] (front) and surfaces[i+1] (back)
                s_front = surfaces[i]
                s_back = surfaces[i + 1]
                cthick = s_front.thickness
                diam = max(2.0 * s_front.semi_diameter, 2.0 * s_back.semi_diameter, 0.005)
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
                i = i + 1  # Next iteration will examine s_back (which has air, so it will advance to i+2)
            elseif num_glass == 2
                # Cemented Doublet: surfaces i, i+1, i+2
                s1 = surfaces[i]
                s2 = surfaces[i + 1]
                s3 = surfaces[i + 2]
                diam = max(2.0 * s1.semi_diameter, 2.0 * s2.semi_diameter, 2.0 * s3.semi_diameter, 0.005)
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
                # Cemented Triplet: surfaces i, i+1, i+2, i+3
                s1 = surfaces[i]
                s2 = surfaces[i + 1]
                s3 = surfaces[i + 2]
                s4 = surfaces[i + 3]
                diam = max(2.0 * s1.semi_diameter, 2.0 * s2.semi_diameter, 2.0 * s3.semi_diameter, 2.0 * s4.semi_diameter, 0.005)
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
                # More than 3 cemented surfaces: treat first 2 as doublet and continue
                s1 = surfaces[i]
                s2 = surfaces[i + 1]
                s3 = surfaces[i + 2]
                diam = max(2.0 * s1.semi_diameter, 2.0 * s2.semi_diameter, 2.0 * s3.semi_diameter, 0.005)
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

        # 4. Air gap or coordinate break / dummy surface
        i += 1
    end

    # 5. Image surface / Detector
    if n_surfs >= start_idx
        img_surf = surfaces[end]
        det_diam = max(2.0 * img_surf.semi_diameter, 0.01)
        push!(elements, ZmxDetector(
            name = "Detector",
            surface = img_surf,
            diameter = det_diam,
            axial_position = pos_y[end]
        ))
    end

    return elements
end

