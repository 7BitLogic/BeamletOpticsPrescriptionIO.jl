"""
Zemax .zmx file parser.
Supports UTF-16LE, UTF-16BE, UTF-8, and ISO-8859-1 encodings.
"""

using StringEncodings

"""
    read_zmx_string(filepath::AbstractString)::String

Reads a Zemax .zmx file and decodes it according to its byte-order mark (BOM) or fallback encodings.
"""
function read_zmx_string(filepath::AbstractString)::String
    raw = read(filepath)
    if length(raw) >= 2 && raw[1] == 0xff && raw[2] == 0xfe
        return decode(raw, "UTF-16LE")
    elseif length(raw) >= 2 && raw[1] == 0xfe && raw[2] == 0xff
        return decode(raw, "UTF-16BE")
    elseif length(raw) >= 3 && raw[1] == 0xef && raw[2] == 0xbb && raw[3] == 0xbf
        return decode(raw, "UTF-8")
    end

    # Try UTF-8 first
    try
        return String(raw)
    catch
    end

    # Try UTF-16LE (many Zemax files omit BOM but are UTF-16LE)
    try
        return decode(raw, "UTF-16LE")
    catch
    end

    # Fallback to ISO-8859-1
    return decode(raw, "ISO-8859-1")
end

"""
    parse_zmx(filepath::AbstractString)::ZmxSystem

Parses a .zmx file from the filesystem into a `ZmxSystem`.
"""
function parse_zmx(filepath::AbstractString; configuration::Int=1)::ZmxSystem
    content = read_zmx_string(filepath)
    return parse_zmx_content(content; configuration=configuration)
end

"""
    parse_zmx_content(content::AbstractString; configuration::Int=1)::ZmxSystem

Parses the decoded text content of a .zmx file into a `ZmxSystem`.
Applies multi-configuration overrides (`THIC`) for the specified `configuration` (default: 1).
"""
function parse_zmx_content(content::AbstractString; configuration::Int=1)::ZmxSystem
    title = ""
    unit_sym = :mm
    scale_to_m = 1e-3
    version = ""
    enpd_val = 0.0
    obna_val = 0.0
    
    wavelengths = Float64[]
    spectral_weights = Float64[]
    primary_wvl_idx = 1
    
    surfaces = ZmxSurface[]
    thic_map = Dict{Tuple{Int, Int}, Float64}() # (surface_idx, conf_idx) => thickness_m
    raw_fields = Tuple{Float64, Float64, Float64}[]
    field_type = :angle
    expected_num_fields = 0
    
    # Temporary parsing state for current surface
    cur_surf_idx = -1
    cur_type = :STANDARD
    cur_comment = ""
    cur_curv = 0.0
    cur_disz = 0.0
    cur_glas = ""
    cur_nd = 1.0
    cur_vd = 0.0
    cur_diam = 0.0
    cur_coni = 0.0
    cur_parms = Float64[]
    cur_is_stop = false
    cur_dec_x = 0.0
    cur_dec_y = 0.0
    cur_tilt_x = 0.0
    cur_tilt_y = 0.0
    cur_tilt_z = 0.0

    function finish_current_surface!()
        if cur_surf_idx >= 0
            curv_si = cur_curv / scale_to_m
            rad_si = iszero(cur_curv) ? Inf : (1.0 / curv_si)
            
            disz_si = isinf(cur_disz) ? Inf : (cur_disz * scale_to_m)
            semi_diam_si = cur_diam * scale_to_m
            
            dec_x = cur_dec_x
            dec_y = cur_dec_y
            tilt_x = cur_tilt_x
            tilt_y = cur_tilt_y
            tilt_z = cur_tilt_z
            if cur_type == :COORDBRK && length(cur_parms) >= 5
                dec_x = cur_parms[1] * scale_to_m
                dec_y = cur_parms[2] * scale_to_m
                tilt_x = deg2rad(cur_parms[3])
                tilt_y = deg2rad(cur_parms[4])
                tilt_z = deg2rad(cur_parms[5])
            end

            scaled_parms = copy(cur_parms)
            if cur_type == :EVENASPH
                for i in eachindex(scaled_parms)
                    exp_pow = 2 * i - 1
                    scaled_parms[i] = scaled_parms[i] / (scale_to_m^exp_pow)
                end
            end

            is_mirror = (uppercase(strip(cur_glas)) == "MIRROR")

            surf = ZmxSurface(
                index = cur_surf_idx,
                surface_type = cur_type,
                comment = cur_comment,
                curvature = curv_si,
                radius = rad_si,
                thickness = disz_si,
                glass_name = cur_glas,
                nd = cur_nd,
                vd = cur_vd,
                semi_diameter = semi_diam_si,
                conic = cur_coni,
                parms = scaled_parms,
                is_stop = cur_is_stop,
                is_mirror = is_mirror,
                decenter_x = dec_x,
                decenter_y = dec_y,
                tilt_x = tilt_x,
                tilt_y = tilt_y,
                tilt_z = tilt_z
            )
            push!(surfaces, surf)
        end
    end

    for line in split(content, '\n')
        line = strip(line)
        isempty(line) && continue
        
        parts = split(line)
        isempty(parts) && continue
        cmd = uppercase(parts[1])
        
        if cmd == "VERS"
            version = length(parts) >= 2 ? parts[2] : ""
        elseif cmd == "NAME"
            title_match = match(r"NAME\s+(.*)", line)
            if title_match !== nothing
                title = strip(title_match.captures[1], [' ', '\t', '\"'])
            end
        elseif cmd == "UNIT"
            if length(parts) >= 2
                u = uppercase(parts[2])
                if u == "MM"
                    unit_sym = :mm
                    scale_to_m = 1e-3
                elseif u == "IN" || u == "INCH"
                    unit_sym = :inch
                    scale_to_m = 25.4e-3
                elseif u == "CM"
                    unit_sym = :cm
                    scale_to_m = 1e-2
                elseif u == "M" || u == "METER"
                    unit_sym = :meter
                    scale_to_m = 1.0
                end
            end
        elseif cmd == "ENPD" || cmd == "EPD"
            if length(parts) >= 2
                v = tryparse(Float64, parts[2])
                if v !== nothing
                    enpd_val = v * scale_to_m
                end
            end
        elseif cmd == "OBNA"
            if length(parts) >= 2
                v = tryparse(Float64, parts[2])
                if v !== nothing
                    obna_val = v
                end
            end
        elseif cmd == "FTYP"
            if length(parts) >= 2
                ft = tryparse(Int, parts[2])
                if ft !== nothing
                    if ft == 0
                        field_type = :angle
                    elseif ft == 1
                        field_type = :object_height
                    elseif ft == 2
                        field_type = :paraxial_image_height
                    elseif ft == 3
                        field_type = :real_image_height
                    else
                        field_type = :angle
                    end
                end
            end
            if length(parts) >= 4
                nf = tryparse(Int, parts[4])
                if nf !== nothing && nf >= 1
                    expected_num_fields = nf
                end
            end
        elseif cmd == "NFLD"
            if length(parts) >= 2
                nf = tryparse(Int, parts[2])
                if nf !== nothing && nf >= 1
                    expected_num_fields = nf
                end
            end
        elseif cmd == "XFLN" || cmd == "XFLM"
            vals = [tryparse(Float64, p) for p in parts[2:end]]
            for (i, v) in enumerate(vals)
                v === nothing && continue
                while length(raw_fields) < i
                    push!(raw_fields, (0.0, 0.0, 1.0))
                end
                raw_fields[i] = (v, raw_fields[i][2], raw_fields[i][3])
            end
        elseif cmd == "YFLN" || cmd == "YFLM"
            vals = [tryparse(Float64, p) for p in parts[2:end]]
            for (i, v) in enumerate(vals)
                v === nothing && continue
                while length(raw_fields) < i
                    push!(raw_fields, (0.0, 0.0, 1.0))
                end
                raw_fields[i] = (raw_fields[i][1], v, raw_fields[i][3])
            end
        elseif cmd == "FWGN" || cmd == "FWGM"
            vals = [tryparse(Float64, p) for p in parts[2:end]]
            for (i, v) in enumerate(vals)
                v === nothing && continue
                while length(raw_fields) < i
                    push!(raw_fields, (0.0, 0.0, 1.0))
                end
                raw_fields[i] = (raw_fields[i][1], raw_fields[i][2], v)
            end
        elseif cmd == "THIC"
            if length(parts) >= 4
                s_idx = tryparse(Int, parts[2])
                c_idx = tryparse(Int, parts[3])
                val = tryparse(Float64, parts[4])
                if s_idx !== nothing && c_idx !== nothing && val !== nothing
                    thic_map[(s_idx, c_idx)] = val * scale_to_m
                end
            end
        elseif cmd == "WAVM" || cmd == "WAVS"
            if length(parts) >= 3
                idx = tryparse(Int, parts[2])
                wvl_um = tryparse(Float64, parts[3])
                wt = length(parts) >= 4 ? tryparse(Float64, parts[4]) : 1.0
                if idx !== nothing && wvl_um !== nothing && idx >= 1
                    wvl_m = wvl_um * 1e-6
                    while length(wavelengths) < idx
                        push!(wavelengths, 550e-9)
                        push!(spectral_weights, 1.0)
                    end
                    wavelengths[idx] = wvl_m
                    spectral_weights[idx] = wt !== nothing ? wt : 1.0
                end
            end
        elseif cmd == "WAVL"
            empty!(wavelengths)
            for p in parts[2:end]
                v = tryparse(Float64, p)
                if v !== nothing
                    push!(wavelengths, v * 1e-6)
                end
            end
        elseif cmd == "WWGT"
            empty!(spectral_weights)
            for p in parts[2:end]
                v = tryparse(Float64, p)
                if v !== nothing
                    push!(spectral_weights, v)
                end
            end
        elseif cmd == "PWAV"
            if length(parts) >= 2
                p_idx = tryparse(Int, parts[2])
                if p_idx !== nothing && p_idx >= 1
                    primary_wvl_idx = p_idx
                end
            end
        elseif cmd == "SURF"
            finish_current_surface!()
            
            cur_surf_idx = length(parts) >= 2 ? (tryparse(Int, parts[2]) === nothing ? cur_surf_idx + 1 : tryparse(Int, parts[2])) : (cur_surf_idx + 1)
            cur_type = :STANDARD
            cur_comment = ""
            cur_curv = 0.0
            cur_disz = 0.0
            cur_glas = ""
            cur_nd = 1.0
            cur_vd = 0.0
            cur_diam = 0.0
            cur_coni = 0.0
            empty!(cur_parms)
            cur_is_stop = false
            cur_dec_x = 0.0
            cur_dec_y = 0.0
            cur_tilt_x = 0.0
            cur_tilt_y = 0.0
            cur_tilt_z = 0.0
        elseif cur_surf_idx >= 0
            if cmd == "TYPE"
                if length(parts) >= 2
                    cur_type = Symbol(uppercase(parts[2]))
                end
            elseif cmd == "COMM"
                comm_match = match(r"COMM\s+(.*)", line)
                if comm_match !== nothing
                    cur_comment = strip(comm_match.captures[1], [' ', '\t', '\"'])
                end
            elseif cmd == "CURV"
                if length(parts) >= 2
                    v = tryparse(Float64, parts[2])
                    cur_curv = v !== nothing ? v : 0.0
                end
            elseif cmd == "DISZ"
                if length(parts) >= 2
                    p_up = uppercase(parts[2])
                    if p_up == "INFINITY" || p_up == "INF"
                        cur_disz = Inf
                    else
                        v = tryparse(Float64, parts[2])
                        cur_disz = v !== nothing ? v : 0.0
                    end
                end
            elseif cmd == "GLAS"
                if length(parts) >= 2
                    raw_glas = parts[2]
                    if startswith(raw_glas, "___BLANK")
                        cur_glas = "___BLANK"
                    elseif raw_glas != "AIR" && raw_glas != "0"
                        cur_glas = raw_glas
                    else
                        cur_glas = "AIR"
                    end
                    
                    # Zemax GLAS format: GLAS <Name> <SolveType> <SolveParm> <Nd> <Vd> <DpgF> ...
                    # If parts has at least 5 elements, parts[5] is Nd
                    if length(parts) >= 5
                        val_nd = tryparse(Float64, parts[5])
                        if val_nd !== nothing && val_nd > 0.0
                            cur_nd = val_nd
                        end
                    end
                    if length(parts) >= 6
                        val_vd = tryparse(Float64, parts[6])
                        if val_vd !== nothing && val_vd > 0.0
                            cur_vd = val_vd
                        end
                    end
                    
                    # Fallback to scanning numbers if Nd is still <= 1.0
                    if cur_nd <= 1.0
                        nums = [tryparse(Float64, p) for p in parts[3:end]]
                        filter!(x -> x !== nothing, nums)
                        if length(nums) >= 4 && nums[3] > 1.0
                            cur_nd = nums[3]
                            cur_vd = length(nums) >= 4 ? nums[4] : 0.0
                        elseif length(nums) >= 2 && nums[1] > 1.0
                            cur_nd = nums[1]
                            cur_vd = nums[2]
                        end
                    end
                end
            elseif cmd == "DIAM" || cmd == "SDMA"
                if length(parts) >= 2
                    v = tryparse(Float64, parts[2])
                    if v !== nothing
                        cur_diam = v
                    end
                end
            elseif cmd == "CONI"
                if length(parts) >= 2
                    v = tryparse(Float64, parts[2])
                    cur_coni = v !== nothing ? v : 0.0
                end
            elseif cmd == "PARM"
                if length(parts) >= 3
                    pidx = tryparse(Int, parts[2])
                    pval = tryparse(Float64, parts[3])
                    if pidx !== nothing && pval !== nothing && pidx >= 1
                        while length(cur_parms) < pidx
                            push!(cur_parms, 0.0)
                        end
                        cur_parms[pidx] = pval
                    end
                end
            elseif cmd == "STOP"
                cur_is_stop = true
            end
        end
    end
    
    finish_current_surface!()

    # Apply multi-configuration thickness overrides for the requested configuration
    for i in 1:length(surfaces)
        s = surfaces[i]
        if haskey(thic_map, (s.index, configuration))
            new_thic = thic_map[(s.index, configuration)]
            surfaces[i] = ZmxSurface(
                index = s.index,
                surface_type = s.surface_type,
                comment = s.comment,
                curvature = s.curvature,
                radius = s.radius,
                thickness = new_thic,
                glass_name = s.glass_name,
                nd = s.nd,
                vd = s.vd,
                semi_diameter = s.semi_diameter,
                conic = s.conic,
                parms = s.parms,
                is_stop = s.is_stop,
                is_mirror = s.is_mirror,
                decenter_x = s.decenter_x,
                decenter_y = s.decenter_y,
                tilt_x = s.tilt_x,
                tilt_y = s.tilt_y,
                tilt_z = s.tilt_z
            )
        end
    end

    while length(wavelengths) > 1 && isapprox(wavelengths[end], 550e-9, atol=1e-12)
        pop!(wavelengths)
        if !isempty(spectral_weights)
            pop!(spectral_weights)
        end
    end

    if isempty(wavelengths)
        wavelengths = [587.56e-9]
        spectral_weights = [1.0]
    end

    if expected_num_fields > 0 && length(raw_fields) >= expected_num_fields
        raw_fields = raw_fields[1:expected_num_fields]
    end

    if isempty(raw_fields)
        push!(raw_fields, (0.0, 0.0, 1.0))
    end

    # Scale field coordinates from file units to meters if they represent physical heights
    if field_type in (:object_height, :paraxial_image_height, :real_image_height)
        raw_fields = [(f[1] * scale_to_m, f[2] * scale_to_m, f[3]) for f in raw_fields]
    end

    obj_dist = Inf
    if !isempty(surfaces) && surfaces[1].index == 0
        obj_dist = surfaces[1].thickness
    end

    return ZmxSystem(
        title = title,
        unit = unit_sym,
        scale_to_m = scale_to_m,
        version = version,
        wavelengths = wavelengths,
        spectral_weights = spectral_weights,
        primary_wavelength_idx = clamp(primary_wvl_idx, 1, length(wavelengths)),
        enpd = enpd_val,
        obna = obna_val,
        object_distance = obj_dist,
        fields = raw_fields,
        field_type = field_type,
        surfaces = surfaces
    )
end

