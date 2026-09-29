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
function parse_zmx(filepath::AbstractString)::ZmxSystem
    content = read_zmx_string(filepath)
    return parse_zmx_content(content)
end

"""
    parse_zmx_content(content::AbstractString)::ZmxSystem

Parses the decoded text content of a .zmx file into a `ZmxSystem`.
"""
function parse_zmx_content(content::AbstractString)::ZmxSystem
    title = ""
    unit_sym = :mm
    scale_to_m = 1e-3
    version = ""
    
    wavelengths = Float64[]
    spectral_weights = Float64[]
    primary_wvl_idx = 1
    
    surfaces = ZmxSurface[]
    
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
        elseif cmd == "WAVM"
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
                        nums = [tryparse(Float64, p) for p in parts[3:end]]
                        filter!(x -> x !== nothing, nums)
                        if length(nums) >= 4
                            cur_nd = nums[3]
                            cur_vd = nums[4]
                        elseif length(nums) >= 2
                            cur_nd = nums[1]
                            cur_vd = nums[2]
                        end
                    elseif raw_glas != "AIR" && raw_glas != "0"
                        cur_glas = raw_glas
                        nums = [tryparse(Float64, p) for p in parts[3:end]]
                        filter!(x -> x !== nothing, nums)
                        if length(nums) >= 2
                            cur_nd = nums[end-1]
                            cur_vd = nums[end]
                        end
                    else
                        cur_glas = "AIR"
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

    return ZmxSystem(
        title = title,
        unit = unit_sym,
        scale_to_m = scale_to_m,
        version = version,
        wavelengths = wavelengths,
        spectral_weights = spectral_weights,
        primary_wavelength_idx = clamp(primary_wvl_idx, 1, length(wavelengths)),
        surfaces = surfaces
    )
end

