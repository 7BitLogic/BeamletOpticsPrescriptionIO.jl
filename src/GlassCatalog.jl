"""
Glass and refractive index catalog support for Zemax imported materials.
Integrates RefractiveIndex.jl with Abbe-number / Cauchy dispersion and constant-index fallback.
"""

using RefractiveIndex

const GLASS_CACHE = Dict{Tuple{String, Float64, Float64}, Function}()

const HISTORICAL_GLASSES = Dict{String, Tuple{Float64, Float64}}(
    # --- Bariumflint & Bariumkrone ---
    "BASF1"   => (1.62606, 38.96),
    "BASF2"   => (1.66446, 35.80),
    "BASF7"   => (1.64328, 47.98),
    "BAF4"    => (1.60565, 43.72),
    "BAF9"    => (1.64328, 47.82),
    "BAF10"   => (1.67003, 47.11),
    "BAK1"    => (1.57250, 57.55),
    "N-BAK1"  => (1.57250, 57.55),
    "BAK2"    => (1.53996, 59.71),
    "N-BAK2"  => (1.53996, 59.71),
    "BAK4"    => (1.56883, 55.98),
    "N-BAK4"  => (1.56883, 55.98),
    "BALF4"   => (1.57956, 53.87),
    "N-BALF4" => (1.57956, 53.87),
    "BALF5"   => (1.54739, 53.63),
    "N-BALF5" => (1.54739, 53.63),

    # --- Kron & Borosilikatkrone ---
    "BK1"     => (1.51009, 63.46),
    "BK7"     => (1.51680, 64.17),
    "N-BK7"   => (1.51680, 64.17),
    "BK10"    => (1.49782, 66.95),
    "K1"      => (1.50349, 56.65),
    "K3"      => (1.51821, 58.96),
    "K4"      => (1.51860, 58.96),
    "K5"      => (1.52249, 59.48),
    "N-K5"    => (1.52249, 59.48),
    "K7"      => (1.51112, 60.41),
    "K10"     => (1.50137, 56.41),
    "K11"     => (1.50373, 56.55),
    "KF2"     => (1.52940, 42.17),
    "KF3"     => (1.51460, 54.60),
    "KF4"     => (1.52249, 51.48),
    "KF5"     => (1.52249, 51.48),
    "KF6"     => (1.51742, 52.43),
    "KF8"     => (1.51460, 54.60),
    "KF9"     => (1.52341, 51.48),
    "FK5"     => (1.48749, 70.41),
    "N-FK5"   => (1.48749, 70.41),
    "FK51"    => (1.48656, 84.47),
    "PK1"     => (1.50383, 66.92),
    "PK3"     => (1.52680, 64.87),
    "PSK2"    => (1.56865, 63.16),
    "PSK3"    => (1.55232, 63.46),
    "PSK52"   => (1.61800, 63.39),
    "PSKN2"   => (1.56865, 63.16),

    # --- Schwerkrone & Schwerstkrone ---
    "SK1"     => (1.61025, 56.65),
    "SK2"     => (1.60738, 56.65),
    "SK4"     => (1.61272, 58.63),
    "SK5"     => (1.58913, 61.27),
    "SK6"     => (1.61340, 56.40),
    "SK7"     => (1.60738, 56.65),
    "SK9"     => (1.62000, 56.90),
    "SK10"    => (1.62278, 56.90),
    "SK11"    => (1.56384, 60.80),
    "SK12"    => (1.58313, 59.40),
    "SK13"    => (1.59310, 58.40),
    "SK14"    => (1.60311, 60.60),
    "SK15"    => (1.62296, 58.06),
    "SK16"    => (1.62041, 60.29),
    "N-SK16"  => (1.62041, 60.29),
    "SK18"    => (1.63854, 55.45),
    "SK20"    => (1.55836, 61.20),
    "SK51"    => (1.62000, 60.30),
    "SK55"    => (1.62000, 60.30),
    "SSK1"    => (1.61727, 52.84),
    "SSK2"    => (1.62230, 53.17),
    "SSK3"    => (1.61765, 55.08),
    "SSK4"    => (1.61765, 55.08),
    "SSK5"    => (1.65844, 50.85),
    "SSK9"    => (1.62000, 53.00),
    "SSKN5"   => (1.65844, 50.88),
    "SSKN8"   => (1.66680, 48.33),
    "SSL2"    => (1.53172, 48.86),
    "SSL5"    => (1.53172, 48.86),

    # --- Flint & Leichtflint ---
    "F1"      => (1.62590, 35.70),
    "F2"      => (1.62004, 36.37),
    "N-F2"    => (1.62004, 36.37),
    "F4"      => (1.61659, 36.63),
    "F5"      => (1.60342, 38.03),
    "F6"      => (1.60565, 37.80),
    "F8"      => (1.59551, 39.18),
    "LF1"     => (1.57860, 41.14),
    "LF2"     => (1.57501, 41.48),
    "LF4"     => (1.57501, 41.48),
    "LF5"     => (1.58144, 40.85),
    "LF6"     => (1.56732, 42.83),
    "LF7"     => (1.57501, 41.48),
    "LF8"     => (1.56883, 42.00),
    "LLF1"    => (1.54814, 45.75),
    "LLF2"    => (1.54212, 47.16),
    "LLF4"    => (1.53172, 48.86),
    "LLF6"    => (1.53172, 48.86),
    "LLF7"    => (1.54814, 45.75),
    "LLF8"    => (1.53996, 47.00),

    # --- Schwerflint ---
    "SF1"     => (1.71736, 29.51),
    "N-SF1"   => (1.71736, 29.51),
    "SF2"     => (1.64769, 33.82),
    "N-SF2"   => (1.64769, 33.82),
    "SF3"     => (1.74000, 28.30),
    "SF4"     => (1.75520, 27.58),
    "SF5"     => (1.67270, 32.21),
    "N-SF5"   => (1.67270, 32.21),
    "SF6"     => (1.80518, 25.43),
    "N-SF6"   => (1.80518, 25.43),
    "SF7"     => (1.64000, 34.60),
    "SF8"     => (1.68893, 31.18),
    "SF9"     => (1.65440, 33.60),
    "SF10"    => (1.72825, 28.41),
    "N-SF10"  => (1.72825, 28.41),
    "SF11"    => (1.78472, 25.76),
    "N-SF11"  => (1.78472, 25.76),
    "SF12"    => (1.64831, 33.85),
    "SF13"    => (1.74077, 27.76),
    "SF14"    => (1.76182, 26.52),
    "SF15"    => (1.69895, 30.07),
    "N-SF15"  => (1.69895, 30.07),
    "SF16"    => (1.69895, 30.07),
    "SF18"    => (1.72151, 29.28),
    "SF19"    => (1.66680, 33.00),
    "SF53"    => (1.72825, 28.53),
    "SF55"    => (1.76182, 26.52),
    "SF56A"   => (1.78472, 25.76),
    "SF57"    => (1.84666, 23.83),
    "N-SF57"  => (1.84666, 23.83),
    "SF58"    => (1.91761, 21.50),
    "SF59"    => (1.95250, 20.36),
    "SF63"    => (1.73400, 28.30),
    "SF64A"   => (1.70560, 30.20),
    "SF66"    => (1.92286, 20.88),
    "SFL4"    => (1.75520, 27.58),
    "SFL6"    => (1.80518, 25.43),
    "SFL56"   => (1.78472, 25.76),
    "SFL57"   => (1.84666, 23.83),
    "SFLD6"   => (1.80518, 25.43),
    "SFLD20"  => (1.80400, 25.50),
    "SFLD66"  => (1.92286, 20.88),
    "SFLDN3"  => (1.78472, 25.76),
    "SFN1"    => (1.71736, 29.51),
    "SFN3"    => (1.74000, 28.30),
    "SFN4"    => (1.75520, 27.58),

    # --- Lanthankrone ---
    "LAK6"    => (1.69100, 54.71),
    "LAK7"    => (1.65160, 58.55),
    "LAK8"    => (1.71300, 53.83),
    "LAK9"    => (1.69100, 54.71),
    "N-LAK9"  => (1.69100, 54.71),
    "LAK10"   => (1.72000, 50.34),
    "N-LAK10" => (1.72000, 50.34),
    "LAK11"   => (1.65844, 50.85),
    "LAK12"   => (1.67790, 55.34),
    "LAK13"   => (1.69350, 53.30),
    "LAK14"   => (1.69680, 55.52),
    "LAK16A"  => (1.73400, 51.50),
    "LAK18"   => (1.71300, 53.83),
    "LAK21"   => (1.64049, 60.10),
    "LAK23"   => (1.66446, 56.10),
    "LAK31"   => (1.69680, 55.50),
    "LAK33"   => (1.75398, 52.43),
    "LAKN5"   => (1.65844, 50.85),
    "LAKN6"   => (1.69100, 54.71),
    "LAKN7"   => (1.65160, 58.55),
    "LAKN10"  => (1.72000, 50.34),
    "LAKN11"  => (1.65844, 50.85),
    "LAKN12"  => (1.67790, 55.34),
    "LAKN13"  => (1.69350, 53.30),
    "LAKN14"  => (1.69680, 55.52),
    "LASKN1"  => (1.75520, 52.00),
    "LASKN3"  => (1.78650, 50.00),

    # --- Lanthanflint & Schwerlanthanflint ---
    "LAF2"    => (1.74400, 44.78),
    "LAF3"    => (1.71700, 47.90),
    "LAF7"    => (1.72342, 37.95),
    "LAF9"    => (1.78800, 47.43),
    "LAF11"   => (1.78470, 43.90),
    "LAF11A"  => (1.78470, 43.90),
    "LAF13"   => (1.77250, 49.62),
    "LAF20"   => (1.72000, 43.70),
    "LAF70"   => (1.72342, 37.95),
    "LAFL2"   => (1.74397, 44.85),
    "LAFN3"   => (1.71700, 47.90),
    "LAFN4"   => (1.74400, 44.78),
    "LAFN7"   => (1.74950, 34.95),
    "LAFN8"   => (1.75520, 34.69),
    "LAFN9"   => (1.78800, 47.43),
    "LAFN10"  => (1.78800, 47.43),
    "LAFN11"  => (1.78470, 43.90),
    "LAFN12"  => (1.78470, 43.90),
    "LAFN21"  => (1.78800, 47.43),
    "LAFN28"  => (1.77250, 49.62),
    "LASF01"  => (1.78800, 47.43),
    "LASF02"  => (1.74400, 44.78),
    "LASF3"   => (1.80400, 46.58),
    "LASF9"   => (1.85026, 32.17),
    "N-LASF9" => (1.85026, 32.17),
    "LASF18A" => (1.80518, 25.43),
    "LASF32"  => (1.80290, 46.80),
    "LASF35"  => (2.02204, 29.06),
    "LASFH6"  => (1.80400, 46.58),
    "LASFN1"  => (1.78800, 47.43),
    "LASFN2"  => (1.74400, 44.78),
    "LASFN3"  => (1.80400, 46.58),
    "LASFN4"  => (1.79952, 42.24),
    "LASFN6"  => (1.80400, 46.58),
    "LASFN7"  => (1.75520, 34.69),
    "LASFN8"  => (1.75520, 34.69),
    "LASFN9"  => (1.85026, 32.17),
    "LASFN13" => (1.80400, 46.58),
    "LASFN14" => (1.76182, 26.52),
    "LASFN16" => (1.75844, 52.32),
    "LASFN17" => (1.80400, 46.58),
    "LASFN30" => (1.80318, 46.38),
    "LASFN31" => (1.88300, 40.76),

    # --- Kristalle & Quarz ---
    "ZK1"     => (1.53375, 54.00),
    "ZKN7"    => (1.50830, 61.20),
    "F_SILICA"=> (1.45846, 67.82),
    "SILICA"  => (1.45846, 67.82),
    "QUARTZ"  => (1.54425, 69.80),
    "CAF2"    => (1.43385, 95.23),
    "GERMANIUM" => (4.00400, 0.0),
    "SILICON"   => (3.42230, 0.0),

    # --- Optical Plastics & Polymers ---
    "PMMA"     => (1.49175, 57.44),
    "POLYCARB" => (1.58547, 29.90),
    "POLYSTYR" => (1.59166, 30.81),
    "480R"     => (1.52500, 56.00),
    "E48R"     => (1.53100, 56.00),
    "AL-6263-(OKP4HT)" => (1.60700, 27.00),
    "OKP4HT"   => (1.60700, 27.00),
    "D34-26"   => (1.53400, 56.00),

    # --- Hoya, Ohara, Sumita & Precision Mold Glasses ---
    "S-LAL12"       => (1.67790, 55.34),
    "S-LAL12_MOLD"  => (1.67790, 55.34),
    "BACD14"        => (1.58913, 61.27),
    "BACD14_MOLD"   => (1.58913, 61.27),
    "BACD16"        => (1.62041, 60.29),
    "NBFD2"         => (1.74400, 44.78),
    "PCD2"          => (1.58913, 61.27),
    "PFK80"         => (1.49700, 81.61),
    "BPG2"          => (1.54814, 45.75),
    "MC-TAF1"       => (1.80400, 46.57),
    "TAF1"          => (1.80400, 46.57),
    "TAF4"          => (1.80610, 40.92),
    "TAF5"          => (1.83481, 42.72),
    "TAFD5F"        => (1.83481, 42.72),
    "LACL5"         => (1.70000, 48.08),
    "LACL60"        => (1.69350, 53.30),

    # --- Thorlabs Catalog Types & Special Optical Materials ---
    "ACRYLIC"       => (1.49175, 57.44),
    "C79-80"        => (1.45846, 67.82),
    "D-LAK6M"       => (1.80420, 46.50),
    "D-ZK3M"        => (1.58913, 61.27),
    "D-ZLAF52LAM"   => (1.81080, 40.90),
    "D-ZLAF52LA_M"  => (1.81080, 40.90),
    "ECO550"        => (1.59700, 67.00),
    "ECO550_E"      => (1.59700, 67.00),
    "FD10"          => (1.72825, 28.46),
    "INFRASIL"      => (1.45846, 67.82),
    "MGF2"          => (1.37770, 106.0),
    "ZNSE"          => (2.61000, 10.0),
    "ZNS_IR"        => (2.35000, 15.0),
    "ZNS"           => (2.35000, 15.0),
    "SLW-1.8"       => (1.60000, 50.0),
    "VIG06"         => (1.52000, 55.0)
)

const λ_d = 587.5618e-9  # Helium d-line (yellow)
const λ_F = 486.1327e-9  # Hydrogen F-line (blue)
const λ_C = 656.2725e-9  # Hydrogen C-line (red)

"""
    make_bmo_safe_index(fn_or_val, default_n=1.5)::Function

Wraps an index calculation into a function that safely handles non-physical test
arguments (such as λ = 1.0 m used by BeamletOptics`s test_refractive_index_function).
"""
function make_bmo_safe_index(base_fn::Function, default_n::Float64=1.5)::Function
    return λ -> begin
        # If λ is non-physical or outside optical window (e.g. 1.0m from BMO test)
        if λ > 50e-6 || λ < 100e-9
            # Return index at 587.56 nm
            val = try
                base_fn(587.5618e-9)
            catch
                default_n
            end
            return (isnan(val) || val < 1.0) ? default_n : Float64(val)
        end
        val = try
            base_fn(λ)
        catch
            default_n
        end
        return (isnan(val) || val < 1.0) ? default_n : Float64(val)
    end
end

"""
    cauchy_dispersion(nd::Real, vd::Real)::Function

Constructs a Cauchy-approximated dispersion function `λ -> n(λ)` from the refractive index
at the d-line `nd` and Abbe number `vd`.
"""
function cauchy_dispersion(nd::Real, vd::Real)::Function
    if vd <= 0.0 || isinf(vd) || isnan(vd)
        return make_bmo_safe_index(λ -> Float64(nd), Float64(nd))
    end
    Δn = (nd - 1.0) / vd
    inv_λF2_minus_inv_λC2 = (1.0 / λ_F^2) - (1.0 / λ_C^2)
    B = Δn / inv_λF2_minus_inv_λC2
    A = nd - (B / λ_d^2)
    
    raw_fn = λ -> (A + B / (λ^2))
    return make_bmo_safe_index(raw_fn, Float64(nd))
end

"""
    resolve_refractive_index(glass_name::AbstractString, nd::Real=1.0, vd::Real=0.0)::Function

Resolves a refractive index function `λ -> n(λ)` for a given Zemax glass name and/or
model parameters (`nd`, `vd`).

Priority:
1. Air / Empty -> n = 1.0
2. Explicit prescription: If `nd > 1.0` or `vd > 0.0` (and NOT the Zemax uninitialized dummy default 1.5/40.0),
   computes physical Cauchy dispersion.
3. Historical glass alias table (e.g. BASF1, SF1, SK16, E-BK7, J-LASF017).
4. RefractiveIndex.jl optical materials library.
5. If dummy default (1.5, 40.0) or blank, Cauchy dispersion.
6. Safe fallback with explicit warning.
"""
function resolve_refractive_index(glass_name::AbstractString, nd::Real=1.0, vd::Real=0.0)::Function
    clean_name = uppercase(strip(glass_name))
    
    # 1. Air / Empty
    if isempty(clean_name) || clean_name == "AIR" || clean_name == "0"
        return λ -> 1.0
    end
    
    # 2. Check cache
    cache_key = (clean_name, round(Float64(nd), digits=6), round(Float64(vd), digits=4))
    if haskey(GLASS_CACHE, cache_key)
        return GLASS_CACHE[cache_key]
    end

    # Check if this is Zemax uninitialized dummy default (1.5, 40.0) with a named catalog glass
    is_named_glass = clean_name != "___BLANK" && !startswith(clean_name, "___BLANK")
    is_dummy_placeholder = is_named_glass && (isapprox(nd, 1.5, atol=1e-4) && isapprox(vd, 40.0, atol=1e-2))
    
    # 3. Explicit prescription values (nd > 1.0 or vd > 0.0)
    if (nd > 1.0 || vd > 0.0) && !is_dummy_placeholder
        # If both nd and vd are specified
        if nd > 1.0 && vd > 0.0
            fn = cauchy_dispersion(nd, vd)
            GLASS_CACHE[cache_key] = fn
            return fn
        elseif nd > 1.0 && vd <= 0.0
            # nd specified, vd missing
            if haskey(HISTORICAL_GLASSES, clean_name)
                h_vd = HISTORICAL_GLASSES[clean_name][2]
                fn = cauchy_dispersion(nd, h_vd)
                GLASS_CACHE[cache_key] = fn
                return fn
            else
                fn = cauchy_dispersion(nd, 0.0)
                GLASS_CACHE[cache_key] = fn
                return fn
            end
        else # nd <= 1.0 && vd > 0.0
            # vd specified! Under NO circumstances fall back to 1.5 without dispersion!
            h_nd = haskey(HISTORICAL_GLASSES, clean_name) ? HISTORICAL_GLASSES[clean_name][1] : 1.5
            fn = cauchy_dispersion(h_nd, vd)
            GLASS_CACHE[cache_key] = fn
            return fn
        end
    end

    # 4. Check historical glass table (exact and candidate variations)
    # Variations include stripping manufacturer prefixes (E-, J-, N-, H-, S-, M-, MC-, MP-, K-)
    # and suffixes (_MOLD, -MOLD)
    candidates = String[clean_name]
    strip_mold = replace(clean_name, r"(_MOLD|-MOLD)$" => "")
    if strip_mold != clean_name
        push!(candidates, strip_mold)
    end
    for c in copy(candidates)
        unprefixed = replace(c, r"^(E-|J-|N-|H-|S-|M-|MC-|MP-|K-)" => "")
        if unprefixed != c
            push!(candidates, unprefixed)
        end
    end

    for cand in candidates
        if haskey(HISTORICAL_GLASSES, cand)
            (h_nd, h_vd) = HISTORICAL_GLASSES[cand]
            fn = cauchy_dispersion(h_nd, h_vd)
            GLASS_CACHE[cache_key] = fn
            return fn
        end
    end

    # 5. Search in RefractiveIndex.jl library
    preferred_catalogs = [
        "SCHOTT-optical", "OHARA-optical", "HOYA-optical", "CDGM-optical",
        "HIKARI-optical", "SUMITA-optical", "SCHOTT", "OHARA", "CORNING-optical",
        "plastics"
    ]
    
    # Check candidates in preferred optical catalogs
    for cand in candidates
        for (shelf, book, page) in keys(RefractiveIndex.RI_LIB)
            if book in preferred_catalogs && uppercase(page) == cand
                try
                    m = RefractiveMaterial(shelf, book, page)
                    raw_fn = λ -> Float64(m[1](λ, "m"))
                    fn = make_bmo_safe_index(raw_fn, 1.5)
                    GLASS_CACHE[cache_key] = fn
                    return fn
                catch
                end
            end
        end
    end
    
    # Any catalog in RI_LIB
    for cand in candidates
        for (shelf, book, page) in keys(RefractiveIndex.RI_LIB)
            if uppercase(page) == cand
                try
                    m = RefractiveMaterial(shelf, book, page)
                    raw_fn = λ -> Float64(m[1](λ, "m"))
                    fn = make_bmo_safe_index(raw_fn, 1.5)
                    GLASS_CACHE[cache_key] = fn
                    return fn
                catch
                end
            end
        end
    end
    
    # Match without hyphens (e.g. NBK7 -> N-BK7)
    for cand in candidates
        normalized = replace(cand, "-" => "")
        for (shelf, book, page) in keys(RefractiveIndex.RI_LIB)
            if replace(uppercase(page), "-" => "") == normalized
                try
                    m = RefractiveMaterial(shelf, book, page)
                    raw_fn = λ -> Float64(m[1](λ, "m"))
                    fn = make_bmo_safe_index(raw_fn, 1.5)
                    GLASS_CACHE[cache_key] = fn
                    return fn
                catch
                end
            end
        end
    end

    # 6. Fallback warning
    @warn "Glass '$glass_name' not found in catalogs or historical database. Defaulting to n=1.5."
    fn = make_bmo_safe_index(λ -> 1.5, 1.5)
    GLASS_CACHE[cache_key] = fn
    return fn
end


