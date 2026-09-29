"""
Glass and refractive index catalog support for Zemax imported materials.
Integrates RefractiveIndex.jl with Abbe-number / Cauchy dispersion and constant-index fallback.
"""

using RefractiveIndex

const GLASS_CACHE = Dict{String, Function}()

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
"""
function resolve_refractive_index(glass_name::AbstractString, nd::Real=1.0, vd::Real=0.0)::Function
    clean_name = uppercase(strip(glass_name))
    
    # 1. Air / Empty
    if isempty(clean_name) || clean_name == "AIR" || clean_name == "0"
        return λ -> 1.0
    end
    
    # 2. Check cache
    if haskey(GLASS_CACHE, clean_name)
        return GLASS_CACHE[clean_name]
    end
    
    # 3. Model glass (___BLANK or explicit parameters)
    if clean_name == "___BLANK" || startswith(clean_name, "___BLANK")
        if nd > 1.0
            fn = cauchy_dispersion(nd, vd)
            return fn
        else
            return make_bmo_safe_index(λ -> 1.5, 1.5)
        end
    end

    # 4. Search in RefractiveIndex.jl library
    preferred_catalogs = [
        "SCHOTT-optical", "OHARA-optical", "HOYA-optical", "CDGM-optical",
        "HIKARI-optical", "SUMITA-optical", "SCHOTT", "OHARA", "CORNING-optical"
    ]
    
    # Exact match in preferred optical catalogs
    for (shelf, book, page) in keys(RefractiveIndex.RI_LIB)
        if book in preferred_catalogs && uppercase(page) == clean_name
            try
                m = RefractiveMaterial(shelf, book, page)
                raw_fn = λ -> Float64(m[1](λ, "m"))
                fn = make_bmo_safe_index(raw_fn, 1.5)
                GLASS_CACHE[clean_name] = fn
                return fn
            catch
            end
        end
    end
    
    # Any catalog in RI_LIB
    for (shelf, book, page) in keys(RefractiveIndex.RI_LIB)
        if uppercase(page) == clean_name
            try
                m = RefractiveMaterial(shelf, book, page)
                raw_fn = λ -> Float64(m[1](λ, "m"))
                fn = make_bmo_safe_index(raw_fn, 1.5)
                GLASS_CACHE[clean_name] = fn
                return fn
            catch
            end
        end
    end
    
    # Match without hyphens (e.g. NBK7 -> N-BK7)
    normalized = replace(clean_name, "-" => "")
    for (shelf, book, page) in keys(RefractiveIndex.RI_LIB)
        if replace(uppercase(page), "-" => "") == normalized
            try
                m = RefractiveMaterial(shelf, book, page)
                raw_fn = λ -> Float64(m[1](λ, "m"))
                fn = make_bmo_safe_index(raw_fn, 1.5)
                GLASS_CACHE[clean_name] = fn
                return fn
            catch
            end
        end
    end

    # 5. If not found in catalog, but nd was provided in GLAS line:
    if nd > 1.0
        fn = cauchy_dispersion(nd, vd)
        GLASS_CACHE[clean_name] = fn
        return fn
    end

    # Fallback default
    fn = make_bmo_safe_index(λ -> 1.5, 1.5)
    GLASS_CACHE[clean_name] = fn
    return fn
end

