"""
Glass and refractive index catalog support for Zemax imported materials.
Integrates RefractiveIndex.jl with Abbe-number / Cauchy dispersion and constant-index fallback.
"""

using RefractiveIndex

# In-memory lookup cache
const GLASS_CACHE = Dict{String, Function}()

# Wavelength constants for Fraunhofer lines [m]
const λ_d = 587.5618e-9  # Helium d-line (yellow)
const λ_F = 486.1327e-9  # Hydrogen F-line (blue)
const λ_C = 656.2725e-9  # Hydrogen C-line (red)

"""
    cauchy_dispersion(nd::Real, vd::Real)::Function

Constructs a Cauchy-approximated dispersion function `λ -> n(λ)` from the refractive index
at the d-line `nd` and Abbe number `vd`.
Formula: n(λ) = A + B / λ^2
"""
function cauchy_dispersion(nd::Real, vd::Real)::Function
    if vd <= 0.0 || isinf(vd) || isnan(vd)
        return λ -> Float64(nd)
    end
    # n_F - n_C = (nd - 1) / vd
    Δn = (nd - 1.0) / vd
    inv_λF2_minus_inv_λC2 = (1.0 / λ_F^2) - (1.0 / λ_C^2)
    B = Δn / inv_λF2_minus_inv_λC2
    A = nd - (B / λ_d^2)
    
    return λ -> Float64(A + B / (λ^2))
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
    
    # 3. Model glass (___BLANK or explicit parameters without known catalog name)
    if clean_name == "___BLANK" || startswith(clean_name, "___BLANK")
        if nd > 1.0
            fn = cauchy_dispersion(nd, vd)
            return fn
        else
            return λ -> 1.5
        end
    end

    # 4. Search in RefractiveIndex.jl library
    preferred_catalogs = [
        "SCHOTT-optical", "OHARA-optical", "HOYA-optical", "CDGM-optical",
        "HIKARI-optical", "SUMITA-optical", "SCHOTT", "OHARA", "CORNING-optical"
    ]
    
    # First try exact match in preferred optical catalogs
    for (shelf, book, page) in keys(RefractiveIndex.RI_LIB)
        if book in preferred_catalogs && uppercase(page) == clean_name
            try
                m = RefractiveMaterial(shelf, book, page)
                fn = λ -> Float64(m[1](λ, "m"))
                GLASS_CACHE[clean_name] = fn
                return fn
            catch
            end
        end
    end
    
    # Try any catalog in RI_LIB
    for (shelf, book, page) in keys(RefractiveIndex.RI_LIB)
        if uppercase(page) == clean_name
            try
                m = RefractiveMaterial(shelf, book, page)
                fn = λ -> Float64(m[1](λ, "m"))
                GLASS_CACHE[clean_name] = fn
                return fn
            catch
            end
        end
    end
    
    # Also try without hyphens (e.g. NBK7 -> N-BK7 or N-BK7 -> NBK7)
    normalized = replace(clean_name, "-" => "")
    for (shelf, book, page) in keys(RefractiveIndex.RI_LIB)
        if replace(uppercase(page), "-" => "") == normalized
            try
                m = RefractiveMaterial(shelf, book, page)
                fn = λ -> Float64(m[1](λ, "m"))
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
    @warn "Glass '$glass_name' not found in database; using standard n = 1.5"
    fn = λ -> 1.5
    GLASS_CACHE[clean_name] = fn
    return fn
end

