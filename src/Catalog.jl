"""
Catalog integration for standard manufacturer Zemax catalogs (e.g., Thorlabs, Edmund Optics).
Provides high-level utilities to discover, query, and directly instantiate optical components
as native BeamletOptics.jl lens objects.
"""

using BeamletOptics
using LinearAlgebra

"""
Registered standard catalog directories.
"""
const DEFAULT_CATALOGS = Dict{Symbol, String}(
    :thorlabs => normpath(joinpath(@__DIR__, "..", "test_data", "thorlabs", "zmx"))
)

"""
    register_catalog!(name::Symbol, path::AbstractString)

Registers or updates the directory path for an optical component catalog `name`.
"""
function register_catalog!(name::Symbol, path::AbstractString)
    clean_path = normpath(abspath(path))
    DEFAULT_CATALOGS[name] = clean_path
    return clean_path
end

"""
    get_catalog_path(catalog::Union{Symbol, AbstractString})::String

Resolves the filesystem directory for `catalog`. Supports registered symbols (e.g. `:thorlabs`),
environment variable overrides (e.g. `THORLABS_CATALOG_PATH`), and direct directory paths.
"""
function get_catalog_path(catalog::Symbol)::String
    # 1. Registered catalogs
    if haskey(DEFAULT_CATALOGS, catalog)
        path = DEFAULT_CATALOGS[catalog]
        if isdir(path)
            return path
        end
    end

    # 2. Environment variable override (e.g. THORLABS_CATALOG_PATH)
    env_var = string(uppercase(string(catalog)), "_CATALOG_PATH")
    if haskey(ENV, env_var) && isdir(ENV[env_var])
        return normpath(ENV[env_var])
    end

    # 3. Standard fallback relative paths
    fallbacks = [
        normpath(joinpath(@__DIR__, "..", "test_data", string(catalog), "zmx")),
        normpath(joinpath(@__DIR__, "..", "test_data", string(catalog))),
        normpath(joinpath(pwd(), "test_data", string(catalog), "zmx")),
        normpath(joinpath(pwd(), "test_data", string(catalog))),
    ]
    for fb in fallbacks
        if isdir(fb)
            return fb
        end
    end

    registered = get(DEFAULT_CATALOGS, catalog, "<none>")
    throw(ArgumentError("Catalog :$catalog directory not found (registered at '$registered'). Use `register_catalog!(:$catalog, path)` to set the directory."))
end

function get_catalog_path(catalog::AbstractString)::String
    sym = Symbol(lowercase(catalog))
    if haskey(DEFAULT_CATALOGS, sym) || sym == :thorlabs
        try
            return get_catalog_path(sym)
        catch
        end
    end
    if isdir(catalog)
        return normpath(abspath(catalog))
    end
    throw(ArgumentError("Catalog directory '$catalog' does not exist."))
end

"""
    list_catalog_lenses(catalog::Union{Symbol, AbstractString})::Vector{String}

Returns a sorted list of all available lens item names (without `.zmx` extension) in `catalog`.
"""
function list_catalog_lenses(catalog::Union{Symbol, AbstractString})::Vector{String}
    cat_dir = get_catalog_path(catalog)
    entries = String[]
    for (root, _, files) in walkdir(cat_dir)
        for f in files
            if endswith(lowercase(f), ".zmx")
                push!(entries, f[1:end-4])
            end
        end
    end
    return sort!(unique(entries))
end

"""
    find_catalog_lenses(catalog::Union{Symbol, AbstractString}, pattern::Union{Regex, AbstractString})::Vector{String}

Finds all lens item names in `catalog` matching `pattern` (regular expression or case-insensitive substring).
"""
function find_catalog_lenses(catalog::Union{Symbol, AbstractString}, pattern::Regex)::Vector{String}
    all_lenses = list_catalog_lenses(catalog)
    return filter(name -> occursin(pattern, name), all_lenses)
end

function find_catalog_lenses(catalog::Union{Symbol, AbstractString}, pattern::AbstractString)::Vector{String}
    all_lenses = list_catalog_lenses(catalog)
    p_lower = lowercase(pattern)
    return filter(name -> occursin(p_lower, lowercase(name)), all_lenses)
end

"""
    find_catalog_lens_file(catalog::Union{Symbol, AbstractString}, lens_name::AbstractString)::String

Resolves the full filepath of `lens_name` in `catalog`. Handles exact matches, case-insensitive
matching, optional `.zmx` extension, and fuzzy normalization.
"""
function find_catalog_lens_file(catalog::Union{Symbol, AbstractString}, lens_name::AbstractString)::String
    cat_dir = get_catalog_path(catalog)
    clean_name = endswith(lowercase(lens_name), ".zmx") ? lens_name[1:end-4] : lens_name

    # 1. Exact match with .zmx
    direct = joinpath(cat_dir, clean_name * ".zmx")
    if isfile(direct)
        return direct
    end

    # 2. Case-insensitive exact match
    target_lower = lowercase(clean_name * ".zmx")
    for (root, _, files) in walkdir(cat_dir)
        for f in files
            if lowercase(f) == target_lower
                return joinpath(root, f)
            end
        end
    end

    # 3. Normalized alphanumeric match (ignores hyphens, underscores, spaces)
    target_norm = replace(lowercase(clean_name), r"[^a-z0-9]" => "")
    for (root, _, files) in walkdir(cat_dir)
        for f in files
            if endswith(lowercase(f), ".zmx")
                fn_norm = replace(lowercase(f[1:end-4]), r"[^a-z0-9]" => "")
                if fn_norm == target_norm
                    return joinpath(root, f)
                end
            end
        end
    end

    # If not found, provide helpful close suggestions if any
    matches = find_catalog_lenses(catalog, clean_name)
    hint = isempty(matches) ? "" : " Did you mean one of: $(join(first(matches, 5), ", "))?"
    throw(ArgumentError("Lens '$(lens_name)' not found in catalog '$catalog'.$hint"))
end

"""
    load_lens_from_zmx_cat(
        catalog::Union{Symbol, AbstractString},
        lens_name::AbstractString;
        position::Union{Nothing, Real, AbstractVector{<:Real}} = nothing,
        configuration::Int = 1,
        return_result::Bool = false
    )

Loads a specific lens from an optical component catalog (e.g. `:thorlabs`) and constructs
the native BeamletOptics.jl lens object (`Lens`, `DoubletLens`, `TripletLens`).

# Arguments
- `catalog`: Catalog identifier (`:thorlabs` or a direct directory path).
- `lens_name`: Item name or part number (e.g. `"AC254-050-A"`, `"LA1951"`, case-insensitive).
- `position`: Optional position offset. If a `Real` is provided, sets the front vertex axial
  position along the optical axis (+Y). If a 3-element vector `[x, y, z]` is provided, translates
  the front vertex to that coordinate.
- `configuration`: Multi-configuration index (default: 1).
- `return_result`: If `true`, returns a tuple `(lens, res::BMOImportResult)`. If `false` (default),
  returns just the BeamletOptics lens object.

# Returns
A concrete `BeamletOptics.AbstractObject` representing the imported lens.
"""
function load_lens_from_zmx_cat(
    catalog::Union{Symbol, AbstractString},
    lens_name::AbstractString;
    position::Union{Nothing, Real, AbstractVector{<:Real}} = nothing,
    configuration::Int = 1,
    return_result::Bool = false
)
    filepath = find_catalog_lens_file(catalog, lens_name)
    res = import_zmx(filepath; configuration = configuration)

    # Extract optical lenses (exclude detectors and stops)
    optical_lenses = [obj for obj in res.system.objects if !(obj isa BeamletOptics.AbstractDetector)]

    if isempty(optical_lenses)
        throw(ErrorException("No optical refractive element found in $(basename(filepath))."))
    end

    lens_obj = length(optical_lenses) == 1 ? optical_lenses[1] : optical_lenses

    # Optional translation to requested position
    if position !== nothing
        # Determine the initial front vertex axial position (+Y)
        init_y = 0.0
        for el in res.elements
            if !(el isa ZmxStop) && !(el isa ZmxDetector)
                init_y = el.axial_position
                break
            end
        end

        delta = if position isa Real
            [0.0, Float64(position) - init_y, 0.0]
        elseif position isa AbstractVector{<:Real}
            if length(position) != 3
                throw(ArgumentError("Position vector must have length 3 [x, y, z]. Got length $(length(position))."))
            end
            [Float64(position[1]), Float64(position[2]) - init_y, Float64(position[3])]
        else
            throw(ArgumentError("Position must be a Real (Y position) or 3-element vector [x, y, z]."))
        end

        if !iszero(norm(delta))
            if lens_obj isa AbstractVector
                for lo in lens_obj
                    translate3d!(lo, delta)
                end
            else
                translate3d!(lens_obj, delta)
            end
        end
    end

    if return_result
        return (lens_obj, res)
    else
        return lens_obj
    end
end
