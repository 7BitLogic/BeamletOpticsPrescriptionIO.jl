module BeamletOpticsZMX

using BeamletOptics
using LinearAlgebra
using Printf
using StringEncodings
using RefractiveIndex

include("Types.jl")
include("Parser.jl")
include("GlassCatalog.jl")
include("ElementGrouper.jl")
include("Builder.jl")
include("CodeGen.jl")

# Exports
export import_zmx, parse_zmx, parse_zmx_content, group_elements, build_bmo_system, build_surface, build_bmo_object
export generate_bmo_code, generate_bmo_script
export resolve_refractive_index, cauchy_dispersion
export ZmxSurface, ZmxSystem, ZmxElement, ZmxSinglet, ZmxDoublet, ZmxTriplet, ZmxMirror, ZmxStop, ZmxDetector, BMOImportResult

end # module BeamletOpticsZMX
