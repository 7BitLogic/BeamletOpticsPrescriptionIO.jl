"""
Types and data structures for Zemax .zmx optical system representation.
"""

"""
    ZmxSurface

Represents a single optical interface parsed from a Zemax .zmx file.
All physical length quantities (radius, thickness, semi-diameter) are stored in SI units (meters).
"""
Base.@kwdef struct ZmxSurface
    index::Int
    surface_type::Symbol = :STANDARD # :STANDARD, :EVENASPH, :COORDBRK, :TOROIDAL, etc.
    comment::String = ""
    curvature::Float64 = 0.0          # Curvature c in 1/m (c = 1/R)
    radius::Float64 = Inf             # Radius of curvature R = 1/c in m (Inf if flat)
    thickness::Float64 = 0.0          # Axial thickness to next surface in m
    glass_name::String = ""           # Glass identifier (e.g. "N-BK7", "AIR", "___BLANK", "")
    nd::Float64 = 1.0                 # Refractive index at d-line (587.56 nm)
    vd::Float64 = 0.0                 # Abbe number (0.0 if not specified)
    semi_diameter::Float64 = 0.0      # Clear semi-diameter in m
    conic::Float64 = 0.0              # Conic constant (kappa / k)
    parms::Vector{Float64} = Float64[]# Parameter array (e.g. aspheric coeffs or decenter/tilt)
    is_stop::Bool = false             # Aperture stop flag
    is_mirror::Bool = false           # Reflective surface flag
    decenter_x::Float64 = 0.0         # Coordinate break X decenter [m]
    decenter_y::Float64 = 0.0         # Coordinate break Y decenter [m]
    tilt_x::Float64 = 0.0             # Coordinate break X tilt [rad]
    tilt_y::Float64 = 0.0             # Coordinate break Y tilt [rad]
    tilt_z::Float64 = 0.0             # Coordinate break Z tilt [rad]
end

"""
    ZmxSystem

Contains the parsed metadata and sequential surface list of a Zemax file.
"""
Base.@kwdef struct ZmxSystem
    title::String = ""
    unit::Symbol = :mm                # :mm, :inch, :cm, :meter
    scale_to_m::Float64 = 1e-3        # Multiplier to convert file units to meters
    version::String = ""
    wavelengths::Vector{Float64} = [587.56e-9] # Wavelengths in meters
    spectral_weights::Vector{Float64} = [1.0]
    primary_wavelength_idx::Int = 1
    enpd::Float64 = 0.0               # Entrance Pupil Diameter [m]
    obna::Float64 = 0.0               # Object Numerical Aperture
    object_distance::Float64 = Inf    # Distance from Object (Surf 0) to first optical surface [m]
    fields::Vector{Tuple{Float64, Float64, Float64}} = Tuple{Float64, Float64, Float64}[] # (x, y, weight)
    field_type::Symbol = :angle       # :angle (degrees), :object_height, :paraxial_image_height, :real_image_height
    surfaces::Vector{ZmxSurface} = ZmxSurface[]
end

# High-level BMO Element abstractions
abstract type ZmxElement end

"""
    ZmxSinglet <: ZmxElement

Represents a single optical lens formed by two surfaces (front and back) enclosing a glass medium.
"""
Base.@kwdef struct ZmxSinglet <: ZmxElement
    name::String = "Singlet"
    front_surface::ZmxSurface
    back_surface::ZmxSurface
    center_thickness::Float64   # [m]
    diameter::Float64           # [m]
    glass_name::String = ""
    nd::Float64 = 1.0
    vd::Float64 = 0.0
    axial_position::Float64 = 0.0 # Front vertex position along BMO optical axis (+Y) [m]
end

"""
    ZmxDoublet <: ZmxElement

Represents a cemented doublet lens (3 surfaces, 2 glass media).
"""
Base.@kwdef struct ZmxDoublet <: ZmxElement
    name::String = "Doublet"
    surface1::ZmxSurface
    surface2::ZmxSurface
    surface3::ZmxSurface
    thickness1::Float64
    thickness2::Float64
    diameter::Float64
    glass1::String = ""
    glass2::String = ""
    nd1::Float64 = 1.0
    vd1::Float64 = 0.0
    nd2::Float64 = 1.0
    vd2::Float64 = 0.0
    axial_position::Float64 = 0.0 # Front vertex position along +Y [m]
end

"""
    ZmxTriplet <: ZmxElement

Represents a cemented triplet lens (4 surfaces, 3 glass media).
"""
Base.@kwdef struct ZmxTriplet <: ZmxElement
    name::String = "Triplet"
    surface1::ZmxSurface
    surface2::ZmxSurface
    surface3::ZmxSurface
    surface4::ZmxSurface
    thickness1::Float64
    thickness2::Float64
    thickness3::Float64
    diameter::Float64
    glass1::String = ""
    glass2::String = ""
    glass3::String = ""
    nd1::Float64 = 1.0
    vd1::Float64 = 0.0
    nd2::Float64 = 1.0
    vd2::Float64 = 0.0
    nd3::Float64 = 1.0
    vd3::Float64 = 0.0
    axial_position::Float64 = 0.0 # Front vertex position along +Y [m]
end

"""
    ZmxMirror <: ZmxElement

Represents a reflective optical mirror.
"""
Base.@kwdef struct ZmxMirror <: ZmxElement
    name::String = "Mirror"
    surface::ZmxSurface
    diameter::Float64
    axial_position::Float64 = 0.0 # Vertex position along +Y [m]
end

"""
    ZmxStop <: ZmxElement

Represents an aperture stop.
"""
Base.@kwdef struct ZmxStop <: ZmxElement
    name::String = "Stop"
    surface::ZmxSurface
    diameter::Float64
    axial_position::Float64 = 0.0 # Position along +Y [m]
end

"""
    ZmxDetector <: ZmxElement

Represents the image plane / detector.
"""
Base.@kwdef struct ZmxDetector <: ZmxElement
    name::String = "Detector"
    surface::ZmxSurface
    diameter::Float64
    axial_position::Float64 = 0.0 # Position along +Y [m]
end

"""
    BMOImportResult

Result of importing a Zemax file into BeamletOptics.jl.
"""
struct BMOImportResult
    system::Any                   # BeamletOptics.System
    detector::Union{Nothing, Any} # BeamletOptics.Detector (if created)
    stop::Union{Nothing, ZmxStop} # Aperture Stop element (if present)
    elements::Vector{ZmxElement}
    zmx_system::ZmxSystem
end

BMOImportResult(system, detector, elements, zmx_system) = BMOImportResult(system, detector, nothing, elements, zmx_system)

