#!/usr/bin/env julia

"""
Example: Loading & Analyzing Lenses from a Local Catalog (e.g. Thorlabs)
========================================================================

Demonstrates:
  1. Specifying / registering the local catalog directory on disk
  2. Discovering and searching available lenses in the catalog
  3. Loading an achromatic doublet and automatically placing the detector
  4. Tracing a collimated beam and calculating the RMS spot size
  5. Optimizing focus with `find_best_focus`
"""

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using BeamletOptics
using BeamletOpticsPrescriptionIO
using Printf

# ==============================================================================
# Step 1: Tell BeamletOptics where the catalog is stored on your disk
# ==============================================================================
# You have three flexible options:
#
# Option A: Register your custom folder once for the session:
#   register_catalog!(:thorlabs, "/path/to/your/thorlabs/zmx")
#
# Option B: Pass your folder path directly as a string to functions:
#   load_lens_from_zmx_cat("/path/to/your/thorlabs/zmx", "AC254-050-A")
#
# Option C: Use the default repository path (populated by scripts/download_test_data.jl):
#   test_data/thorlabs/zmx

catalog_path = joinpath(@__DIR__, "..", "test_data", "thorlabs", "zmx")

if !isdir(catalog_path)
    println("""
Catalog directory not found at: $catalog_path
To download the Thorlabs catalog files, run:
    julia --project=. scripts/download_test_data.jl
Or register an existing folder on your disk:
    register_catalog!(:thorlabs, "/path/to/your/thorlabs/zmx")
""")
    exit(1)
end

# Register the path
register_catalog!(:thorlabs, catalog_path)
println("Using Thorlabs catalog at: ", catalog_path)

# ==============================================================================
# Step 2: Search for available lenses in the catalog
# ==============================================================================
println("\nSearching for 50 mm achromatic doublets (AC254-050)...")
matches = find_catalog_lenses(:thorlabs, "AC254-050")
for m in matches
    println("  - $m")
end

# ==============================================================================
# Step 3: Load the lens and its nominal detector from the prescription
# ==============================================================================
# By passing `return_result=true`, we get:
#   - `lens`: The 3D physical BeamletOptics DoubletLens object
#   - `res`:  The BMOImportResult containing the nominal detector at the back focal length
part_number = "AC254-050-A"
println("\nLoading lens '$part_number' at position y = 100 mm...")
lens, res = load_lens_from_zmx_cat(:thorlabs, part_number, position=0.10, return_result=true)

println("Lens type:        ", typeof(lens))
println("Nominal detector: y = ", @sprintf("%.2f mm", position(res.detector)[2] * 1e3))

# ==============================================================================
# Step 4: Trace a collimated ray bundle
# ==============================================================================
# Recommend a ray source matching the lens diameter (25.4 mm clear aperture, 587.6 nm d-line)
source = suggest_source(res; num_rays=120, num_rings=5)
solve_system!(res.system, source)

hits = spot_diagram(res.detector)
println("Rays intercepted: ", length(hits))

# ==============================================================================
# Step 5: Focus analysis & RMS spot radius
# ==============================================================================
foc = find_best_focus(res.detector)
println("\nFocus Analysis:")
println("  Nominal detector y: ", @sprintf("%.3f mm (RMS = %.2f µm)", foc.y_nom * 1e3, foc.rms_nom * 1e6))
println("  Optimal focus y:    ", @sprintf("%.3f mm (RMS = %.2f µm)", foc.y_opt * 1e3, foc.rms_opt * 1e6))
println("  Defocus offset Δy:  ", @sprintf("%+.2f µm", foc.delta_y * 1e6))

# Apply optimal focus
refocus!(res)
println("Refocused detector to optimal focal plane!")
