#!/usr/bin/env julia

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using BeamletOptics
using BeamletOpticsPrescriptionIO

include("common.jl")

const zmx_file = joinpath(@__DIR__, "..", "test_data", "dan_reiley", "microscope-objectives", "JP1985-063512-1.zmx")

println("Loading Microscope Objective design...")
res = import_zmx(zmx_file)

display_system_summary(res, "Microscope Objectives")

# Microscope Objective: Surface 1 is aperture stop (y=0) with offset DISZ = -15.12 mm.
# First lens begins at y = -15.12 mm. Source placed upstream at y = -20 mm.
λ = res.zmx_system.wavelengths[res.zmx_system.primary_wavelength_idx]
beam_diam = 0.0022 # 2.2 mm (matches aperture / object field)
source_y = -0.020  # 20 mm before stop, 5 mm before front lens
source = CollimatedSource([0.0, source_y, 0.0], [0.0, 1.0, 0.0], beam_diam, λ, num_rays=150, num_rings=5)

println("Starting ray trace (150 rays)...")
solve_system!(res.system, source)

if res.detector !== nothing
    plot_bmo_spot_diagram(res, source, joinpath(@__DIR__, "output", "03_microscope_objective_spot.png");
                          title="JP1985-063512-1 Microscope Objective", category="Microscope Objectives")
end
