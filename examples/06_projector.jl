#!/usr/bin/env julia

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using BeamletOptics
using BeamletOpticsPrescriptionIO

include("common.jl")

const zmx_file = joinpath(@__DIR__, "..", "test_data", "dan_reiley", "projectors", "US03126786-1.zmx")

println("Loading Projector Lens design (US Patent 3,126,786)...")
res = import_zmx(zmx_file)

display_system_summary(res, "Projectors")

# Projector: Large aperture stop (Diam = 39.9 mm), front lens diameter 55 mm
λ = res.zmx_system.wavelengths[res.zmx_system.primary_wavelength_idx]
beam_diam = 0.025 # 25 mm
source_y = -0.015 # 15 mm before front lens
source = CollimatedSource([0.0, source_y, 0.0], [0.0, 1.0, 0.0], beam_diam, λ, num_rays=200, num_rings=6)

println("Starting ray trace (200 rays)...")
solve_system!(res.system, source)

if res.detector !== nothing
    plot_bmo_spot_diagram(res, source, joinpath(@__DIR__, "output", "06_projector_spot.png");
                          title="US03126786-1 Projector", category="Projectors")
end
