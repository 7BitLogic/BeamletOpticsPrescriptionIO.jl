#!/usr/bin/env julia

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using BeamletOptics
using BeamletOpticsPrescriptionIO

include("common.jl")

const zmx_file = joinpath(@__DIR__, "..", "test_data", "dan_reiley", "telescopes", "Figure1.zmx")

println("Loading Reflecting Telescope design (Schmidt-Cassegrain / Mangin System)...")
res = import_zmx(zmx_file)

display_system_summary(res, "Telescopes")

# Reflecting Telescope (Dennis Gabor catadioptric system with meniscus corrector and spherical mirror):
# ZMX header: ENPD = 78.4 mm. Entrance aperture / corrector lens Diam ~50 mm, primary mirror Diam 64 mm.
λ = res.zmx_system.wavelengths[res.zmx_system.primary_wavelength_idx]
beam_diam = 0.075 # 75 mm
source_y = -0.020 # 20 mm before entrance aperture
source = CollimatedSource([0.0, source_y, 0.0], [0.0, 1.0, 0.0], beam_diam, λ, num_rays=200, num_rings=6)

println("Starting ray trace (200 rays)...")
solve_system!(res.system, source)

if res.detector !== nothing
    plot_bmo_spot_diagram(res, source, joinpath(@__DIR__, "output", "09_telescope_spot.png");
                          title="Figure1 Reflecting Telescope", category="Telescopes")
end
