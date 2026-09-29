#!/usr/bin/env julia

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using BeamletOptics
using BeamletOpticsZMX

include("common.jl")

const zmx_file = joinpath(@__DIR__, "..", "test_data", "dan_reiley", "telescopes", "Figure1.zmx")

println("Lade Spiegelteleskop (Schmidt-Cassegrain / Mangin System)...")
res = import_zmx(zmx_file)

display_system_summary(res, "Teleskope (Telescopes)")

λ = res.zmx_system.wavelengths[res.zmx_system.primary_wavelength_idx]
beam_diam = 0.050 # 50 mm
source_y = -0.02
source = CollimatedSource([0.0, source_y, 0.0], [0.0, 1.0, 0.0], beam_diam, λ, num_rays=200, num_rings=6)

println("Starte Strahlverfolgung (200 Strahlen)...")
solve_system!(res.system, source)

if res.detector !== nothing
    plot_bmo_spot_diagram(res, source, joinpath(@__DIR__, "output", "09_telescope_spot.png");
                          title="Figure1 Spiegelteleskop", category="Teleskope")
end
