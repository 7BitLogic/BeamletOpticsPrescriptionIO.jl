#!/usr/bin/env julia

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using BeamletOptics
using BeamletOpticsZMX

include("common.jl")

const zmx_file = joinpath(@__DIR__, "..", "test_data", "dan_reiley", "projectors", "US03126786-1.zmx")

println("Lade Projektor-Objektiv (US Patent 3,126,786)...")
res = import_zmx(zmx_file)

display_system_summary(res, "Projektoren (Projectors)")

# Projektor: Große Aperturblende (Diam = 39.9 mm), Frontlinsendurchmesser 55 mm
λ = res.zmx_system.wavelengths[res.zmx_system.primary_wavelength_idx]
beam_diam = 0.025 # 25 mm
source_y = -0.015 # 15 mm vor Frontlinse
source = CollimatedSource([0.0, source_y, 0.0], [0.0, 1.0, 0.0], beam_diam, λ, num_rays=200, num_rings=6)

println("Starte Strahlverfolgung (200 Strahlen)...")
solve_system!(res.system, source)

if res.detector !== nothing
    plot_bmo_spot_diagram(res, source, joinpath(@__DIR__, "output", "06_projector_spot.png");
                          title="US03126786-1 Projektor", category="Projektoren")
end
