#!/usr/bin/env julia

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using BeamletOptics
using BeamletOpticsZMX

include("common.jl")

const zmx_file = joinpath(@__DIR__, "..", "test_data", "dan_reiley", "photographic-lenses-zoom", "US05867325-1.zmx")

println("Lade Foto-Zoomobjektiv (US Patent 5,867,325)...")
res = import_zmx(zmx_file)

display_system_summary(res, "Foto-Zoomobjektive (Photographic Zooms)")

λ = res.zmx_system.wavelengths[res.zmx_system.primary_wavelength_idx]
beam_diam = 0.015 # 15 mm
source_y = -0.01
source = CollimatedSource([0.0, source_y, 0.0], [0.0, 1.0, 0.0], beam_diam, λ, num_rays=200, num_rings=6)

println("Starte Strahlverfolgung (200 Strahlen)...")
solve_system!(res.system, source)

if res.detector !== nothing
    hits = spot_diagram(res.detector)
    print_ascii_spot_diagram(hits)
end
