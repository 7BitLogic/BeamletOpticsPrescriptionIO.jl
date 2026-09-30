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

# Foto-Zoom: ZMX-Spezifikation definiert ENPD = 7.0 mm (Eintrittspupille)
λ = res.zmx_system.wavelengths[res.zmx_system.primary_wavelength_idx]
beam_diam = 0.007 # 7.0 mm (entspricht exakt ENPD)
source_y = -0.005 # 5 mm vor Frontlinse
source = CollimatedSource([0.0, source_y, 0.0], [0.0, 1.0, 0.0], beam_diam, λ, num_rays=200, num_rings=6)

println("Starte Strahlverfolgung (200 Strahlen)...")
solve_system!(res.system, source)

if res.detector !== nothing
    plot_bmo_spot_diagram(res, source, joinpath(@__DIR__, "output", "05_photographic_zoom_spot.png");
                          title="US05867325-1 Foto-Zoom", category="Foto-Zoomobjektive")
end
