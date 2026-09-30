#!/usr/bin/env julia

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using BeamletOptics
using BeamletOpticsZMX

include("common.jl")

const zmx_file = joinpath(@__DIR__, "..", "test_data", "dan_reiley", "eyepieces", "UK565851-1.zmx")

println("Lade Okular-Design (Eyepiece)...")
res = import_zmx(zmx_file)

display_system_summary(res, "Okulare (Eyepieces)")

# Okular: In Zemax rückwärts gerechnet (vom Auge zur Zwischenbildebene)
# ENPD = 8.0 mm, Augenabstand (Eye Relief) = 5.08 mm
λ = res.zmx_system.wavelengths[res.zmx_system.primary_wavelength_idx]
beam_diam = 0.008 # 8.0 mm (entspricht ENPD)
source_y = -0.005 # 5 mm vor Augenpupille
source = CollimatedSource([0.0, source_y, 0.0], [0.0, 1.0, 0.0], beam_diam, λ, num_rays=150, num_rings=5)

println("Starte Strahlverfolgung (150 Strahlen)...")
solve_system!(res.system, source)

if res.detector !== nothing
    plot_bmo_spot_diagram(res, source, joinpath(@__DIR__, "output", "02_eyepiece_spot.png");
                          title="UK565851-1 Okular", category="Okulare")
end
