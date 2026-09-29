#!/usr/bin/env julia

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using BeamletOptics
using BeamletOpticsZMX

include("common.jl")

const zmx_file = joinpath(@__DIR__, "..", "test_data", "dan_reiley", "scan-lenses", "US05025268-1.ZMX")

println("Lade Scan-Objektiv (f-theta Scan Lens)...")
res = import_zmx(zmx_file)

display_system_summary(res, "Scan-Objektive (Scan Lenses)")

λ = res.zmx_system.wavelengths[res.zmx_system.primary_wavelength_idx]
beam_diam = 0.006 # 6 mm
source_y = -0.01
source = CollimatedSource([0.0, source_y, 0.0], [0.0, 1.0, 0.0], beam_diam, λ, num_rays=150, num_rings=5)

println("Starte Strahlverfolgung (150 Strahlen)...")
solve_system!(res.system, source)

if res.detector !== nothing
    plot_bmo_spot_diagram(res, source, joinpath(@__DIR__, "output", "07_scan_lens_spot.png");
                          title="US05025268-1 Scan-Objektiv", category="Scan-Objektive")
end
