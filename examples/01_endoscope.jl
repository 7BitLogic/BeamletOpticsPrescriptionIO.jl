#!/usr/bin/env julia

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using BeamletOptics
using BeamletOpticsZMX

include("common.jl")

const zmx_file = joinpath(@__DIR__, "..", "test_data", "dan_reiley", "endoscopes", "US05724190-1.ZMX")

println("Lade Endoskop-Design...")
res = import_zmx(zmx_file)

display_system_summary(res, "Endoskope (Endoscopes)")

# Endoskop: kleiner Eingangsstrahldurchmesser (z.B. 0.8 mm)
λ = res.zmx_system.wavelengths[res.zmx_system.primary_wavelength_idx]
beam_diam = 0.0008 # 0.8 mm
source_y = -0.002 # 2 mm vor erster Fläche
source = CollimatedSource([0.0, source_y, 0.0], [0.0, 1.0, 0.0], beam_diam, λ, num_rays=150, num_rings=5)

println("Starte Strahlverfolgung (150 Strahlen)...")
solve_system!(res.system, source)

if res.detector !== nothing
    plot_bmo_spot_diagram(res, source, joinpath(@__DIR__, "output", "01_endoscope_spot.png");
                          title="US05724190-1 Endoskop", category="Endoskope")
end
