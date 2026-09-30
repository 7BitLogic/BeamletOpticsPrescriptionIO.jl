#!/usr/bin/env julia

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using BeamletOptics
using BeamletOpticsZMX

include("common.jl")

const zmx_file = joinpath(@__DIR__, "..", "test_data", "dan_reiley", "microscope-objectives", "JP1985-063512-1.zmx")

println("Lade Mikroskop-Objektiv...")
res = import_zmx(zmx_file)

display_system_summary(res, "Mikroskop-Objektive (Microscope Objectives)")

# Mikroskop-Objektiv: Surface 1 ist Aperturblende (y=0) mit Versatz DISZ = -15.12 mm.
# Erste Linse beginnt bei y = -15.12 mm. Quelle muss davor liegen (y = -20 mm).
λ = res.zmx_system.wavelengths[res.zmx_system.primary_wavelength_idx]
beam_diam = 0.0022 # 2.2 mm (entspricht Apertur/Objektfeld)
source_y = -0.020  # 20 mm vor Blende, 5 mm vor Frontlinse
source = CollimatedSource([0.0, source_y, 0.0], [0.0, 1.0, 0.0], beam_diam, λ, num_rays=150, num_rings=5)

println("Starte Strahlverfolgung (150 Strahlen)...")
solve_system!(res.system, source)

if res.detector !== nothing
    plot_bmo_spot_diagram(res, source, joinpath(@__DIR__, "output", "03_microscope_objective_spot.png");
                          title="JP1985-063512-1 Mikroskop-Objektiv", category="Mikroskop-Objektive")
end
