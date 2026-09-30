#!/usr/bin/env julia

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using BeamletOptics
using BeamletOpticsPrescriptionIO

include("common.jl")

const zmx_file = joinpath(@__DIR__, "..", "test_data", "dan_reiley", "spectro", "US05644396-1.ZMX")

println("Loading Spectrometer design (US Patent 5,644,396)...")
res = import_zmx(zmx_file)

display_system_summary(res, "Spectrometers")

# Spectrometer: Folded catadioptric Littrow system (double-pass)
# OBNA = 0.41, finite conjugate distance Surf 0 DISZ = 24.04 mm
λ = res.zmx_system.wavelengths[res.zmx_system.primary_wavelength_idx]
beam_diam = 0.004 # 4 mm beam bundle
source_y = -0.010 # 10 mm before front lens
source = CollimatedSource([0.0, source_y, 0.0], [0.0, 1.0, 0.0], beam_diam, λ, num_rays=150, num_rings=5)

println("Starting ray trace (150 rays)...")
solve_system!(res.system, source)

if res.detector !== nothing
    plot_bmo_spot_diagram(res, source, joinpath(@__DIR__, "output", "08_spectrometer_spot.png");
                          title="US05644396-1 Spectrometer", category="Spectrometers")
end
