#!/usr/bin/env julia

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using BeamletOptics
using BeamletOpticsPrescriptionIO

include("common.jl")

const zmx_file = joinpath(@__DIR__, "..", "test_data", "dan_reiley", "eyepieces", "UK565851-1.zmx")

println("Loading Eyepiece design...")
res = import_zmx(zmx_file)

display_system_summary(res, "Eyepieces")

# Eyepiece: Computed in reverse in Zemax (from eye pupil to intermediate image plane)
# ENPD = 8.0 mm, Eye Relief = 5.08 mm
λ = res.zmx_system.wavelengths[res.zmx_system.primary_wavelength_idx]
beam_diam = 0.008 # 8.0 mm (matches ENPD)
source_y = -0.005 # 5 mm before eye pupil
source = CollimatedSource([0.0, source_y, 0.0], [0.0, 1.0, 0.0], beam_diam, λ, num_rays=150, num_rings=5)

println("Starting ray trace (150 rays)...")
solve_system!(res.system, source)

if res.detector !== nothing
    plot_bmo_spot_diagram(res, source, joinpath(@__DIR__, "output", "02_eyepiece_spot.png");
                          title="UK565851-1 Eyepiece", category="Eyepieces")
end
