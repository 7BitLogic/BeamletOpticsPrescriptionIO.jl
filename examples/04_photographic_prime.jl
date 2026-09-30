#!/usr/bin/env julia

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using BeamletOptics
using BeamletOpticsPrescriptionIO

include("common.jl")

const zmx_file = joinpath(@__DIR__, "..", "test_data", "dan_reiley", "photographic-lenses-prime", "US00583336-2-scaled.zmx")

println("Loading Photographic Prime design (Rudolph Double Gauss)...")
res = import_zmx(zmx_file)

display_system_summary(res, "Photographic Primes")

# Photographic Prime (Double Gauss):
# Aperture stop is located centrally between the two doublets (Diam = 4.2 mm)
λ = res.zmx_system.wavelengths[res.zmx_system.primary_wavelength_idx]
beam_diam = 0.005 # 5 mm
source_y = -0.005 # 5 mm before front lens
source = CollimatedSource([0.0, source_y, 0.0], [0.0, 1.0, 0.0], beam_diam, λ, num_rays=200, num_rings=6)

println("Starting ray trace (200 rays)...")
solve_system!(res.system, source)

if res.detector !== nothing
    plot_bmo_spot_diagram(res, source, joinpath(@__DIR__, "output", "04_photographic_prime_spot.png");
                          title="US00583336-2 Rudolph Double Gauss", category="Photographic Primes")
end
