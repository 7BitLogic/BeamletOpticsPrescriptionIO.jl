#!/usr/bin/env julia

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using BeamletOptics
using BeamletOpticsZMX

include("common.jl")

const zmx_file = joinpath(@__DIR__, "..", "test_data", "dan_reiley", "photographic-lenses-prime", "US00583336-2-scaled.zmx")

println("Lade Foto-Festbrennweite (Rudolph Double Gauss)...")
res = import_zmx(zmx_file)

display_system_summary(res, "Foto-Festbrennweiten (Photographic Primes)")

# Foto-Festbrennweite (Double Gauss):
# Aperturblende liegt im Zentrum zwischen den beiden Dubletts (Diam = 4.2 mm)
λ = res.zmx_system.wavelengths[res.zmx_system.primary_wavelength_idx]
beam_diam = 0.005 # 5 mm
source_y = -0.005 # 5 mm vor Frontlinse
source = CollimatedSource([0.0, source_y, 0.0], [0.0, 1.0, 0.0], beam_diam, λ, num_rays=200, num_rings=6)

println("Starte Strahlverfolgung (200 Strahlen)...")
solve_system!(res.system, source)

if res.detector !== nothing
    plot_bmo_spot_diagram(res, source, joinpath(@__DIR__, "output", "04_photographic_prime_spot.png");
                          title="US00583336-2 Rudolph Double Gauss", category="Foto-Festbrennweiten")
end
