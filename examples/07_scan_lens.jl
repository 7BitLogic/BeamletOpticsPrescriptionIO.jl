#!/usr/bin/env julia

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using BeamletOptics
using BeamletOpticsPrescriptionIO

include("common.jl")

const zmx_file = joinpath(@__DIR__, "..", "test_data", "dan_reiley", "scan-lenses", "US05025268-1.ZMX")

println("Loading Scan Lens design (f-theta Scan Lens)...")
res = import_zmx(zmx_file)

display_system_summary(res, "Scan Lenses")

# Scan Lens: ENPD = 1.0 mm (laser beam diameter on scanning mirror)
# Surface 1 is aperture stop at y = 0 with Diam = 0.5 mm
λ = res.zmx_system.wavelengths[res.zmx_system.primary_wavelength_idx]
beam_diam = 0.001 # 1.0 mm (matches ENPD)
source_y = -0.005 # 5 mm before scan stop
source = CollimatedSource([0.0, source_y, 0.0], [0.0, 1.0, 0.0], beam_diam, λ, num_rays=150, num_rings=5)

println("Starting ray trace (150 rays)...")
solve_system!(res.system, source)

if res.detector !== nothing
    plot_bmo_spot_diagram(res, source, joinpath(@__DIR__, "output", "07_scan_lens_spot.png");
                          title="US05025268-1 Scan Lens", category="Scan Lenses")
end
