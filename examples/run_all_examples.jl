#!/usr/bin/env julia

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Printf

const example_scripts = [
    ("Endoscopes", "01_endoscope.jl"),
    ("Eyepieces", "02_eyepiece.jl"),
    ("Microscope Objectives", "03_microscope_objective.jl"),
    ("Photographic Primes", "04_photographic_prime.jl"),
    ("Photographic Zooms", "05_photographic_zoom.jl"),
    ("Projectors", "06_projector.jl"),
    ("Scan Lenses", "07_scan_lens.jl"),
    ("Spectrometers", "08_spectrometer.jl"),
    ("Telescopes", "09_telescope.jl")
]

println("="^80)
println("Running all 9 category examples from Dan Reiley database...")
println("="^80)

results = []

for (cat, script) in example_scripts
    path = joinpath(@__DIR__, script)
    t0 = time()
    try
        # Run script
        include(path)
        dt = time() - t0
        push!(results, (cat, script, "SUCCESS", dt))
    catch e
        dt = time() - t0
        push!(results, (cat, script, "ERROR: $e", dt))
    end
end

println("\n" * "="^80)
println("SUMMARY OF ALL 9 CATEGORIES")
println("="^80)
@printf("%-24s %-28s %-15s %s\n", "Category", "Script", "Status", "Elapsed")
println("-"^80)
for (cat, script, status, dt) in results
    @printf("%-24s %-28s %-15s %.2f s\n", cat, script, status, dt)
end
println("="^80)
