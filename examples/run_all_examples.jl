#!/usr/bin/env julia

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Printf

const example_scripts = [
    ("Endoskope", "01_endoscope.jl"),
    ("Okulare", "02_eyepiece.jl"),
    ("Mikroskop-Objektive", "03_microscope_objective.jl"),
    ("Foto-Festbrennweiten", "04_photographic_prime.jl"),
    ("Foto-Zoomobjektive", "05_photographic_zoom.jl"),
    ("Projektoren", "06_projector.jl"),
    ("Scan-Objektive", "07_scan_lens.jl"),
    ("Spektrometer", "08_spectrometer.jl"),
    ("Teleskope", "09_telescope.jl")
]

println("="^80)
println("Führe alle 9 Kategorie-Beispiele von Dan Reiley aus...")
println("="^80)

results = []

for (cat, script) in example_scripts
    path = joinpath(@__DIR__, script)
    t0 = time()
    try
        # Run script
        include(path)
        dt = time() - t0
        push!(results, (cat, script, "ERFOLGREICH", dt))
    catch e
        dt = time() - t0
        push!(results, (cat, script, "FEHLER: $e", dt))
    end
end

println("\n" * "="^80)
println("ZUSAMMENFASSUNG ALLER 9 KATEGORIEN")
println("="^80)
@printf("%-24s %-28s %-15s %s\n", "Kategorie", "Skript", "Status", "Laufzeit")
println("-"^80)
for (cat, script, status, dt) in results
    @printf("%-24s %-28s %-15s %.2f s\n", cat, script, status, dt)
end
println("="^80)
