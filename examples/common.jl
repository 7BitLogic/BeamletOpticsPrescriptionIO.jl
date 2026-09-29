"""
Shared utilities for rendering and analyzing optical systems in the example scripts.
"""

using Printf
using LinearAlgebra

function display_system_summary(res, category::String)
    zsys = res.zmx_system
    println("="^80)
    println("KATEGORIE:    $category")
    println("MODELL:       $(isempty(zsys.title) ? "Zemax Design" : zsys.title)")
    println("="^80)
    println("Dateieinheit: $(zsys.unit) (Skalierung zu m: $(zsys.scale_to_m))")
    println("Flächen:      $(length(zsys.surfaces))")
    println("Baugruppen:   $(length(res.elements))")
    println("Wellenlängen: $([@sprintf("%.1f nm", w*1e9) for w in zsys.wavelengths])")
    println("-"^80)
    println("Rekonstruierte BMO-Komponenten:")
    for (i, el) in enumerate(res.elements)
        if isa(el, ZmxSinglet)
            r1 = isinf(el.front_surface.radius) ? "Plano" : @sprintf("%.2f mm", el.front_surface.radius*1e3)
            r2 = isinf(el.back_surface.radius) ? "Plano" : @sprintf("%.2f mm", el.back_surface.radius*1e3)
            println(@sprintf("  [%d] Einzellinse:        %-12s (Glas: %-8s d=%6.2f mm, R1=%9s, R2=%9s, y=%6.2f mm)",
                i, el.name, el.glass_name, el.center_thickness*1e3, r1, r2, el.axial_position*1e3))
        elseif isa(el, ZmxDoublet)
            println(@sprintf("  [%d] Verkittetes Dublett: %-12s (Gläser: %s + %s, y=%6.2f mm)",
                i, el.name, el.glass1, el.glass2, el.axial_position*1e3))
        elseif isa(el, ZmxTriplet)
            println(@sprintf("  [%d] Verkittetes Triplett:%-12s (Gläser: %s + %s + %s, y=%6.2f mm)",
                i, el.name, el.glass1, el.glass2, el.glass3, el.axial_position*1e3))
        elseif isa(el, ZmxMirror)
            r = isinf(el.surface.radius) ? "Plan" : @sprintf("%.2f mm", el.surface.radius*1e3)
            println(@sprintf("  [%d] Spiegel:             %-12s (Radius: %s, Diam=%6.2f mm, y=%6.2f mm)",
                i, el.name, r, el.diameter*1e3, el.axial_position*1e3))
        elseif isa(el, ZmxStop)
            println(@sprintf("  [%d] Aperturblende:       y=%6.2f mm, Durchmesser=%6.2f mm",
                i, el.axial_position*1e3, el.diameter*1e3))
        elseif isa(el, ZmxDetector)
            println(@sprintf("  [%d] Bildsensor/Detektor: y=%6.2f mm, Durchmesser=%6.2f mm",
                i, el.axial_position*1e3, el.diameter*1e3))
        end
    end
    println("="^80)
end

function print_ascii_spot_diagram(hits; width=60, height=18)
    if isempty(hits)
        println("  (Keine Strahlen auf dem Detektor aufgetroffen)")
        return
    end
    
    xs = [h[1] for h in hits]
    ys = [h[2] for h in hits]
    
    # Statistics
    n = length(hits)
    mean_x = sum(xs) / n
    mean_y = sum(ys) / n
    rms_radius = sqrt(sum((x - mean_x)^2 + (y - mean_y)^2 for (x, y) in zip(xs, ys)) / n)
    max_radius = maximum(norm(h) for h in hits)

    min_x, max_x = minimum(xs), maximum(xs)
    min_y, max_y = minimum(ys), maximum(ys)
    
    span_x = max(max_x - min_x, 1e-7)
    span_y = max(max_y - min_y, 1e-7)
    
    grid = fill(' ', height, width)
    
    for (x, y) in zip(xs, ys)
        col = clamp(round(Int, (x - min_x) / span_x * (width - 1)) + 1, 1, width)
        row = clamp(round(Int, (y - min_y) / span_y * (height - 1)) + 1, 1, height)
        grid[row, col] = '*'
    end
    
    # Mark centroid with '+'
    c_col = clamp(round(Int, (mean_x - min_x) / span_x * (width - 1)) + 1, 1, width)
    c_row = clamp(round(Int, (mean_y - min_y) / span_y * (height - 1)) + 1, 1, height)
    if grid[c_row, c_col] == ' '
        grid[c_row, c_col] = '+'
    end

    border = "+" * repeat("-", width) * "+"
    println("\n--- Spot-Diagramm (ASCII-Visualisierung) ---")
    println(border)
    for r in height:-1:1
        println("|" * String(grid[r, :]) * "|")
    end
    println(border)
    @printf("  Statistik (%d Strahlen):\n", n)
    @printf("    RMS Spot-Radius: %.3f µm (%.4f mm)\n", rms_radius * 1e6, rms_radius * 1e3)
    @printf("    Max Spot-Radius: %.3f µm (%.4f mm)\n", max_radius * 1e6, max_radius * 1e3)
    @printf("    Schwerpunkt (X, Y): (%.3f µm, %.3f µm)\n", mean_x * 1e6, mean_y * 1e6)
    @printf("    Spot-Ausdehnung (ΔX × ΔY): %.2f µm × %.2f µm\n", span_x * 1e6, span_y * 1e6)
    println()
end
