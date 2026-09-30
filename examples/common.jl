"""
Shared utilities for rendering and analyzing optical systems in the example scripts.
Utilizes BeamletOptics' native spot_diagram(detector) tool, statistical metrics,
and CairoMakie 3D/2D visualization.
"""

using Printf
using LinearAlgebra
using BeamletOptics
using CairoMakie

function display_system_summary(res, category::String)
    zsys = res.zmx_system
    println("="^80)
    println("CATEGORY:     $category")
    println("MODEL:        $(isempty(zsys.title) ? "Zemax Design" : zsys.title)")
    println("="^80)
    println("File Unit:    $(zsys.unit) (Scale to m: $(zsys.scale_to_m))")
    println("Surfaces:     $(length(zsys.surfaces))")
    println("Assemblies:   $(length(res.elements))")
    println("Wavelengths:  $([@sprintf("%.1f nm", w*1e9) for w in zsys.wavelengths])")
    println("-"^80)
    println("Reconstructed BMO Components:")
    for (i, el) in enumerate(res.elements)
        if isa(el, ZmxSinglet)
            r1 = isinf(el.front_surface.radius) ? "Plano" : @sprintf("%.2f mm", el.front_surface.radius*1e3)
            r2 = isinf(el.back_surface.radius) ? "Plano" : @sprintf("%.2f mm", el.back_surface.radius*1e3)
            println(@sprintf("  [%d] Singlet Lens:       %-12s (Glass: %-8s d=%6.2f mm, R1=%9s, R2=%9s, y=%6.2f mm)",
                i, el.name, el.glass_name, el.center_thickness*1e3, r1, r2, el.axial_position*1e3))
        elseif isa(el, ZmxDoublet)
            println(@sprintf("  [%d] Cemented Doublet:   %-12s (Glasses: %s + %s, y=%6.2f mm)",
                i, el.name, el.glass1, el.glass2, el.axial_position*1e3))
        elseif isa(el, ZmxTriplet)
            println(@sprintf("  [%d] Cemented Triplet:   %-12s (Glasses: %s + %s + %s, y=%6.2f mm)",
                i, el.name, el.glass1, el.glass2, el.glass3, el.axial_position*1e3))
        elseif isa(el, ZmxMirror)
            r = isinf(el.surface.radius) ? "Plano" : @sprintf("%.2f mm", el.surface.radius*1e3)
            println(@sprintf("  [%d] Mirror:             %-12s (Radius: %s, Diam=%6.2f mm, y=%6.2f mm)",
                i, el.name, r, el.diameter*1e3, el.axial_position*1e3))
        elseif isa(el, ZmxStop)
            println(@sprintf("  [%d] Aperture Stop:      y=%6.2f mm, Diameter=%6.2f mm",
                i, el.axial_position*1e3, el.diameter*1e3))
        elseif isa(el, ZmxDetector)
            println(@sprintf("  [%d] Image Sensor:       y=%6.2f mm, Diameter=%6.2f mm",
                i, el.axial_position*1e3, el.diameter*1e3))
        end
    end
    println("="^80)
end

"""
    analyze_bmo_spot_diagram(detector)

Uses BeamletOptics' native `spot_diagram(detector)` function to extract hit coordinates
and calculates optical spot metrics: Centroid, RMS radius, GEO (maximum) radius,
and spatial extent.
"""
function analyze_bmo_spot_diagram(detector)
    spots = spot_diagram(detector) # Vector{Point2{Float64}} in local detector coords (meters)
    n = length(spots)
    if n == 0
        println("  (BMO spot_diagram: No rays intercepted the detector)")
        return (spots=spots, n=0, rms=0.0, geo=0.0, centroid=(0.0, 0.0), dx=0.0, dz=0.0)
    end

    xs = [p[1] for p in spots]
    zs = [p[2] for p in spots]

    cx = sum(xs) / n
    cz = sum(zs) / n

    # RMS spot radius
    rms_r = sqrt(sum((x - cx)^2 + (z - cz)^2 for (x, z) in zip(xs, zs)) / n)
    # Geometric (max) spot radius from centroid
    geo_r = maximum(sqrt((x - cx)^2 + (z - cz)^2) for (x, z) in zip(xs, zs))

    dx = maximum(xs) - minimum(xs)
    dz = maximum(zs) - minimum(zs)

    println("\n--- BMO Spot Diagram Analysis (`spot_diagram(detector)`) ---")
    @printf("  Detected Rays:         %d\n", n)
    @printf("  RMS Spot Radius:       %.3f µm (%.4f mm)\n", rms_r * 1e6, rms_r * 1e3)
    @printf("  GEO Spot Radius (max): %.3f µm (%.4f mm)\n", geo_r * 1e6, geo_r * 1e3)
    @printf("  Centroid (X, Z):       (%.3f µm, %.3f µm)\n", cx * 1e6, cz * 1e6)
    @printf("  Span (ΔX × ΔZ):        %.2f µm × %.2f µm\n", dx * 1e6, dz * 1e6)
    println("-"^60)

    return (spots=spots, n=n, rms=rms_r, geo=geo_r, centroid=(cx, cz), dx=dx, dz=dz)
end

"""
    plot_bmo_spot_diagram(res, source, outpath; title="BMO Spot Diagram", category="")

Creates a comprehensive visualization combining:
1. 3D ray trace and optical components rendered with BeamletOptics Makie extension
2. 2D Spot diagram using hit points from `spot_diagram(detector)` with RMS & GEO circles
Saves the figure to `outpath`.
"""
function plot_bmo_spot_diagram(res, source, outpath; title="BMO Spot Diagram", category="")
    if res.detector === nothing
        @warn "No detector present in optical system, skipping spot diagram plot."
        return nothing
    end

    stats = analyze_bmo_spot_diagram(res.detector)
    spots = stats.spots

    fig = Figure(size=(1100, 520), fontsize=13)

    # 1. 3D Ray Trace (aligned side view: light travels horizontally from left to right)
    ax1 = Axis3(fig[1, 1], aspect=:data,
                azimuth=0.0, elevation=0.15,
                title="3D Ray Trace $(isempty(category) ? "" : "($category)")",
                xlabel="X [m]", ylabel="Y [m]", zlabel="Z [m]")
    
    # Render system components safely
    for obj in res.system.objects
        try
            render!(ax1, obj; alpha=0.35)
        catch
        end
    end

    total_y = res.detector !== nothing ? position(res.detector)[2] : 0.05
    render_every = max(1, div(length(source.beams), 25))
    render!(ax1, source; render_every=render_every, flen=max(0.001, total_y * 0.02), color=:royalblue)

    # 2. 2D Spot Diagram
    ax2 = Axis(fig[1, 2], aspect=DataAspect(),
               xlabel="x [µm]", ylabel="z [µm]",
               title="$title\n(RMS: $(round(stats.rms*1e6, digits=2)) µm, GEO: $(round(stats.geo*1e6, digits=2)) µm)")

    if !isempty(spots)
        xs_um = [p[1] * 1e6 for p in spots]
        zs_um = [p[2] * 1e6 for p in spots]
        scatter!(ax2, xs_um, zs_um, markersize=5, color=:royalblue, label="$(stats.n) Rays")

        # Plot centroid
        cx_um = stats.centroid[1] * 1e6
        cz_um = stats.centroid[2] * 1e6
        scatter!(ax2, [cx_um], [cz_um], marker=:cross, markersize=14, color=:black, label="Centroid")

        # RMS circle
        θ = LinRange(0, 2π, 150)
        rms_um = stats.rms * 1e6
        geo_um = stats.geo * 1e6
        lines!(ax2, cx_um .+ rms_um .* cos.(θ), cz_um .+ rms_um .* sin.(θ),
               color=:firebrick, linestyle=:dash, linewidth=1.5, label="RMS Radius")
        lines!(ax2, cx_um .+ geo_um .* cos.(θ), cz_um .+ geo_um .* sin.(θ),
               color=:gray50, linestyle=:dot, linewidth=1.2, label="GEO Radius")

        axislegend(ax2, position=:rt, labelsize=10)
    else
        text!(ax2, 0, 0, text="No rays intercepted the detector", align=(:center, :center))
    end

    mkpath(dirname(outpath))
    save(outpath, fig, px_per_unit=2)
    println("  -> Saved spot diagram plot to: $outpath")

    return fig
end

function print_ascii_spot_diagram(hits; width=60, height=18)
    if isempty(hits)
        println("  (No rays intercepted the detector)")
        return
    end
    
    xs = [h[1] for h in hits]
    ys = [h[2] for h in hits]
    
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
    
    c_col = clamp(round(Int, (mean_x - min_x) / span_x * (width - 1)) + 1, 1, width)
    c_row = clamp(round(Int, (mean_y - min_y) / span_y * (height - 1)) + 1, 1, height)
    if grid[c_row, c_col] == ' '
        grid[c_row, c_col] = '+'
    end

    border = "+" * repeat("-", width) * "+"
    println("\n--- Spot Diagram (ASCII Preview) ---")
    println(border)
    for r in height:-1:1
        println("|" * String(grid[r, :]) * "|")
    end
    println(border)
end

"""
    verify_airy_criterion(res, source; max_ratio=5.0)

Calculates the numerical aperture NA of the focused ray bundle on the detector:
  NA = sin(θ_marginal)
Calculates the theoretical diffraction-limited Airy radius:
  r_Airy = 0.61 * λ / NA
And compares it against the geometric RMS spot radius:
  ratio = r_RMS / r_Airy
Checks whether ratio <= max_ratio (default: 5.0).
Returns a named tuple `(passed=Bool, r_rms=Float64, r_airy=Float64, ratio=Float64, na=Float64)`.
"""
function verify_airy_criterion(res, source; max_ratio=5.0)
    if res.detector === nothing || res.detector.hits === nothing || isempty(res.detector.hits)
        @warn "No detector hits available for Airy criterion verification."
        return (passed=false, r_rms=NaN, r_airy=NaN, ratio=NaN, na=NaN)
    end

    spots = spot_diagram(res.detector)
    n = length(spots)
    if n == 0
        @warn "No detector hits available for Airy criterion verification."
        return (passed=false, r_rms=NaN, r_airy=NaN, ratio=NaN, na=NaN)
    end

    # 1. Spot Centroid and RMS in detector plane
    xs = [p[1] for p in spots]
    zs = [p[2] for p in spots]
    cx = sum(xs) / n
    cz = sum(zs) / n
    r_rms = sqrt(sum((x - cx)^2 + (z - cz)^2 for (x, z) in zip(xs, zs)) / n)

    # 2. Convergence angle and NA
    hits = res.detector.hits
    sin_thetas = [sqrt(h.ray.dir[1]^2 + h.ray.dir[3]^2) for h in hits]
    na = maximum(sin_thetas)
    if na <= 1e-6
        na = 0.01
    end

    # 3. Wavelength
    λ = hits[1].ray.λ
    if λ <= 0.0
        λ = res.zmx_system.wavelengths[res.zmx_system.primary_wavelength_idx]
    end

    # 4. Airy radius
    r_airy = 0.61 * λ / na
    ratio = r_rms / r_airy

    passed = (ratio <= max_ratio)

    # 5. Optimal Focus Analysis
    foc = find_best_focus(res.detector)
    ratio_opt = foc.rms_opt / r_airy
    passed_opt = (ratio_opt <= max_ratio)

    println("\n=== Airy Radius vs. RMS Spot Radius Verification ===")
    @printf("  Wavelength λ:            %.1f nm\n", λ * 1e9)
    @printf("  Numerical Aperture (NA): %.4f (Half-cone angle θ ≈ %.2f°)\n", na, rad2deg(asin(clamp(na, 0.0, 1.0))))
    @printf("  Theor. Airy Radius:      %.3f µm\n", r_airy * 1e6)
    @printf("  Nominal RMS Spot:        %.3f µm  (y = %.4f mm)\n", r_rms * 1e6, position(res.detector)[2] * 1e3)
    @printf("  Ratio RMS / Airy:        %.2f  (Criterion: <= %.1f) -> %s\n", 
            ratio, max_ratio, passed ? "PASSED" : "ATTENTION")
    @printf("  Optimal Focus (y_opt):   %.4f mm  (Δy = %+.2f µm)\n", 
            foc.y_opt * 1e3, foc.delta_y * 1e6)
    @printf("  Refocused RMS Spot:      %.3f µm  (RMS/Airy: %.2f) -> %s\n", 
            foc.rms_opt * 1e6, ratio_opt, passed_opt ? "PASSED (Near diffraction limit)" : "ATTENTION")
    println("="^60)

    return (passed=passed, r_rms=r_rms, r_airy=r_airy, ratio=ratio, na=na, 
            y_opt=foc.y_opt, delta_y=foc.delta_y, rms_opt=foc.rms_opt, ratio_opt=ratio_opt, passed_opt=passed_opt)
end

