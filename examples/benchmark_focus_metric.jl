#!/usr/bin/env julia

"""
Benchmark Focus Metric & Airy Ratio Analysis
Evaluates across ZMX optical models:
  Kategorie / Modell, Nominal y, Optimal y, Defokus Δy, RMS (nominal),
  RMS (refokussiert), r_Airy, RMS_opt/r_Airy, Status (<= 5.0)

Generates:
  1. CSV export (examples/output/focus_metrics.csv)
  2. Summary statistics: Mean, Median, Skewness, Std Dev
  3. Quality distribution histogram (examples/output/focus_metrics_histogram.png)
  4. Vector graphics (PDF) for ray layout and spot diagram in each model's directory
"""

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using BeamletOptics
using BeamletOpticsPrescriptionIO
using LinearAlgebra
using Printf
using Statistics
using CairoMakie
using Logging

# Suppress glass warning spam during massive batch benchmarking
Logging.disable_logging(Logging.Warn)

# Parse command-line arguments
# Usage: julia --project=. examples/benchmark_focus_metric.jl [--max-per-cat=N] [--all] [--format=pdf|svg|both]
max_per_cat = 10 # Default: 10 models per category (~65 models across all 9 categories)
export_format = "pdf" # "pdf", "svg", or "both"

for arg in ARGS
    if startswith(arg, "--max-per-cat=")
        global max_per_cat = parse(Int, split(arg, "=")[2])
    elseif arg == "--all"
        global max_per_cat = typemax(Int)
    elseif startswith(arg, "--format=")
        global export_format = lowercase(split(arg, "=")[2])
    end
end

const base_dir = joinpath(@__DIR__, "..", "test_data", "dan_reiley")
if !isdir(base_dir)
    error("Test data directory not found: $base_dir")
end

const out_dir = joinpath(@__DIR__, "output")
mkpath(out_dir)
const csv_path = joinpath(out_dir, "focus_metrics.csv")
const plot_path = joinpath(out_dir, "focus_metrics_histogram.png")

struct ModelMetric
    category::String
    model_name::String
    y_nom::Float64       # [m]
    y_opt::Float64       # [m]
    delta_y::Float64     # [m]
    rms_nom::Float64     # [m]
    rms_opt::Float64     # [m]
    r_airy::Float64      # [m]
    ratio_opt::Float64   # dimensionless
    passed::Bool
    diagram_path::String
end

"""
    save_vector_ray_and_spot_diagram(res, source, foc, r_airy, outpath; title="", category="")

Generates a compact, pure 2D vector graphic (PDF / SVG) containing:
1. Meridian (Y-Z) ray tracing cross-section with exact curved lens profiles and rays.
2. 2D Spot diagram on the detector with centroid, RMS circle and GEO circle.
"""
function save_vector_ray_and_spot_diagram(res, source, foc, r_airy, outpath; title="", category="")
    fig = Figure(size=(1100, 480), fontsize=12)

    # 1. Left Panel: 2D Meridian (Y-Z) Optical Layout & Ray Paths
    ax1 = Axis(fig[1, 1],
               title="2D Strahlengang (Y-Z Schnitt) $(isempty(category) ? "" : "[$category]")",
               xlabel="Optische Achse Y [mm]", ylabel="Z [mm]")

    function draw_surface_arc!(ax, y_vert_m, R_m, diam_m; color=(:deepskyblue4, 0.9), lw=1.3)
        r_max = diam_m / 2.0 * 1e3
        rs = LinRange(-r_max, r_max, 50)
        y_v_mm = y_vert_m * 1e3
        ys = if isinf(R_m) || abs(R_m) > 1e5
            fill(y_v_mm, length(rs))
        else
            R_mm = R_m * 1e3
            y_v_mm .+ (R_mm .- sign(R_mm) .* sqrt.(max.(0.0, R_mm^2 .- rs.^2)))
        end
        lines!(ax, ys, rs, color=color, linewidth=lw)
        return ys
    end

    # Draw optical components
    for el in res.elements
        if el isa ZmxSinglet
            ys1 = draw_surface_arc!(ax1, el.axial_position, el.front_surface.radius, el.diameter)
            ys2 = draw_surface_arc!(ax1, el.axial_position + el.center_thickness, el.back_surface.radius, el.diameter)
            r_max = el.diameter / 2.0 * 1e3
            lines!(ax1, [ys1[1], ys2[1]], [-r_max, -r_max], color=(:deepskyblue4, 0.9), linewidth=1.0)
            lines!(ax1, [ys1[end], ys2[end]], [r_max, r_max], color=(:deepskyblue4, 0.9), linewidth=1.0)
        elseif el isa ZmxDoublet
            ys1 = draw_surface_arc!(ax1, el.axial_position, el.surface1.radius, el.diameter)
            ys2 = draw_surface_arc!(ax1, el.axial_position + el.thickness1, el.surface2.radius, el.diameter; color=(:gray50, 0.8), lw=1.0)
            ys3 = draw_surface_arc!(ax1, el.axial_position + el.thickness1 + el.thickness2, el.surface3.radius, el.diameter)
            r_max = el.diameter / 2.0 * 1e3
            lines!(ax1, [ys1[1], ys3[1]], [-r_max, -r_max], color=(:deepskyblue4, 0.9), linewidth=1.0)
            lines!(ax1, [ys1[end], ys3[end]], [r_max, r_max], color=(:deepskyblue4, 0.9), linewidth=1.0)
        elseif el isa ZmxTriplet
            ys1 = draw_surface_arc!(ax1, el.axial_position, el.surface1.radius, el.diameter)
            ys4 = draw_surface_arc!(ax1, el.axial_position + el.thickness1 + el.thickness2 + el.thickness3, el.surface4.radius, el.diameter)
            r_max = el.diameter / 2.0 * 1e3
            lines!(ax1, [ys1[1], ys4[1]], [-r_max, -r_max], color=(:deepskyblue4, 0.9), linewidth=1.0)
            lines!(ax1, [ys1[end], ys4[end]], [r_max, r_max], color=(:deepskyblue4, 0.9), linewidth=1.0)
        elseif el isa ZmxMirror
            draw_surface_arc!(ax1, el.axial_position, el.surface.radius, el.diameter; color=:firebrick, lw=2.5)
        elseif el isa ZmxStop
            r = el.diameter / 2.0 * 1e3
            y = el.axial_position * 1e3
            lines!(ax1, [y, y], [r, r * 1.3], color=:black, linewidth=2.5)
            lines!(ax1, [y, y], [-r * 1.3, -r], color=:black, linewidth=2.5)
        elseif el isa ZmxDetector
            r = el.diameter / 2.0 * 1e3
            y = el.axial_position * 1e3
            lines!(ax1, [y, y], [-r, r], color=:firebrick, linewidth=2.5, label="Detektor")
        end
    end

    # Draw ray paths
    for b in source.beams
        for r in b.rays
            p1 = position(r)
            inter = BeamletOptics.intersection(r)
            len = inter === nothing ? 0.05 : length(inter)
            p2 = p1 + len * BeamletOptics.direction(r)
            lines!(ax1, [p1[2] * 1e3, p2[2] * 1e3], [p1[3] * 1e3, p2[3] * 1e3],
                   color=(:royalblue, 0.6), linewidth=0.7)
        end
    end

    # 2. Right Panel: 2D Spot Diagram on Detector
    ax2 = Axis(fig[1, 2], aspect=DataAspect(),
               title="$title\n(RMS: $(round(foc.rms_nom * 1e6, digits=2)) µm, Opt: $(round(foc.rms_opt * 1e6, digits=2)) µm, Airy: $(round(r_airy * 1e6, digits=2)) µm)",
               xlabel="x [µm]", ylabel="z [µm]")

    spots = spot_diagram(res.detector)
    if !isempty(spots)
        xs_um = [p[1] * 1e6 for p in spots]
        zs_um = [p[2] * 1e6 for p in spots]
        scatter!(ax2, xs_um, zs_um, markersize=4.5, color=:royalblue, label="$(length(spots)) Strahlen")

        # Centroid
        cx_um = sum(xs_um) / length(xs_um)
        cz_um = sum(zs_um) / length(zs_um)
        scatter!(ax2, [cx_um], [cz_um], marker=:cross, markersize=12, color=:black, label="Schwerpunkt")

        # RMS & Airy Circles
        θ_circ = LinRange(0, 2π, 120)
        rms_um = foc.rms_nom * 1e6
        airy_um = r_airy * 1e6
        lines!(ax2, cx_um .+ rms_um .* cos.(θ_circ), cz_um .+ rms_um .* sin.(θ_circ),
               color=:firebrick, linestyle=:dash, linewidth=1.5, label="RMS Radius")
        lines!(ax2, cx_um .+ airy_um .* cos.(θ_circ), cz_um .+ airy_um .* sin.(θ_circ),
               color=:forestgreen, linestyle=:dot, linewidth=1.5, label="Airy Radius")
        axislegend(ax2, position=:rt, labelsize=9)
    else
        text!(ax2, 0, 0, text="No rays intercepted the detector", align=(:center, :center))
    end

    mkpath(dirname(outpath))
    save(outpath, fig)
    return fig
end

println("="^135)
println("SYSTEMATIC FOCUS METRIC & AIRY RATIO ANALYSIS (INCL. VECTOR GRAPHICS)")
println("Models per category: $(max_per_cat == typemax(Int) ? "ALL" : string(max_per_cat)) | Format: $export_format")
println("="^135)

category_folders = Tuple{String, String}[]
for root in [joinpath(@__DIR__, "..", "test_data", "dan_reiley"), joinpath(@__DIR__, "..", "test_data", "sample")]
    if isdir(root)
        if basename(root) == "sample"
            push!(category_folders, ("sample", root))
        else
            for d in readdir(root)
                dp = joinpath(root, d)
                if isdir(dp)
                    push!(category_folders, (d, dp))
                end
            end
        end
    end
end
sort!(category_folders, by=first)

metrics = ModelMetric[]
t0 = time()
skipped_count = 0

for (cat, cat_dir) in category_folders
    files = filter(f -> endswith(lowercase(f), ".zmx"), readdir(cat_dir))
    sort!(files)
    selected_files = files[1:min(max_per_cat, length(files))]

    println("\n--> Analyzing category: $cat ($(length(selected_files)) / $(length(files)) files)")

    for (idx, f) in enumerate(selected_files)
        if idx % 25 == 0 || idx == length(selected_files)
            print("  [$cat] Progress: $idx / $(length(selected_files))\n")
            flush(stdout)
        end
        path = joinpath(cat_dir, f)
        try
            res = import_zmx(path)
            if res.detector === nothing
                global skipped_count += 1
                continue
            end

            # Generate rays according to ZMX prescription (finite/infinite aware)
            src = suggest_source(res; num_rays=60, num_rings=3)
            solve_system!(res.system, src)

            hits = res.detector.hits
            if hits === nothing || isempty(hits)
                global skipped_count += 1
                continue
            end

            # Analytical optimal focus calculation
            foc = find_best_focus(res.detector)
            y_nom = position(res.detector)[2]
            y_opt = foc.y_opt
            delta_y = foc.delta_y
            rms_nom = foc.rms_nom
            rms_opt = foc.rms_opt

            # Numerical aperture & Airy radius
            dirs = [h.ray.dir for h in hits]
            sin_thetas = [sqrt(d[1]^2 + d[3]^2) for d in dirs]
            na = maximum(sin_thetas)
            if na <= 1e-6
                na = 0.01
            end

            λ = hits[1].ray.λ
            if λ <= 0.0
                λ = res.zmx_system.wavelengths[res.zmx_system.primary_wavelength_idx]
            end

            r_airy = 0.61 * λ / na
            ratio_opt = rms_opt / r_airy
            passed = (ratio_opt <= 5.0)

            # Filter unphysical degenerate cases
            if !isnan(ratio_opt) && !isinf(ratio_opt) && ratio_opt >= 0.0
                # Generate and save vector graphics (PDF/SVG) in the model directory
                model_stem = splitext(f)[1]
                vector_path = ""
                if export_format in ("pdf", "both")
                    pdf_out = joinpath(cat_dir, "$(model_stem)_diagram.pdf")
                    try
                        save_vector_ray_and_spot_diagram(res, src, foc, r_airy, pdf_out; title=f, category=cat)
                        vector_path = pdf_out
                    catch
                    end
                end
                if export_format in ("svg", "both")
                    svg_out = joinpath(cat_dir, "$(model_stem)_diagram.svg")
                    try
                        save_vector_ray_and_spot_diagram(res, src, foc, r_airy, svg_out; title=f, category=cat)
                        vector_path = svg_out
                    catch
                    end
                end

                push!(metrics, ModelMetric(
                    cat,
                    f,
                    y_nom,
                    y_opt,
                    delta_y,
                    rms_nom,
                    rms_opt,
                    r_airy,
                    ratio_opt,
                    passed,
                    vector_path
                ))
            else
                global skipped_count += 1
            end
        catch e
            global skipped_count += 1
        end
    end
end

total_time = time() - t0
N = length(metrics)
println("\n" * "="^135)
println("ANALYSIS RESULTS: $N models successfully computed ($skipped_count skipped) in $(round(total_time, digits=1)) s")
println("="^135)

# 1. Print formatted console table
@printf("%-24s | %-22s | %10s | %10s | %11s | %10s | %10s | %8s | %12s | %s\n",
        "Category", "Model", "Nominal y", "Optimal y", "Defocus Δy", "RMS (nom)", "RMS (opt)", "r_Airy", "RMS_opt/Airy", "Status (<= 5.0)")
println("-"^135)

for m in metrics
    dy_str = abs(m.delta_y * 1e3) < 1.0 ? @sprintf("%+8.2f µm", m.delta_y * 1e6) : @sprintf("%+8.3f mm", m.delta_y * 1e3)
    status_str = m.passed ? "PASSED" : "ATTENTION"
    @printf("%-24s | %-22s | %8.3f mm | %8.3f mm | %11s | %8.2f µm | %8.2f µm | %6.2f µm | %12.2f | %s\n",
            m.category,
            length(m.model_name) > 22 ? m.model_name[1:19] * "..." : m.model_name,
            m.y_nom * 1e3,
            m.y_opt * 1e3,
            dy_str,
            m.rms_nom * 1e6,
            m.rms_opt * 1e6,
            m.r_airy * 1e6,
            m.ratio_opt,
            status_str)
end
println("="^135)

# 2. Write CSV export
open(csv_path, "w") do io
    println(io, "Category,Model,Nominal_y_mm,Optimal_y_mm,Defocus_dy_um,RMS_nom_um,RMS_opt_um,r_Airy_um,RMS_opt_over_r_Airy,Status_le_5,Vector_Graphic_Path")
    for m in metrics
        println(io, @sprintf("%s,%s,%.4f,%.4f,%.3f,%.3f,%.3f,%.3f,%.4f,%s,%s",
            m.category, m.model_name,
            m.y_nom * 1e3, m.y_opt * 1e3, m.delta_y * 1e6,
            m.rms_nom * 1e6, m.rms_opt * 1e6, m.r_airy * 1e6,
            m.ratio_opt, m.passed ? "PASSED" : "ATTENTION", m.diagram_path))
    end
end
println("\n[CSV export saved]: $csv_path")

# 3. Statistical Analysis
ratios = [m.ratio_opt for m in metrics]
passed_count = count(m -> m.passed, metrics)
passed_pct = (passed_count / N) * 100.0

mean_val = mean(ratios)
median_val = median(ratios)
std_val = std(ratios)
skew_val = if std_val > 0 && N > 2
    (N / ((N - 1) * (N - 2))) * sum(((r - mean_val) / std_val)^3 for r in ratios)
else
    0.0
end

println("\n" * "="^60)
println("STATISTICAL SUMMARY METRICS FOR RMS_opt / r_Airy")
println("="^60)
@printf("  Evaluated Models (N):           %d\n", N)
@printf("  Passing Criterion (<= 5.0):     %d (%.1f %%)\n", passed_count, passed_pct)
@printf("  Mean (μ):                       %.3f\n", mean_val)
@printf("  Median:                         %.3f\n", median_val)
@printf("  Standard Deviation (σ):         %.3f\n", std_val)
@printf("  Skewness (γ₁):                  %.3f\n", skew_val)
println("="^70)

# 4. Plot Histogram
fig = Figure(size=(1200, 700), fontsize=13)

ax1 = Axis(fig[1, 1],
           title="Linear Distribution: RMS_opt / r_Airy (Range 0 - 25)",
           xlabel="Ratio RMS_opt / r_Airy",
           ylabel="Number of Models")

ax2 = Axis(fig[1, 2],
           title="Full Distribution (Log10 Scale)",
           xlabel="log10(RMS_opt / r_Airy)",
           ylabel="Number of Models")

ratios_capped = filter(r -> r <= 30.0, ratios)
hist!(ax1, ratios_capped, bins=25, color=(:royalblue, 0.65), strokewidth=1.2, strokecolor=:navy)
vlines!(ax1, [5.0], color=:firebrick, linestyle=:dash, linewidth=2.0, label="Criterion (<= 5.0)")
vlines!(ax1, [median_val], color=:forestgreen, linestyle=:dot, linewidth=2.0, label=@sprintf("Median: %.2f", median_val))
vlines!(ax1, [min(25.0, mean_val)], color=:darkorange, linestyle=:dashdot, linewidth=2.0, label=@sprintf("Mean: %.2f", mean_val))
axislegend(ax1, position=:rt)

log_ratios = [log10(max(1e-3, r)) for r in ratios]
hist!(ax2, log_ratios, bins=25, color=(:seagreen, 0.65), strokewidth=1.2, strokecolor=:darkgreen)
vlines!(ax2, [log10(5.0)], color=:firebrick, linestyle=:dash, linewidth=2.0, label="log10(5.0) = 0.70")
vlines!(ax2, [log10(max(1e-3, median_val))], color=:forestgreen, linestyle=:dot, linewidth=2.0, label=@sprintf("log10(Med): %.2f", log10(max(1e-3, median_val))))
vlines!(ax2, [log10(max(1e-3, mean_val))], color=:darkorange, linestyle=:dashdot, linewidth=2.0, label=@sprintf("log10(Mean): %.2f", log10(max(1e-3, mean_val))))
axislegend(ax2, position=:rt)

stats_text = @sprintf("Sample: N = %d  |  Pass Rate (<= 5.0): %.1f%%  |  Mean: %.2f  |  Median: %.2f  |  Std: %.2f  |  Skew: %.2f",
                      N, passed_pct, mean_val, median_val, std_val, skew_val)
Label(fig[0, 1:2], "Quality Benchmark: Ratio of RMS Spot to Airy Radius after Refocusing\n$stats_text",
      font=:bold, fontsize=14, color=:gray15)

save(plot_path, fig, px_per_unit=2)
println("\n[Histogram saved]: $plot_path")

# Copy to artifacts directory
const artifact_plot = "/home/agent/.gemini/antigravity/brain/f1e8176e-4884-456e-add4-ba811deace17/focus_metrics_histogram.png"
try
    cp(plot_path, artifact_plot, force=true)
    println("[Histogram copied to artifact directory]: $artifact_plot")
catch
end
