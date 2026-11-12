#!/usr/bin/env julia

"""
Benchmark Focus Metric & Airy Ratio Analysis for Thorlabs Catalog Lenses
Evaluates across Thorlabs catalog lenses:
  Category / Model, Nominal y, Optimal y, Defocus Δy, RMS (nominal),
  RMS (refocused), r_Airy, RMS_opt/r_Airy, Status (<= 5.0)

Generates:
  1. CSV export (examples/output/thorlabs_focus_metrics.csv)
  2. Summary statistics & breakdown per lens category (Achromats, Aspheres, Singlets)
  3. Quality distribution histogram (examples/output/thorlabs_focus_histogram.png)
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
# Usage: julia --project=. examples/benchmark_thorlabs.jl [--category=all|achromats|aspheres|singlets] [--limit=N]
category_filter = "all"
max_limit = typemax(Int)

for arg in ARGS
    if startswith(arg, "--category=")
        global category_filter = lowercase(split(arg, "=")[2])
    elseif startswith(arg, "--limit=")
        global max_limit = parse(Int, split(arg, "=")[2])
    end
end

const thorlabs_dir = normpath(joinpath(@__DIR__, "..", "test_data", "thorlabs", "zmx"))
if !isdir(thorlabs_dir)
    error("Thorlabs catalog directory not found: $thorlabs_dir")
end

const out_dir = joinpath(@__DIR__, "output")
mkpath(out_dir)
const csv_path = joinpath(out_dir, "thorlabs_focus_metrics.csv")
const plot_path = joinpath(out_dir, "thorlabs_focus_histogram.png")

struct ThorlabsMetric
    category::String
    model_name::String
    focal_length_est::Float64 # [mm]
    y_nom::Float64            # [m]
    y_opt::Float64            # [m]
    delta_y::Float64          # [m]
    rms_nom::Float64          # [m]
    rms_opt::Float64          # [m]
    r_airy::Float64           # [m]
    ratio_opt::Float64        # dimensionless
    passed::Bool
end

# Categorize Thorlabs lens by prefix
function categorize_thorlabs_file(fname::String)
    upper = uppercase(fname)
    if startswith(upper, "AC") || startswith(upper, "PAC")
        return "Achromats (Doublets/Triplets)"
    elseif startswith(upper, "AL") || startswith(upper, "354") || startswith(upper, "C") || startswith(upper, "A")
        return "Aspheric Lenses"
    elseif startswith(upper, "LA")
        return "Plano-Convex Singlets (LA)"
    elseif startswith(upper, "LB")
        return "Bi-Convex Singlets (LB)"
    elseif startswith(upper, "LF")
        return "Positive Meniscus (LF)"
    else
        return "Other Focusing Lenses"
    end
end

# Check if category matches user filter
function matches_category(cat_name::String, filter_mode::String)
    if filter_mode == "all"
        return true
    elseif filter_mode == "achromats"
        return occursin("Achromat", cat_name)
    elseif filter_mode == "aspheres"
        return occursin("Aspher", cat_name)
    elseif filter_mode == "singlets"
        return occursin("Singlet", cat_name) || occursin("Meniscus", cat_name)
    else
        return occursin(lowercase(filter_mode), lowercase(cat_name))
    end
end

println("="^110)
println("THORLABS LENS CATALOG: SPOT RADIUS & DIFFRACTION LIMIT BENCHMARK")
println("Category Filter: $category_filter  |  Limit: $(max_limit == typemax(Int) ? "No Limit" : max_limit)")
println("="^110)

all_files = filter(f -> endswith(lowercase(f), ".zmx"), readdir(thorlabs_dir))
println("Total ZMX files in catalog: $(length(all_files))")

# Filter files that match category and exclude pure concave/cylindrical/flat types
candidate_files = Tuple{String, String}[]
for f in all_files
    # Exclude concave/diverging lenses, axicons, and pure flats
    # LC = plano-concave, LK = bi-concave, LD = concave, LE = negative meniscus
    # LJ = cylindrical, AX = axicon, WG/WW = windows
    upper = uppercase(f)
    if startswith(upper, "LC") || startswith(upper, "LK") || startswith(upper, "LD") ||
       startswith(upper, "LE") || startswith(upper, "LJ") || startswith(upper, "AX") ||
       startswith(upper, "WG") || startswith(upper, "WW")
        continue
    end

    cat = categorize_thorlabs_file(f)
    if matches_category(cat, category_filter)
        push!(candidate_files, (cat, f))
    end
end

println("Candidate focusing lenses matching filter '$category_filter': $(length(candidate_files))")

metrics = ThorlabsMetric[]
skipped_count = 0
t0 = time()

# Limit candidate files if requested
eval_files = length(candidate_files) > max_limit ? candidate_files[1:max_limit] : candidate_files

for (idx, (cat, f)) in enumerate(eval_files)
    filepath = joinpath(thorlabs_dir, f)
    model_name = f[1:end-4]
    
    try
        res = import_zmx(filepath)
        
        # Verify detector is downstream (+Y)
        if res.detector === nothing
            global skipped_count += 1
            continue
        end
        det_y = position(res.detector)[2]
        if det_y <= 0.001
            # Virtual focus (diverging / negative)
            global skipped_count += 1
            continue
        end

        src = suggest_source(res; num_rays=80, num_rings=4)
        solve_system!(res.system, src)
        
        hits = res.detector.hits
        if hits === nothing || length(hits) < 10
            global skipped_count += 1
            continue
        end

        opt = find_best_focus(res.detector)
        λ = res.zmx_system.wavelengths[res.zmx_system.primary_wavelength_idx]
        dirs = [h.ray.dir for h in hits]
        sin_thetas = [sqrt(d[1]^2 + d[3]^2) for d in dirs]
        na = maximum(sin_thetas)
        if na <= 1e-6
            na = 0.005
        end

        r_airy = 0.61 * λ / na
        ratio_opt = opt.rms_opt / r_airy
        passed = ratio_opt <= 5.0

        push!(metrics, ThorlabsMetric(
            cat,
            model_name,
            det_y * 1e3,
            det_y,
            opt.y_opt,
            opt.delta_y,
            opt.rms_nom,
            opt.rms_opt,
            r_airy,
            ratio_opt,
            passed
        ))

        if idx % 100 == 0 || idx == length(eval_files)
            @printf("  [%4d/%4d] Processed: %-22s -> RMS_opt/r_Airy = %6.2f [%s]\n",
                    idx, length(eval_files), length(model_name) > 22 ? model_name[1:22] : model_name,
                    ratio_opt, passed ? "PASS" : "FAIL")
        end
    catch e
        global skipped_count += 1
    end
end

total_time = time() - t0
N = length(metrics)
println("\n" * "="^110)
println("COMPLETED: $N lenses successfully evaluated ($skipped_count skipped/virtual) in $(round(total_time, digits=1)) s")
println("="^110)

if N == 0
    println("No lenses could be evaluated.")
    exit(0)
end

# 1. Print formatted console table of sample models
println("\nSample Results (First 25 models):")
@printf("%-26s | %-18s | %10s | %10s | %11s | %10s | %10s | %8s | %12s | %s\n",
        "Category", "Model", "Nominal y", "Optimal y", "Defocus Δy", "RMS (nom)", "RMS (opt)", "r_Airy", "RMS_opt/Airy", "Status (<= 5.0)")
println("-"^135)

for m in metrics[1:min(25, N)]
    dy_str = abs(m.delta_y * 1e3) < 1.0 ? @sprintf("%+8.2f µm", m.delta_y * 1e6) : @sprintf("%+8.3f mm", m.delta_y * 1e3)
    status_str = m.passed ? "PASSED" : "ATTENTION"
    @printf("%-26s | %-18s | %8.2f mm | %8.2f mm | %11s | %8.2f µm | %8.2f µm | %6.2f µm | %12.2f | %s\n",
            m.category,
            length(m.model_name) > 18 ? m.model_name[1:15] * "..." : m.model_name,
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
    println(io, "Category,Model,Focal_Length_Est_mm,Nominal_y_mm,Optimal_y_mm,Defocus_dy_um,RMS_nom_um,RMS_opt_um,r_Airy_um,RMS_opt_over_r_Airy,Status_le_5")
    for m in metrics
        println(io, @sprintf("%s,%s,%.3f,%.4f,%.4f,%.3f,%.3f,%.3f,%.3f,%.4f,%s",
            m.category, m.model_name, m.focal_length_est,
            m.y_nom * 1e3, m.y_opt * 1e3, m.delta_y * 1e6,
            m.rms_nom * 1e6, m.rms_opt * 1e6, m.r_airy * 1e6,
            m.ratio_opt, m.passed ? "PASSED" : "ATTENTION"))
    end
end
println("\n[CSV export saved]: $csv_path")

# 3. Category Breakdown & Statistics
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

println("\n" * "="^70)
println("STATISTICAL SUMMARY METRICS (RMS_opt / r_Airy)")
println("="^70)
@printf("  Evaluated Models (N):           %d\n", N)
@printf("  Passing Criterion (<= 5.0):     %d (%.1f %%)\n", passed_count, passed_pct)
@printf("  Mean (μ):                       %.3f\n", mean_val)
@printf("  Median:                         %.3f\n", median_val)
@printf("  Standard Deviation (σ):         %.3f\n", std_val)
@printf("  Skewness (γ₁):                  %.3f\n", skew_val)
println("="^70)

# Per-category summary table
cats = unique([m.category for m in metrics])
println("\nBreakdown by Lens Type:")
@printf("%-30s | %6s | %10s | %10s | %10s\n",
        "Lens Type", "Count", "Pass (<=5)", "Median", "Mean")
println("-"^75)
for c in cats
    c_metrics = filter(m -> m.category == c, metrics)
    cN = length(c_metrics)
    c_passed = count(m -> m.passed, c_metrics)
    c_pct = (c_passed / cN) * 100.0
    c_ratios = [m.ratio_opt for m in c_metrics]
    @printf("%-30s | %6d | %9.1f%% | %10.2f | %10.2f\n",
            c, cN, c_pct, median(c_ratios), mean(c_ratios))
end
println("="^75)

# 4. Generate Histogram Plot
fig = Figure(size=(1200, 700), fontsize=13)

ax1 = Axis(fig[1, 1],
           title="Linear Distribution: RMS_opt / r_Airy (Range 0 - 25)",
           xlabel="Ratio RMS_opt / r_Airy",
           ylabel="Number of Lenses")

ax2 = Axis(fig[1, 2],
           title="Full Distribution (Log10 Scale)",
           xlabel="log10(RMS_opt / r_Airy)",
           ylabel="Number of Lenses")

ratios_capped = filter(r -> r <= 30.0, ratios)
hist!(ax1, ratios_capped, bins=30, color=(:royalblue, 0.65), strokewidth=1.2, strokecolor=:navy)
vlines!(ax1, [5.0], color=:firebrick, linestyle=:dash, linewidth=2.0, label="Criterion (<= 5.0)")
vlines!(ax1, [min(30.0, median_val)], color=:forestgreen, linestyle=:dot, linewidth=2.0, label=@sprintf("Median: %.2f", median_val))
axislegend(ax1, position=:rt)

log_ratios = [log10(max(1e-3, r)) for r in ratios]
hist!(ax2, log_ratios, bins=30, color=(:seagreen, 0.65), strokewidth=1.2, strokecolor=:darkgreen)
vlines!(ax2, [log10(5.0)], color=:firebrick, linestyle=:dash, linewidth=2.0, label="log10(5.0) = 0.70")
vlines!(ax2, [log10(max(1e-3, median_val))], color=:forestgreen, linestyle=:dot, linewidth=2.0, label=@sprintf("log10(Med): %.2f", log10(max(1e-3, median_val))))
axislegend(ax2, position=:rt)

stats_text = @sprintf("Thorlabs Catalog (N = %d)  |  Pass Rate (<= 5.0): %.1f%%  |  Median: %.2f  |  Mean: %.2f",
                      N, passed_pct, median_val, mean_val)
Label(fig[0, 1:2], "Thorlabs Lens Catalog: Diffraction Limit Benchmark (RMS_opt / r_Airy)\n$stats_text",
      font=:bold, fontsize=14, color=:gray15)

save(plot_path, fig, px_per_unit=2)
println("\n[Histogram saved]: $plot_path")

# Copy to artifacts directory
const artifact_plot = "/home/agent/.gemini/antigravity/brain/f1e8176e-4884-456e-add4-ba811deace17/thorlabs_focus_histogram.png"
try
    cp(plot_path, artifact_plot, force=true)
    println("[Histogram copied to artifact storage]: $artifact_plot")
catch
end
