#!/usr/bin/env julia

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using BeamletOpticsPrescriptionIO
using Printf

function print_usage()
    println("""
Usage:
  julia zmx2bmo.jl <input.zmx> [output.jl] [options]

Arguments:
  input.zmx     Path to the Zemax .zmx file to import
  output.jl     (Optional) Path to save the generated BeamletOptics Julia script

Options:
  --info, -i    Print system metadata and element summary only
  --help, -h    Show this help message
""")
end

function main()
    args = ARGS
    if isempty(args) || "-h" in args || "--help" in args
        print_usage()
        return
    end

    info_only = "-i" in args || "--info" in args
    filtered_args = filter(a -> a != "-i" && a != "--info" && a != "-h" && a != "--help", args)

    if isempty(filtered_args)
        println(stderr, "Error: No input .zmx file specified.")
        print_usage()
        exit(1)
    end

    input_file = filtered_args[1]
    if !isfile(input_file)
        println(stderr, "Error: File not found: $input_file")
        exit(1)
    end

    zmx_sys = parse_zmx(input_file)
    elements = group_elements(zmx_sys)

    println("================================================================================")
    println("Zemax Import: $(isempty(zmx_sys.title) ? basename(input_file) : zmx_sys.title)")
    println("================================================================================")
    println("File units:   $(zmx_sys.unit) (scale to m: $(zmx_sys.scale_to_m))")
    println("Surfaces:     $(length(zmx_sys.surfaces))")
    println("Elements:     $(length(elements))")
    println("Wavelengths:  $([@sprintf("%.1fnm", w*1e9) for w in zmx_sys.wavelengths])")
    println("--------------------------------------------------------------------------------")
    println("Reconstructed Elements:")
    for (idx, el) in enumerate(elements)
        if isa(el, ZmxSinglet)
            r1 = isinf(el.front_surface.radius) ? "Inf" : @sprintf("%.2fmm", el.front_surface.radius*1e3)
            r2 = isinf(el.back_surface.radius) ? "Inf" : @sprintf("%.2fmm", el.back_surface.radius*1e3)
            println("  [$idx] Singlet: $(el.name) (Glass: $(el.glass_name), d=$(@sprintf("%.2fmm", el.center_thickness*1e3)), R1=$r1, R2=$r2, y=$(@sprintf("%.2fmm", el.axial_position*1e3)))")
        elseif isa(el, ZmxDoublet)
            println("  [$idx] Cemented Doublet: $(el.name) (Glasses: $(el.glass1) + $(el.glass2), y=$(@sprintf("%.2fmm", el.axial_position*1e3)))")
        elseif isa(el, ZmxTriplet)
            println("  [$idx] Cemented Triplet: $(el.name) (Glasses: $(el.glass1) + $(el.glass2) + $(el.glass3), y=$(@sprintf("%.2fmm", el.axial_position*1e3)))")
        elseif isa(el, ZmxMirror)
            r = isinf(el.surface.radius) ? "Plano" : @sprintf("%.2fmm", el.surface.radius*1e3)
            println("  [$idx] Mirror: $(el.name) (Radius: $r, diam=$(@sprintf("%.2fmm", el.diameter*1e3)), y=$(@sprintf("%.2fmm", el.axial_position*1e3)))")
        elseif isa(el, ZmxStop)
            println("  [$idx] Aperture Stop: y=$(@sprintf("%.2fmm", el.axial_position*1e3)), diam=$(@sprintf("%.2fmm", el.diameter*1e3))")
        elseif isa(el, ZmxDetector)
            println("  [$idx] Image Detector: y=$(@sprintf("%.2fmm", el.axial_position*1e3)), diam=$(@sprintf("%.2fmm", el.diameter*1e3))")
        end
    end
    println("================================================================================")

    if !info_only
        output_file = length(filtered_args) >= 2 ? filtered_args[2] : replace(input_file, r"\.zmx$"i => ".jl")
        if output_file == input_file
            output_file = input_file * ".jl"
        end
        generate_bmo_script(input_file, output_file)
        println("Generated executable BMO Julia script saved to: $output_file")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
