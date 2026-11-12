#!/usr/bin/env julia

"""
Cross-Platform Test Data Download & Setup Script for BeamletOpticsPrescriptionIO.jl
=======================================================================

Downloads and unpacks external reference optical design databases used for
benchmarking and validation:
  1. Dan Reiley Optical Patent Database (~980 models across 9 categories)
  2. Thorlabs Zemax Component Catalog (4,217 catalog lenses)

LEGAL DISCLAIMER:
  The external files fetched by this script are hosted and maintained by third parties.
  BeamletOpticsPrescriptionIO.jl does NOT own, redistribute, or license this third-party data.
  By running this script, you acknowledge that you are downloading publicly available
  reference data for personal, academic, or evaluation purposes, and you agree to comply
  with any applicable third-party terms of service, patent rights, and copyright laws.
"""

using Downloads
using Printf

const ROOT_DIR = normpath(joinpath(@__DIR__, ".."))
const TEST_DATA_DIR = joinpath(ROOT_DIR, "test_data")

# URLs
const DAN_REILEY_URL = "https://github.com/DanReiley/OpticalDesignNotebook/archive/refs/heads/master.zip"
const THORLABS_URL = "https://media.thorlabs.com/contentassets/0ee24bcce2224cc6ae338b2bdaa9e212/zemaxcatalog.zip"

function print_disclaimer(auto_yes::Bool)
    println("="^80)
    println("BeamletOpticsPrescriptionIO.jl - Reference Test Data Download Helper")
    println("="^80)
    println("""
LEGAL NOTICE & DISCLAIMER:
This script downloads external reference optical design files from third-party
sources for testing and validation purposes.

The authors of BeamletOpticsPrescriptionIO.jl do NOT own or host this data. You are solely
responsible for complying with the respective third-party terms of service and
applicable copyright laws.
""")
    println("="^80)
    if !auto_yes
        print("Do you wish to proceed with downloading the reference data? (y/N): ")
        ans = readline()
        if lowercase(strip(ans)) != "y"
            println("Download aborted by user.")
            exit(0)
        end
    end
end

function extract_zip(zip_path, dest_dir)
    try
        if Sys.iswindows()
            try
                run(`tar -xf $zip_path -C $dest_dir`)
            catch
                run(`powershell -Command "Expand-Archive -Path '$zip_path' -DestinationPath '$dest_dir' -Force"`)
            end
        else
            run(`unzip -q -o $zip_path -d $dest_dir`)
        end
    catch
        println("  [Notice]: System unzip failed. Falling back to Python zipfile...")
        try
            run(`python3 -c "import zipfile, sys; zipfile.ZipFile(sys.argv[1], 'r').extractall(sys.argv[2])" $zip_path $dest_dir`)
        catch e
            println("  [Error]: Failed to extract $zip_path automatically. Please extract it manually to $dest_dir.")
            rethrow(e)
        end
    end
end

function download_and_extract_dan_reiley()
    println("\n[1/2] Processing Dan Reiley Optical Patent Database...")
    target_dir = joinpath(TEST_DATA_DIR, "dan_reiley")
    if isdir(target_dir) && length(readdir(target_dir)) >= 9
        println("  -> Dan Reiley patent library already present at $(relpath(target_dir, ROOT_DIR)). Skipping download.")
        return
    end

    mkpath(TEST_DATA_DIR)
    zip_path = joinpath(TEST_DATA_DIR, "dan_reiley_master.zip")

    println("  Downloading Dan Reiley archive from GitHub...")
    try
        Downloads.download(DAN_REILEY_URL, zip_path)
        println("  Extracting archive...")
        # Cross-platform extraction using system unzip or tar
        extract_zip(zip_path, TEST_DATA_DIR)
        extracted_master = joinpath(TEST_DATA_DIR, "OpticalDesignNotebook-master")
        if isdir(extracted_master)
            # Reorganize into test_data/dan_reiley
            mkpath(target_dir)
            categories = ["endoscopes", "eyepieces", "microscope-objectives",
                          "photographic-lenses-prime", "photographic-lenses-zoom",
                          "projectors", "scan-lenses", "spectro", "telescopes"]
            for cat in categories
                src_cat = joinpath(extracted_master, cat)
                if isdir(src_cat)
                    dst_cat = joinpath(target_dir, cat)
                    rm(dst_cat, recursive=true, force=true)
                    cp(src_cat, dst_cat)
                end
            end
            rm(extracted_master, recursive=true, force=true)
        end
        rm(zip_path, force=true)
        println("  -> Successfully set up Dan Reiley test suite in $(relpath(target_dir, ROOT_DIR))")
    catch e
        println("  [Warning]: Failed to download from GitHub: $e")
        println("  You can manually clone https://github.com/DanReiley/OpticalDesignNotebook into test_data/dan_reiley")
    end
end

# De-obfuscation keystream for Thorlabs binary ZMF files
function decrypt_thorlabs_zmf(data::Vector{UInt8}, a::Real, b::Real)::Vector{UInt8}
    iv = cos(6.0 * a + 3.0 * b)
    iv = cos(655.0 * (pi / 180.0) * iv) + iv
    out = similar(data)
    for p in 0:(length(data)-1)
        val = 13.2 * (iv + sin(17.0 * (p + 3.0))) * (p + 1.0)
        s = @sprintf("%.8e", val)
        # Digits 5 to 7 in scientific string
        k = parse(Int, s[5:7]) & 0xFF
        out[p + 1] = data[p + 1] ⊻ UInt8(k)
    end
    return out
end

function extract_thorlabs_catalog()
    println("\n[2/2] Processing Thorlabs Zemax Component Catalog...")
    zmx_dir = joinpath(TEST_DATA_DIR, "thorlabs", "zmx")
    if isdir(zmx_dir) && length(readdir(zmx_dir)) >= 4000
        println("  -> Thorlabs catalog already unpacked ($(length(readdir(zmx_dir))) lenses in $(relpath(zmx_dir, ROOT_DIR))). Skipping.")
        return
    end

    mkpath(TEST_DATA_DIR)
    zip_path = joinpath(TEST_DATA_DIR, "zemaxcatalog.zip")

    # If file doesn't exist, try downloading from Thorlabs CDN
    if !isfile(zip_path)
        println("  Downloading Thorlabs catalog from media.thorlabs.com...")
        try
            Downloads.download(THORLABS_URL, zip_path)
        catch e
            println("  [Notice]: Automatic download from Thorlabs CDN failed ($e).")
            println("  Please manually download the Zemax Catalog from:")
            println("    https://www.thorlabs.com/software-pages/zemax")
            println("  and place 'zemaxcatalog.zip' into $(relpath(TEST_DATA_DIR, ROOT_DIR))/")
            return
        end
    end

    if isfile(zip_path)
        println("  Extracting Thorlabs catalog archive...")
        thor_base = joinpath(TEST_DATA_DIR, "thorlabs")
        mkpath(thor_base)
        extract_zip(zip_path, thor_base)

        # Locate .ZMF binary files
        mkpath(zmx_dir)
        zmf_files = String[]
        for (root, _, files) in walkdir(thor_base)
            for f in files
                if endswith(uppercase(f), ".ZMF")
                    push!(zmf_files, joinpath(root, f))
                end
            end
        end

        if !isempty(zmf_files)
            println("  De-obfuscating $(length(zmf_files)) Thorlabs ZMF binary archives...")
            for zmf_path in zmf_files
                data = read(zmf_path)
                # Parse Thorlabs ZMF layout
                # Bytes 0-3: format version (1001), then per lens:
                # name (100 bytes), 7 x Int32 (last = data length), keys a, b (2 x Float64), data
                if length(data) > 4
                    offset = 4
                    while offset + 144 <= length(data)
                        # Lens name: 100 bytes null-terminated
                        name_bytes = data[(offset + 1):(offset + 100)]
                        null_idx = findfirst(==(0x00), name_bytes)
                        fn = null_idx !== nothing ? String(name_bytes[1:(null_idx-1)]) : String(name_bytes)
                        fn = strip(fn) * ".zmx"
                        offset += 100
                        file_len = Int(reinterpret(Int32, data[(offset + 25):(offset + 28)])[1])
                        a, b = reinterpret(Float64, data[(offset + 29):(offset + 44)])
                        offset += 44
                        if offset + file_len <= length(data)
                            enc_content = data[(offset + 1):(offset + file_len)]
                            dec_content = decrypt_thorlabs_zmf(enc_content, a, b)
                            if endswith(lowercase(fn), ".zmx")
                                write(joinpath(zmx_dir, fn), dec_content)
                            end
                            offset += file_len
                        else
                            break
                        end
                    end
                end
            end
            println("  -> Successfully unpacked $(length(readdir(zmx_dir))) official Thorlabs .zmx files to $(relpath(zmx_dir, ROOT_DIR))")
        end
    end
end

function main()
    auto_yes = "-y" in ARGS || "--yes" in ARGS
    print_disclaimer(auto_yes)

    do_dan = "--dan-reiley" in ARGS || isempty(ARGS) || ("-y" in ARGS && length(ARGS)==1) || "--all" in ARGS
    do_thor = "--thorlabs" in ARGS || isempty(ARGS) || ("-y" in ARGS && length(ARGS)==1) || "--all" in ARGS

    if do_dan
        download_and_extract_dan_reiley()
    end
    if do_thor
        extract_thorlabs_catalog()
    end

    println("\n" * "="^80)
    println("Setup complete! You can now run the benchmarks:")
    println("  julia --project=. examples/benchmark_focus_metric.jl")
    println("  julia --project=. examples/benchmark_thorlabs.jl --category=achromats")
    println("="^80)
end

main()
