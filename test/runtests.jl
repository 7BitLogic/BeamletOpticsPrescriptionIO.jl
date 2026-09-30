using Test
using BeamletOptics
using BeamletOpticsZMX
using LinearAlgebra

const mm = 1e-3

@testset "BeamletOpticsZMX Test Suite" begin

    @testset "Glass Catalog & Dispersion" begin
        # 1. Catalog glass lookup
        n_bk7 = resolve_refractive_index("N-BK7")
        @test isapprox(n_bk7(587.56e-9), 1.5168, atol=1e-3)
        @test n_bk7(486.13e-9) > n_bk7(656.27e-9) # Normal dispersion: blue index > red index

        n_sk16 = resolve_refractive_index("N-SK16")
        @test isapprox(n_sk16(587.56e-9), 1.6204, atol=1e-3)

        # 2. Air lookup
        n_air = resolve_refractive_index("AIR")
        @test n_air(550e-9) == 1.0

        # 3. Model glass via Cauchy dispersion
        nd_test = 1.6
        vd_test = 40.0
        n_cauchy = cauchy_dispersion(nd_test, vd_test)
        @test isapprox(n_cauchy(587.5618e-9), nd_test, atol=1e-5)
        # Check Abbe number calculation
        nd_calc = n_cauchy(587.5618e-9)
        nF_calc = n_cauchy(486.1327e-9)
        nC_calc = n_cauchy(656.2725e-9)
        vd_calc = (nd_calc - 1.0) / (nF_calc - nC_calc)
        @test isapprox(vd_calc, vd_test, atol=1e-3)
        # 4. Historical glass lookup
        n_basf1 = resolve_refractive_index("BASF1")
        @test isapprox(n_basf1(587.56e-9), 1.62606, atol=1e-3)
        @test isapprox((n_basf1(587.56e-9) - 1.0) / (n_basf1(486.13e-9) - n_basf1(656.27e-9)), 38.96, atol=1e-1)

        # 5. Model glass with Vd priority
        n_custom = resolve_refractive_index("CUSTOM", 1.55, 50.0)
        @test isapprox(n_custom(587.56e-9), 1.55, atol=1e-4)
        @test isapprox((n_custom(587.56e-9) - 1.0) / (n_custom(486.13e-9) - n_custom(656.27e-9)), 50.0, atol=1e-1)
    end

    @testset "Parser on Double Gauss" begin
        path = joinpath(@__DIR__, "..", "test_data", "dan_reiley", "photographic-lenses-prime", "US00583336-2-scaled.zmx")
        @test isfile(path)
        
        sys = parse_zmx(path)
        @test sys.title == "Original Double Gauss"
        @test sys.unit == :mm
        @test sys.scale_to_m == 1e-3
        @test length(sys.surfaces) == 14
        
        # Surface 0 is OBJ
        @test sys.surfaces[1].index == 0
        @test isinf(sys.surfaces[1].thickness)
        
        # Surface 2 has explicit Nd and Vd in GLAS line
        @test isapprox(sys.surfaces[3].nd, 1.572500, atol=1e-4)
        @test isapprox(sys.surfaces[3].vd, 57.549, atol=1e-2)

        # Surface 7 is STOP
        @test sys.surfaces[8].index == 7
        @test sys.surfaces[8].is_stop == true

        # Fields
        @test length(sys.fields) == 3
        @test sys.fields[1] == (0.0, 0.0, 1.0)
        @test isapprox(sys.fields[2][2], 17.0, atol=1e-1)
        @test isapprox(sys.fields[3][2], 24.0, atol=1e-1)
    end

    @testset "Element Grouping and Stop Preservation" begin
        path = joinpath(@__DIR__, "..", "test_data", "dan_reiley", "photographic-lenses-prime", "US00583336-2-scaled.zmx")
        sys = parse_zmx(path)
        elements = group_elements(sys)
        
        @test length(elements) == 6
        @test isa(elements[1], ZmxSinglet)
        @test isa(elements[2], ZmxDoublet)
        @test isa(elements[3], ZmxStop)
        @test isa(elements[4], ZmxDoublet)
        @test isa(elements[5], ZmxSinglet)
        @test isa(elements[6], ZmxDetector)

        # Check geometry of first singlet
        @test isapprox(elements[1].front_surface.radius, 14.161mm, atol=1e-5)
        @test isapprox(elements[1].center_thickness, 1.386mm, atol=1e-5)
        @test isapprox(elements[1].axial_position, 4.5mm, atol=1e-5)

        # Check detector position
        @test isapprox(elements[6].axial_position, 61.882mm, atol=1e-4)

        # Import system and check stop field
        res = import_zmx(path)
        @test res.stop !== nothing
        @test res.stop isa ZmxStop
        @test res.stop.diameter > 0.0
    end

    @testset "suggest_source & generate_zmx_rays" begin
        path = joinpath(@__DIR__, "..", "test_data", "dan_reiley", "photographic-lenses-prime", "US00583336-2-scaled.zmx")
        res = import_zmx(path)

        src = suggest_source(res; num_rays=50, num_rings=3)
        @test src isa BeamletOptics.CollimatedSource
        @test length(src.beams) > 0
        @test isapprox(src.center[2], -0.0055, atol=1e-3)

        zmx_sources = generate_zmx_rays(res; fields=1, wavelengths=:primary, num_rings=2, num_rays=20)
        @test length(zmx_sources) == 1
        @test length(zmx_sources[1].beams) > 0
    end

    @testset "End-to-End Ray Tracing (Double Gauss)" begin
        path = joinpath(@__DIR__, "..", "test_data", "dan_reiley", "photographic-lenses-prime", "US00583336-2-scaled.zmx")
        res = import_zmx(path)
        
        @test res.system isa BeamletOptics.System
        @test res.detector isa BeamletOptics.Detector
        
        source = suggest_source(res; num_rays=200, num_rings=6)
        solve_system!(res.system, source)
        hits = spot_diagram(res.detector)
        
        # Rays should pass through and focus onto the detector
        @test length(hits) >= 150
        # Focused spot RMS
        xs = [h[1] for h in hits]
        zs = [h[2] for h in hits]
        rms = sqrt(sum(x^2 + z^2 for (x, z) in zip(xs, zs)) / length(hits))
        @test rms < 0.001 # < 1 mm
    end

    @testset "Telescope with Mirror (Figure1.zmx)" begin
        path = joinpath(@__DIR__, "..", "test_data", "dan_reiley", "telescopes", "Figure1.zmx")
        if isfile(path)
            res = import_zmx(path)
            @test any(e isa ZmxMirror for e in res.elements)
            mirror = first(filter(e -> e isa ZmxMirror, res.elements))
            @test isapprox(mirror.surface.radius, -216.262mm, atol=1e-3)

            # Trace through telescope with central obscuration absorber
            src = suggest_source(res; num_rays=200, num_rings=6)
            solve_system!(res.system, src)
            hits = spot_diagram(res.detector)
            @test length(hits) > 0
            
            # Focused hits on detector should be tiny (< 50 um RMS!)
            xs = [h[1] for h in hits]
            zs = [h[2] for h in hits]
            rms = sqrt(sum(x^2 + z^2 for (x, z) in zip(xs, zs)) / length(hits))
            @test rms < 50e-6 # Sharp diffraction-limited focus!
        end
    end

    @testset "Code Generation" begin
        path = joinpath(@__DIR__, "..", "test_data", "dan_reiley", "photographic-lenses-prime", "US00583336-2-scaled.zmx")
        code = generate_bmo_code(path)
        
        @test occursin("using BeamletOptics", code)
        @test occursin("SphericalLens", code)
        @test occursin("SphericalDoubletLens", code)
        @test occursin("Detector", code)
        @test occursin("solve_system!", code)
        
        tmp_script = tempname() * ".jl"
        generate_bmo_script(path, tmp_script)
        @test isfile(tmp_script)
        
        # Execute the generated script in isolation
        m = Module(:GeneratedTest)
        Base.include(m, tmp_script)
        @test isdefined(m, :system)
        @test isdefined(m, :detector)
        rm(tmp_script, force=true)
    end

    @testset "Batch Robustness Test across Categories" begin
        base_dir = joinpath(@__DIR__, "..", "test_data", "dan_reiley")
        if isdir(base_dir)
            categories = filter(d -> isdir(joinpath(base_dir, d)), readdir(base_dir))
            tested_count = 0
            for cat in categories
                files = filter(f -> endswith(lowercase(f), ".zmx"), readdir(joinpath(base_dir, cat)))
                sample_files = files[1:min(5, length(files))]
                for f in sample_files
                    filepath = joinpath(base_dir, cat, f)
                    sys = parse_zmx(filepath)
                    @test length(sys.surfaces) >= 2
                    elements = group_elements(sys)
                    @test length(elements) >= 1
                    tested_count += 1
                end
            end
            println("Batch tested $tested_count ZMX files across $(length(categories)) categories without any errors!")
            @test tested_count >= 20
        end
    end

    @testset "Focus Optimization & Refocusing" begin
        dg_file = joinpath(@__DIR__, "..", "test_data", "dan_reiley", "photographic-lenses-prime", "US00583336-2-scaled.zmx")
        res = import_zmx(dg_file)
        src = suggest_source(res; num_rays=60, num_rings=3)
        solve_system!(res.system, src)
        
        opt = find_best_focus(res.detector)
        @test isapprox(opt.y_opt, 62.0e-3, atol=1.0e-3)
        @test opt.rms_opt <= opt.rms_nom
        @test opt.rms_opt < 25e-6 # < 25 um for full 10 mm pupil
        
        # Test refocus! function
        init_y = position(res.detector)[2]
        applied_dy = refocus!(res)
        @test isapprox(applied_dy, opt.delta_y, atol=1e-8)
        new_y = position(res.detector)[2]
        @test isapprox(new_y, init_y + applied_dy, atol=1e-8)
    end

    @testset "Finite vs. Infinite Conjugate Ray Generation" begin
        # 1. Projector: 3m finite object distance & height field
        proj_file = joinpath(@__DIR__, "..", "test_data", "dan_reiley", "projectors", "US03126786-1.zmx")
        res_proj = import_zmx(proj_file)
        @test isapprox(res_proj.zmx_system.object_distance, 3.0, atol=1e-4)
        @test res_proj.zmx_system.field_type == :real_image_height
        src_proj = suggest_source(res_proj)
        @test src_proj isa PointSource
        rays_proj = generate_zmx_rays(res_proj; num_rings=2, num_rays=40)
        @test length(rays_proj) >= 1
        @test rays_proj[1] isa PointSource

        # 2. Spectrometer: OBNA 0.41 & 24mm slit distance
        spec_file = joinpath(@__DIR__, "..", "test_data", "dan_reiley", "spectro", "US05644396-1.ZMX")
        res_spec = import_zmx(spec_file)
        @test isapprox(res_spec.zmx_system.obna, 0.41, atol=1e-4)
        @test isapprox(res_spec.zmx_system.object_distance, 0.024042, atol=1e-5)
        src_spec = suggest_source(res_spec)
        @test src_spec isa PointSource

        # 3. Double Gauss: Infinite conjugate
        dg_file = joinpath(@__DIR__, "..", "test_data", "dan_reiley", "photographic-lenses-prime", "US00583336-2-scaled.zmx")
        res_dg = import_zmx(dg_file)
        @test isinf(res_dg.zmx_system.object_distance)
        src_dg = suggest_source(res_dg)
        @test src_dg isa CollimatedSource
    end

    include("test_catalog.jl")

end
