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
        
        # Surface 7 is STOP
        @test sys.surfaces[8].index == 7
        @test sys.surfaces[8].is_stop == true
    end

    @testset "Element Grouping on Double Gauss" begin
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
    end

    @testset "End-to-End Ray Tracing (Double Gauss)" begin
        path = joinpath(@__DIR__, "..", "test_data", "dan_reiley", "photographic-lenses-prime", "US00583336-2-scaled.zmx")
        res = import_zmx(path)
        
        @test res.system isa BeamletOptics.System
        @test res.detector isa BeamletOptics.Detector
        
        source = CollimatedSource([0.0, -0.01, 0.0], [0.0, 1.0, 0.0], 0.004, 587.56e-9, num_rays=200, num_rings=6)
        solve_system!(res.system, source)
        hits = spot_diagram(res.detector)
        
        # All 200 rays should pass through the 4 lens groups and focus onto the detector
        @test length(hits) == 200
        # Check that rays are focused (centroid spread < 5mm)
        @test all(norm(h) < 0.005 for h in hits)
    end

    @testset "Telescope with Mirror (Figure1.zmx)" begin
        path = joinpath(@__DIR__, "..", "test_data", "dan_reiley", "telescopes", "Figure1.zmx")
        if isfile(path)
            res = import_zmx(path)
            @test any(e isa ZmxMirror for e in res.elements)
            mirror = first(filter(e -> e isa ZmxMirror, res.elements))
            @test isapprox(mirror.surface.radius, -216.262mm, atol=1e-3)
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

end
