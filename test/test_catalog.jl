using Test
using BeamletOptics
using BeamletOpticsPrescriptionIO
using LinearAlgebra

@testset "Catalog Integration (Thorlabs)" begin
    # 1. Listing lenses
    thor_lenses = list_catalog_lenses(:thorlabs)
    @test length(thor_lenses) >= 4000
    @test "AC254-050-A" in thor_lenses
    @test "LA1951" in thor_lenses

    # 2. Searching lenses
    ac_matches = find_catalog_lenses(:thorlabs, r"^AC254-050")
    @test length(ac_matches) >= 5
    @test "AC254-050-A" in ac_matches

    sub_matches = find_catalog_lenses(:thorlabs, "la1951")
    @test "LA1951" in sub_matches

    # 3. Loading Singlet
    singlet = load_lens_from_zmx_cat(:thorlabs, "LA1951")
    @test singlet isa BeamletOptics.Lens
    @test singlet isa BeamletOptics.AbstractObject

    # 4. Loading Doublet
    doublet = load_lens_from_zmx_cat(:thorlabs, "AC254-050-A")
    @test doublet isa BeamletOptics.DoubletLens
    @test doublet isa BeamletOptics.AbstractObject

    # 5. Case-insensitivity & optional .zmx extension
    doublet_lower = load_lens_from_zmx_cat(:thorlabs, "ac254-050-a.zmx")
    @test doublet_lower isa BeamletOptics.DoubletLens

    # 6. Positioning
    # Front vertex translation along optical axis (+Y)
    lens_pos_scalar = load_lens_from_zmx_cat(:thorlabs, "AC254-050-A", position=0.25)
    @test lens_pos_scalar isa BeamletOptics.DoubletLens

    lens_pos_vec = load_lens_from_zmx_cat(:thorlabs, "AC254-050-A", position=[0.01, 0.25, 0.0])
    @test lens_pos_vec isa BeamletOptics.DoubletLens

    # 7. return_result keyword
    (lens_res, res) = load_lens_from_zmx_cat(:thorlabs, "AC254-050-A", return_result=true)
    @test lens_res isa BeamletOptics.DoubletLens
    @test res isa BMOImportResult
    @test length(res.elements) >= 1
    @test res.detector isa BeamletOptics.Detector

    # 8. Error handling
    @test_throws ArgumentError load_lens_from_zmx_cat(:thorlabs, "NON_EXISTENT_LENS_XYZ123")
    @test_throws ArgumentError load_lens_from_zmx_cat(:invalid_catalog_sym, "LA1951")

    # 9. Direct Ray Tracing with loaded catalog lens
    # Construct a custom mini system with the loaded lens and a detector
    sys = System([doublet, Detector(0.025)])
    # AC254-050-A has focal length ~50 mm
    translate3d!(sys.objects[2], [0.0, 0.0547, 0.0])
    src = CollimatedSource([0.0, -0.01, 0.0], [0.0, 1.0, 0.0], 0.015, 587.56e-9; num_rings=3, num_rays=60)
    solve_system!(sys, src)
    hits = spot_diagram(sys.objects[2])
    @test length(hits) >= 40
end
