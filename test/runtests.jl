using Test
using MWFRSRoofPressures
using MWFRSRoofPressures: windward_neg_cp, windward_pos_cp, leeward_cp, zone_cp_a, area_reduction_factor,
                          velocity_pressure, GCPI_ENCLOSED, ZONE_CP_B

@testset "MWFRSRoofPressures" begin

base = (V=100.0, exposure="B", h=21.6, long_dim=30.0, short_dim=25.0, design_code="LRFD")
cases(r) = Dict(r.cases)

@testset "Fig. 27.3-1 table lookups" begin
    @test windward_neg_cp(10.0, 0.25) ≈ -0.7
    @test windward_neg_cp(20.0, 1.0) ≈ -0.7
    @test windward_neg_cp(12.5, 0.25) ≈ -0.6                 # linear in theta
    @test windward_neg_cp(10.0, 0.75) ≈ (-0.9 + -1.3) / 2    # linear in h/L
    @test windward_pos_cp(30.0, 0.5) ≈ 0.2
    @test windward_pos_cp(70.0, 0.5) ≈ 0.7                   # 0.01 theta above 60 deg
    @test leeward_cp(15.0, 0.25) ≈ -0.5
    @test zone_cp_a(0.5) == (-0.9, -0.9, -0.5, -0.3)
    @test area_reduction_factor(100.0) == 1.0 && area_reduction_factor(1000.0) == 0.8
end

@testset "input validation" begin
    @test_throws ErrorException compute_roof_pressures(; base..., roof_slope_rise_per_12=3.0, design_code="WSD")
    @test_throws ErrorException compute_roof_pressures(; base..., roof_slope_rise_per_12=3.0, ridge_parallel_to="x")
    @test_throws ErrorException compute_roof_pressures(; base..., roof_slope_rise_per_12=3.0, roof_type=:monoslope)
end

@testset "gable roofs" begin
    o = compute_roof_pressures(; base..., roof_slope_rise_per_12=3.0)
    @test o.roof_type == :gable
    d1, d2 = o.results
    @test d1.normal_to_ridge && !d1.use_zone_table && !d2.normal_to_ridge && d2.use_zone_table
    @test [l for (l, _) in d1.cases] == ["Windward roof (suction case)", "Windward roof (positive case)", "Leeward roof"]
    @test cases(d1)["Leeward roof"].Cp ≈ leeward_cp(o.theta_deg, 21.6 / 25.0)
    @test all(occursin("Case (", l) for (l, _) in d2.cases)
    # default roof_type is :gable
    @test compute_roof_pressures(; base..., roof_slope_rise_per_12=3.0, roof_type=:gable).results == o.results
    # Eq. 27.3-1 and ASD = 0.6 x LRFD
    c = cases(d1)["Leeward roof"]
    @test c.p_pos ≈ d1.qh * 0.85 * (0.85c.Cp - GCPI_ENCLOSED)
    a = compute_roof_pressures(; base..., roof_slope_rise_per_12=3.0, design_code="ASD")
    @test cases(a.results[1])["Leeward roof"].p_pos ≈ 0.6c.p_pos
end

@testset "hip roofs" begin
    # theta < 10 deg: zone table in both directions, identical to gable
    lo_h = compute_roof_pressures(; base..., roof_slope_rise_per_12=1.5, roof_type=:hip)
    lo_g = compute_roof_pressures(; base..., roof_slope_rise_per_12=1.5)
    @test lo_h.results == lo_g.results && lo_h.roof_type == :hip

    # theta >= 10 deg (3:12): windward/leeward faces in BOTH directions plus side-slope zone table
    o = compute_roof_pressures(; base..., roof_slope_rise_per_12=3.0, roof_type=:hip)
    g = compute_roof_pressures(; base..., roof_slope_rise_per_12=3.0)
    for (r, L) in zip(o.results, (25.0, 30.0))
        cs = cases(r)
        hl = 21.6 / L
        @test !r.use_zone_table && r.zones !== nothing
        @test cs["Windward roof (suction case)"].Cp ≈ windward_neg_cp(o.theta_deg, hl)
        @test cs["Windward roof (positive case)"].Cp ≈ windward_pos_cp(o.theta_deg, hl)
        @test cs["Leeward roof"].Cp ≈ leeward_cp(o.theta_deg, hl)
        @test count(l -> startswith(l, "Side slopes, "), keys(cs)) == 2 * length(r.zones)
        @test all(cs["Side slopes, $(z.label) ft, Case (b)"].Cp == ZONE_CP_B for z in r.zones)
        @test [z.cp_a for z in r.zones] == [z.cp_a for z in MWFRSRoofPressures.zone_table(21.6, L, hl, r.B, o.theta_deg, true)]
    end
    # Direction 1 faces match the gable normal-to-ridge cases; Direction 2 side slopes match the gable zone table
    for (l, c) in g.results[1].cases
        @test cases(o.results[1])[l] == c
    end
    for (l, c) in g.results[2].cases
        @test cases(o.results[2])["Side slopes, $l"] == c
    end
    # steep hip: positive windward case active
    s = compute_roof_pressures(; base..., roof_slope_rise_per_12=12.0, roof_type=:hip)
    @test cases(s.results[2])["Windward roof (positive case)"].Cp ≈ windward_pos_cp(45.0, 21.6 / 30.0)
end

end
