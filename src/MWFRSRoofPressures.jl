module MWFRSRoofPressures

export compute_roof_pressures

# =============================================================================
# ASCE 7-22 MWFRS Roof Wind Pressures -- Directional Procedure (Ch. 26-27)
# Gable roof, enclosed rigid building. Figure 27.3-1 roof Cp, Part 1 (Rigid
# buildings of all heights)
#
# Table values transcribed from ASCE 7-10/16/22 Figure 27.3-1 (unchanged
# across these editions)
#
# Usage from any script whose active environment has this package dev'd in:
#   using MWFRSRoofPressures
#   out = compute_roof_pressures(V=131.0, exposure="C", h=50.49,
#                                 long_dim=236.25, short_dim=121.25,
#                                 roof_slope_rise_per_12=0.25)
#   out.results[1].cases   # Direction 1 pressure cases (windward/leeward or zones)
# =============================================================================

# -----------------------------------------------------------------------
# CONSTANTS
# -----------------------------------------------------------------------
# Terrain exposure constants (alpha, zg [ft]), Table 26.11-1.
const TERRAIN_CONSTANTS = Dict(
    "B" => (alpha=7.5,  zg=3280.0),
    "C" => (alpha=9.8,  zg=2460.0),
    "D" => (alpha=11.5, zg=1935.0),
)

const GCPI_ENCLOSED = 0.18   # Table 26.13-1, enclosed buildings, +/-

# Roof Cp, normal to ridge, theta = 10-45 deg (Fig. 27.3-1), by h/L row.
const HL_ROWS    = (0.25, 0.5, 1.0)
const THETA_BREAKS = (10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 45.0)

const WW_NEG = Dict(   # windward roof, suction-case curve
    0.25 => (-0.70, -0.50, -0.30, -0.20, -0.20,  0.00,  0.00),
    0.5  => (-0.90, -0.70, -0.40, -0.30, -0.20, -0.20,  0.00),
    1.0  => (-1.30, -1.00, -0.70, -0.50, -0.30, -0.20,  0.00),
)
const WW_POS = Dict(   # windward roof, positive-case curve
    0.25 => (-0.18,  0.00,  0.20,  0.30,  0.30,  0.40,  0.40),
    0.5  => (-0.18, -0.18,  0.00,  0.20,  0.20,  0.30,  0.40),
    1.0  => (-0.18, -0.18, -0.18,  0.00,  0.20,  0.20,  0.30),
)

const LEEWARD_THETA_BREAKS = (10.0, 15.0, 20.0)
const LEEWARD = Dict(
    0.25 => (-0.3, -0.5, -0.6),
    0.5  => (-0.5, -0.5, -0.6),
    1.0  => (-0.7, -0.6, -0.6),
)

# Zone table: normal to ridge for theta < 10 deg, OR parallel to ridge for any
# theta (Fig. 27.3-1). Zones by horizontal distance from windward edge: 0-h/2,
# h/2-h, h-2h, >2h. h/L = 1.0 row's tabulated ">h/2" single value is applied to
# all three zones beyond 0-h/2 so both rows share 4 slots for interpolation.
const ZONE_CP_A = Dict(0.5 => (-0.90, -0.90, -0.50, -0.30), 1.0 => (-1.30, -0.70, -0.70, -0.70))
const ZONE_CP_B = -0.18   # constant across all zones and h/L for this table
const ZONE_LABELS = ("0 to h/2", "h/2 to h", "h to 2h", "> 2h")

# Area reduction for the leading (0 to h/2) zone's more-negative Cp value only
# (Fig. 27.3-1 "**" footnote); linear in area, sf.
const AREA_REDUCTION_PTS = [100.0 => 1.0, 200.0 => 0.9, 1000.0 => 0.8]

# -----------------------------------------------------------------------
# FUNCTIONS
# -----------------------------------------------------------------------

"""Kz (or Kh) exposure coefficient (Table 26.10-1 notes): Kz = 2.41*(z/zg)^(2/alpha), z clamped to [15, zg] ft."""
function velocity_pressure_coefficient(z::Real, exposure::AbstractString)
    c = TERRAIN_CONSTANTS[uppercase(exposure)]
    z_eff = clamp(z, 15.0, c.zg)
    return 2.41 * (z_eff / c.zg)^(2 / c.alpha)
end

"""Velocity pressure qz, psf (Eq. 26.10-1). Returns (qz, Kz)."""
function velocity_pressure(z::Real, V::Real, Kzt::Real, Ke::Real, exposure::AbstractString)
    Kz = velocity_pressure_coefficient(z, exposure)
    qz = 0.00256 * Kz * Kzt * Ke * V^2
    return qz, Kz
end

"""Linear interpolation of `vals` at ascending `breaks`; held flat outside [breaks[1], breaks[end]]."""
function lininterp(x::Real, breaks::NTuple{N,Float64}, vals::NTuple{N,Float64}) where N
    x <= breaks[1] && return vals[1]
    x >= breaks[end] && return vals[end]
    for i in 1:(N - 1)
        if breaks[i] <= x <= breaks[i + 1]
            frac = (x - breaks[i]) / (breaks[i + 1] - breaks[i])
            return vals[i] + frac * (vals[i + 1] - vals[i])
        end
    end
    return vals[end]
end

"""Bilinear (theta, h/L) lookup over a windward/leeward Cp grid keyed by HL_ROWS; h/L held flat outside [0.25, 1.0]."""
function interp_theta_hl(theta_deg::Real, h_over_L::Real, theta_breaks::NTuple, table::Dict)
    hl = clamp(h_over_L, HL_ROWS[1], HL_ROWS[end])
    if hl in HL_ROWS
        return lininterp(theta_deg, theta_breaks, table[hl])
    end
    lo_row = maximum(r for r in HL_ROWS if r <= hl)
    hi_row = minimum(r for r in HL_ROWS if r >= hl)
    v_lo = lininterp(theta_deg, theta_breaks, table[lo_row])
    v_hi = lininterp(theta_deg, theta_breaks, table[hi_row])
    frac = (hl - lo_row) / (hi_row - lo_row)
    return v_lo + frac * (v_hi - v_lo)
end

"""Windward roof suction-case Cp (theta 10-45 deg tabulated; 0.0 above 45 deg, Note 2: no value of like sign given)."""
function windward_neg_cp(theta_deg::Real, h_over_L::Real)
    theta_deg > 45.0 && return 0.0
    return interp_theta_hl(theta_deg, h_over_L, THETA_BREAKS, WW_NEG)
end

"""Windward roof positive-case Cp: tabulated 10-45 deg, blended to the theta>=60 deg 0.01*theta formula
over 45-60 deg, then 0.01*theta up to 80 deg, held at 0.8 above 80 deg (Fig. 27.3-1 "#" note, roof acts as a wall)."""
function windward_pos_cp(theta_deg::Real, h_over_L::Real)
    if theta_deg <= 45.0
        return interp_theta_hl(theta_deg, h_over_L, THETA_BREAKS, WW_POS)
    elseif theta_deg < 60.0
        v45 = interp_theta_hl(45.0, h_over_L, THETA_BREAKS, WW_POS)
        return lininterp(theta_deg, (45.0, 60.0), (v45, 0.6))
    elseif theta_deg <= 80.0
        return 0.01 * theta_deg
    else
        return 0.8
    end
end

"""Leeward roof Cp (theta 10-20 deg tabulated; held flat outside that range)."""
leeward_cp(theta_deg::Real, h_over_L::Real) = interp_theta_hl(theta_deg, h_over_L, LEEWARD_THETA_BREAKS, LEEWARD)

"""Zone-table Cp(a) for the 4 windward-edge zones (0-h/2, h/2-h, h-2h, >2h), interpolated in h/L over [0.5, 1.0] (held flat outside)."""
function zone_cp_a(h_over_L::Real)
    hl = clamp(h_over_L, 0.5, 1.0)
    row05, row10 = ZONE_CP_A[0.5], ZONE_CP_A[1.0]
    frac = (hl - 0.5) / (1.0 - 0.5)
    return ntuple(i -> row05[i] + frac * (row10[i] - row05[i]), 4)
end

"""Fig. 27.3-1 "**" area reduction factor for the leading zone's Cp(a), linear between the tabulated area breakpoints (sf)."""
function area_reduction_factor(area::Real)
    pts = AREA_REDUCTION_PTS
    area <= pts[1].first && return pts[1].second
    area >= pts[end].first && return pts[end].second
    for i in 1:(length(pts) - 1)
        lo, hi = pts[i], pts[i + 1]
        if lo.first <= area <= hi.first
            frac = (area - lo.first) / (hi.first - lo.first)
            return lo.second + frac * (hi.second - lo.second)
        end
    end
end

"""Zone rows (label, from, to, cp_a, cp_a_raw, reduction, cp_b, area) for the windward-edge zone table, clipped to L_dim."""
function zone_table(h::Real, L_dim::Real, h_over_L::Real, B::Real, theta_deg::Real, apply_area_reduction::Bool)
    cp_a = zone_cp_a(h_over_L)
    bounds = (0.0, h / 2, h, 2h, Inf)
    zones = NamedTuple[]
    for i in 1:4
        lo = bounds[i]
        lo >= L_dim && break
        hi = min(bounds[i + 1], L_dim)
        width = hi - lo
        area = width * B / cosd(theta_deg)   # sloped roof tributary area for this zone, ft^2
        reduction = 1.0
        cp_a_i = cp_a[i]
        if apply_area_reduction && i == 1
            reduction = area_reduction_factor(area)
            cp_a_i *= reduction
        end
        push!(zones, (label=ZONE_LABELS[i], from=lo, to=hi, cp_a=cp_a_i, cp_a_raw=cp_a[i],
                       reduction=reduction, cp_b=ZONE_CP_B, area=area))
    end
    return zones
end

"""Net design wind pressure (Eq. 27.3-1): p = q*Kd*G*Cp - qi*Kd*GCpi. Roof Cp always uses q = qi = qh."""
design_pressure(qh::Real, Kd::Real, G::Real, Cp::Real, GCpi_signed::Real) = qh * Kd * G * Cp - qh * Kd * GCpi_signed

"""Sec. 27.1.5 minimum roof pressure: magnitude not less than 8 psf (vs. 16 psf for walls)."""
roof_minimum(p::Real) = abs(p) < 8.0 ? copysign(8.0, p) : p

"""One Cp load case -> (Cp, external pressure, +/-GCpi net pressures, controlling sign)."""
function case_pressure(qh::Real, Kd::Real, G::Real, Cp::Real, apply_min::Bool)
    p_ext = qh * Kd * G * Cp
    p_pos = design_pressure(qh, Kd, G, Cp, +GCPI_ENCLOSED)
    p_neg = design_pressure(qh, Kd, G, Cp, -GCPI_ENCLOSED)
    if apply_min
        p_pos, p_neg = roof_minimum(p_pos), roof_minimum(p_neg)
    end
    controls = abs(p_pos) >= abs(p_neg) ? "+GCpi" : "-GCpi"
    return (Cp=Cp, p_ext=p_ext, p_pos=p_pos, p_neg=p_neg, controls=controls)
end

"""Roof pressure cases for one wind direction. B = building width normal to wind; L_dim = building depth parallel to
wind (used for h/L and, in the zone table, measured along this dimension from the windward edge)."""
function direction_pressures(label, B, L_dim, h, theta_deg, normal_to_ridge, V, exposure, Kd, Kzt, Ke, G,
                              apply_roof_minimum, apply_area_reduction)
    h_over_L = h / L_dim
    qh, Kh = velocity_pressure(h, V, Kzt, Ke, exposure)
    gcpi_psf = qh * Kd * GCPI_ENCLOSED

    use_zone_table = !normal_to_ridge || theta_deg < 10.0

    cases = Tuple{String,NamedTuple}[]
    zones = nothing
    if !use_zone_table
        cp_ww_neg = windward_neg_cp(theta_deg, h_over_L)
        cp_ww_pos = windward_pos_cp(theta_deg, h_over_L)
        cp_lw     = leeward_cp(theta_deg, h_over_L)
        push!(cases, ("Windward roof (suction case)",  case_pressure(qh, Kd, G, cp_ww_neg, apply_roof_minimum)))
        push!(cases, ("Windward roof (positive case)", case_pressure(qh, Kd, G, cp_ww_pos, apply_roof_minimum)))
        push!(cases, ("Leeward roof",                   case_pressure(qh, Kd, G, cp_lw,     apply_roof_minimum)))
    else
        zones = zone_table(h, L_dim, h_over_L, B, theta_deg, apply_area_reduction)
        for z in zones
            push!(cases, ("$(z.label) ft, Case (a)", case_pressure(qh, Kd, G, z.cp_a, apply_roof_minimum)))
            push!(cases, ("$(z.label) ft, Case (b)", case_pressure(qh, Kd, G, z.cp_b, apply_roof_minimum)))
        end
    end

    return (label=label, B=B, L=L_dim, h_over_L=h_over_L, qh=qh, Kh=Kh, gcpi_psf=gcpi_psf,
            normal_to_ridge=normal_to_ridge, use_zone_table=use_zone_table, theta_deg=theta_deg,
            cases=cases, zones=zones)
end

"""
    compute_roof_pressures(; V, exposure, h, long_dim, short_dim, roof_slope_rise_per_12,
                             risk_category="II", Kd=0.85, Kzt=1.0, Ke=1.0, G=0.85,
                             ridge_parallel_to="long", apply_roof_minimum=true, apply_area_reduction=true)

Compute ASCE 7-22 MWFRS gable roof pressures (Ch. 27, Directional Procedure,
Part 1: Rigid Buildings of All Heights, Fig. 27.3-1) for both orthogonal wind
directions of an enclosed rigid building.

Returns a NamedTuple:
  .results       -- 2-element vector (Direction 1, Direction 2), each with
                     .cases (label => case_pressure NamedTuple pairs) and,
                     when the zone table applies, .zones (raw zone geometry)
  .h, .theta_deg -- roof reference height and slope angle, deg
plus the resolved inputs for convenience in reporting.
"""
function compute_roof_pressures(; V, exposure, h, long_dim, short_dim, roof_slope_rise_per_12,
                                   risk_category="II", Kd=0.85, Kzt=1.0, Ke=1.0, G=0.85,
                                   ridge_parallel_to="long", apply_roof_minimum=true, apply_area_reduction=true)
    ridge_parallel_to in ("long", "short") || error("ridge_parallel_to must be \"long\" or \"short\"")
    theta_deg = rad2deg(atan(roof_slope_rise_per_12 / 12))

    dir1_normal_to_ridge = ridge_parallel_to == "long"
    dir2_normal_to_ridge = !dir1_normal_to_ridge

    results = [
        direction_pressures("Direction 1 (wind normal to the long wall; B=$long_dim ft, L=$short_dim ft)",
                             long_dim, short_dim, h, theta_deg, dir1_normal_to_ridge,
                             V, exposure, Kd, Kzt, Ke, G, apply_roof_minimum, apply_area_reduction),
        direction_pressures("Direction 2 (wind normal to the short wall; B=$short_dim ft, L=$long_dim ft)",
                             short_dim, long_dim, h, theta_deg, dir2_normal_to_ridge,
                             V, exposure, Kd, Kzt, Ke, G, apply_roof_minimum, apply_area_reduction),
    ]

    return (
        results=results, h=h, theta_deg=theta_deg,
        V=V, risk_category=risk_category, exposure=exposure, Kd=Kd, Kzt=Kzt, Ke=Ke, G=G,
        apply_roof_minimum=apply_roof_minimum, apply_area_reduction=apply_area_reduction,
        long_dim=long_dim, short_dim=short_dim, ridge_parallel_to=ridge_parallel_to,
    )
end

end # module MWFRSRoofPressures
