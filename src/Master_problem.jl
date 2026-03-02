# ─── Master Problem (MP) ─────────────────────────────────────────────────────

MP = Model(Gurobi.Optimizer)
set_optimizer_attribute(MP, "OutputFlag",    0)
set_optimizer_attribute(MP, "LogToConsole",  0)
set_optimizer_attribute(MP, "MIPGap",        1e-3)
set_optimizer_attribute(MP, "BarConvTol",    1e-3)
set_optimizer_attribute(MP, "OptimalityTol", 1e-3)

# ---- Investment variables ------------------------------------------------
@variable(MP, vNewPowGenCap[g in G] >= 0, Int)
@variable(MP, vRetPowGenCap[g in G] >= 0, Int)
@variable(MP, vNewPowStoCap[s in S] >= 0, Int)
@variable(MP, vRetPowStoCap[s in S] >= 0, Int)
@variable(MP, vNewPowTraCap[l in L] >= 0, Int)
@variable(MP, vNewH2GenCap[h in H]  >= 0, Int)
@variable(MP, vRetH2GenCap[h in H]  >= 0, Int)
@variable(MP, vNewH2StoCap[s in Q]  >= 0, Int)
@variable(MP, vRetH2StoCap[s in Q]  >= 0, Int)
@variable(MP, vNewH2Pipe[i in I]    >= 0, Int)
@variable(MP, vRetH2Pipe[i in I]    >= 0, Int)
@variable(MP, vNewH2PipeCompCap[i in I] >= 0, Int)
@variable(MP, vRetH2PipeCompCap[i in I] >= 0, Int)
@variable(MP, vNewH2StoCompCap[s in Q]  >= 0, Int)
@variable(MP, vRetH2StoCompCap[s in Q]  >= 0, Int)

# ---- Emission budget variable --------------------------------------------
@variable(MP, vMaxEmissionByWeek[w in W] >= 0)

# ---- Long-term-duration storage variables --------------------------------
@variable(MP, vH2SOCLast[s in Q, w in W]  >= 0)
@variable(MP, vH2SOCFirst[s in Q, w in W] >= 0)

# ---- Benders theta (recourse approximation) variable --------------------
@variable(MP, theta[w in W] >= 0)

# ---- Capacity expressions ------------------------------------------------
@expression(MP, eTotPowGenCap[g in G],
    pow_gen[g, :existing_cap] .+
    pow_gen[g, :rep_capacity] * (vNewPowGenCap[g] .- vRetPowGenCap[g]))

@expression(MP, eTotPowStoCap[s in S],
    pow_gen[s, :existing_cap] .+
    pow_gen[s, :rep_capacity] * (vNewPowStoCap[s] .- vRetPowStoCap[s]))

@expression(MP, eTotPowTraCap[l in L],
    pow_lines[l, :existing_transmission_cap_mw] .+ vNewPowTraCap[l])

@expression(MP, eTotH2GenCap[h in H],
    hsc_gen[h, :existing_cap_tonne_p_hr] .+
    hsc_gen[h, :rep_capacity] * (vNewH2GenCap[h] .- vRetH2GenCap[h]))

@expression(MP, eTotH2StoCap[s in Q],
    hsc_gen[s, :existing_cap_tonne] + (vNewH2StoCap[s] - vRetH2StoCap[s]))

@expression(MP, eTotH2StoCompCap[s in Q],
    hsc_gen[s, :existing_cap_comp_tonne_hr] + vNewH2StoCompCap[s] - vRetH2StoCompCap[s])

@expression(MP, eTotH2Pipe[i in I],
    hsc_pipelines[i, :existing_num_pipes] + vNewH2Pipe[i] - vRetH2Pipe[i])

@expression(MP, eTotH2PipeCompCap[i in I],
    hsc_pipelines[i, :existing_comp_cap_tonne_hr] + vNewH2PipeCompCap[i] - vRetH2PipeCompCap[i])

# ---- Land-use expressions ------------------------------------------------
@expression(MP, ePowGenLandUse[z in Z],
    sum((vNewPowGenCap[g] - vRetPowGenCap[g]) * pow_gen[g, :rep_capacity] *
        pow_gen[g, :total_land_use_km2_p_cap] * (pow_gen[g, :zone] == z ? 1 : 0)
        for g in G))

@expression(MP, ePowStoLandUse[z in Z],
    sum((vNewPowStoCap[s] - vRetPowStoCap[s]) * pow_gen[s, :rep_capacity] *
        pow_gen[s, :total_land_use_km2_p_cap] * (pow_gen[s, :zone] == z ? 1 : 0)
        for s in S))

@expression(MP, ePowLineLandUse[z in Z],
    0.5 * sum(pow_lines[l, :total_land_use_km2_p_length] * pow_lines[l, :distance_km] *
              vNewPowTraCap[l] * abs(Pow_Network[l, z])
              for l in L))

@expression(MP, eH2GenLandUse[z in Z],
    sum((vNewH2GenCap[h] - vRetH2GenCap[h]) * hsc_gen[h, :rep_capacity] *
        hsc_gen[h, :total_land_use_km2_p_cap] * (hsc_gen[h, :zone] == z ? 1 : 0)
        for h in H))

@expression(MP, eH2StoLandUse[z in Z],
    sum((vNewH2StoCap[s] - vRetH2StoCap[s]) * hsc_gen[s, :total_land_use_km2_p_cap] *
        (hsc_gen[s, :zone] == z ? 1 : 0)
        for s in Q))

@expression(MP, eH2PipeLandUse[z in Z],
    0.5 * sum(hsc_pipelines[i, :total_land_use_km2_p_length] * hsc_pipelines[i, :distance] *
              (vNewH2Pipe[i] - vRetH2Pipe[i]) * abs(H2_Network[i, z])
              for i in I))

@expression(MP, eTotalLandUse[z in Z],
    ePowGenLandUse[z] + ePowStoLandUse[z] + ePowLineLandUse[z] +
    eH2GenLandUse[z]  + eH2StoLandUse[z]  + eH2PipeLandUse[z])

# ---- Investment cost expressions -----------------------------------------
@expression(MP, eCostPowGenInv,
    sum(pow_gen[g, :inv_cost_per_mwyr] .* vNewPowGenCap[g] .* pow_gen[g, :rep_capacity] .+
        pow_gen[g, :fom_cost_per_mwyr] .* eTotPowGenCap[g]
        for g in G))

@expression(MP, eCostPowStoInv,
    sum(pow_gen[s, :inv_cost_per_mwhyr] .* vNewPowStoCap[s] * pow_gen[s, :rep_capacity] .+
        pow_gen[s, :fom_cost_per_mwhyr] .* eTotPowStoCap[s]
        for s in S))

@expression(MP, eCostPowTraInv,
    sum(pow_lines[l, :line_reinforcement_cost_per_mwyr] .* vNewPowTraCap[l] for l in L))

@expression(MP, eCostH2GenInv,
    sum(hsc_gen[h, :inv_cost_tonne_hr_p_yr] * vNewH2GenCap[h] * hsc_gen[h, :rep_capacity] +
        hsc_gen[h, :fom_cost_p_tonne_p_hr_yr] * eTotH2GenCap[h]
        for h in H))

@expression(MP, eCostH2StoInv,
    sum(hsc_gen[s, :inv_cost_tonne_p_yr]       * vNewH2StoCap[s] +
        hsc_gen[s, :inv_cost_comp_tonne_hr_p_yr] * vNewH2StoCompCap[s] +
        hsc_gen[s, :fom_cost_p_tonne_p_yr]      * eTotH2StoCap[s] +
        hsc_gen[s, :fom_cost_comp_tonne_hr_p_yr] * eTotH2StoCompCap[s]
        for s in Q))

@expression(MP, eCostH2TraInv,
    sum(hsc_pipelines[i, :investment_cost_per_length] * hsc_pipelines[i, :distance] * vNewH2Pipe[i] +
        hsc_pipelines[i, :compressor_inv_per_length]  * hsc_pipelines[i, :distance] * vNewH2PipeCompCap[i]
        for i in I))

# ---- Objective -----------------------------------------------------------
MP_obj = eCostPowGenInv .+ eCostPowStoInv .+ eCostPowTraInv .+
         eCostH2GenInv  .+ eCostH2StoInv  .+ eCostH2TraInv

@objective(MP, Min, MP_obj .+ sum(theta[w] for w in W))

# ---- Capacity constraints ------------------------------------------------
@constraint(MP, cMaxPowGenRetCapTher[g in G],
    vRetPowGenCap[g] * pow_gen[g, :rep_capacity] <= pow_gen[g, :existing_cap])
@constraint(MP, cMaxPowGenCap[g in G],   eTotPowGenCap[g] <= pow_gen[g, :max_cap_mw])
@constraint(MP, cPowGenTotCapMin[g in G], pow_gen[g, :min_cap] <= eTotPowGenCap[g])

@constraint(MP, cMaxRetPowSto[s in S],
    vRetPowStoCap[s] * pow_gen[s, :rep_capacity] <= pow_gen[s, :existing_cap_mwh])
@constraint(MP, cMaxPowStoCap[s in S],    eTotPowStoCap[s] <= pow_gen[s, :max_cap_mwh])
@constraint(MP, cPowStoTotCapMin[s in S], pow_gen[s, :min_cap_mwh] .- eTotPowStoCap[s] <= 0)

@constraint(MP, cMaxPowTraCap[l in L],
    vNewPowTraCap[l] <= pow_lines[l, :line_max_reinforcement_mw])

@constraint(MP, cMaxRetH2Cap[h in H],
    vRetH2GenCap[h] * hsc_gen[h, :rep_capacity] <= hsc_gen[h, :existing_cap_tonne_p_hr])
@constraint(MP, cMaxH2GenCap[h in H],  eTotH2GenCap[h] <= hsc_gen[h, :max_cap_tonne_p_hr])
@constraint(MP, cMinH2GenCap[h in H],  hsc_gen[h, :min_cap_tonne_p_hr] <= eTotH2GenCap[h])

@constraint(MP, cMaxRetH2StoCap[s in Q],
    vRetH2StoCap[s] <= hsc_gen[s, :existing_cap_tonne])
@constraint(MP, cMaxH2StoCap[s in Q],
    eTotH2StoCap[s] .- hsc_gen[s, :max_cap_stor_tonne] <= 0)
@constraint(MP, cMinH2StoCap[s in Q],
    hsc_gen[s, :min_cap_stor_tonne] <= eTotH2StoCap[s])

@constraint(MP, cMaxRetH2StorCompCap[s in Q],
    vRetH2StoCompCap[s] <= hsc_gen[s, :existing_cap_comp_tonne_hr])
@constraint(MP, cMaxH2StorCompcCap[s in Q],
    eTotH2StoCompCap[s] <= eTotH2StoCap[s])
@constraint(MP, cMinH2StorCompcCap[s in Q],
    0.01 * eTotH2StoCap[s] <= eTotH2StoCompCap[s])

@constraint(MP, cMaxH2PipeNum[i in I],
    eTotH2Pipe[i] <= hsc_pipelines[i, :max_num_pipes])
@constraint(MP, cMaxRetH2PipeNum[i in I],
    vRetH2Pipe[i] <= hsc_pipelines[i, :existing_num_pipes])
@constraint(MP, cMaxRetH2PipeCompCap[i in I],
    vRetH2PipeCompCap[i] <= hsc_pipelines[i, :existing_comp_cap_tonne_hr])
@constraint(MP, cMaxH2PipeComp[i in I],
    eTotH2PipeCompCap[i] <= eTotH2Pipe[i] * hsc_pipelines[i, :max_pipe_cap_tonne])

# ---- Land-use constraint -------------------------------------------------
@constraint(MP, cLandUse[z in Z],
    eTotalLandUse[z] .- zones[z, :maximum_available_land] <= 0)

# ---- Emission policy constraint ------------------------------------------
@constraint(MP, cZonalEmissionCap,
    sum(vMaxEmissionByWeek[w] for w in W) <=
    0.0 * 0.4 * sum((h2_demand[w][t, z] * 33.3) + pow_demand[w][t, z]
                    for t in T, w in W, z in Z))

# ---- Long-term-duration storage linking constraints ----------------------
@constraint(MP, cH2StoSOCLink[s in Q, w in 2:length(W)],
    vH2SOCFirst[s, w] == vH2SOCLast[s, w-1])
@constraint(MP, cH2StoSOCLinkCycle[s in Q],
    vH2SOCFirst[s, 1] == vH2SOCLast[s, 52])
@constraint(MP, cH2SOCFirstMax[s in Q, w in W],
    vH2SOCFirst[s, w] <= hsc_gen[s, :h2stor_max_level] * eTotH2StoCap[s])
@constraint(MP, cH2SOCLastMax[s in Q, w in W],
    vH2SOCLast[s, w]  <= hsc_gen[s, :h2stor_max_level] * eTotH2StoCap[s])
@constraint(MP, cH2SOCFirstMin[s in Q, w in W],
    hsc_gen[s, :h2stor_min_level] * eTotH2StoCap[s] <= vH2SOCFirst[s, w])
@constraint(MP, cH2SOCLastMin[s in Q, w in W],
    hsc_gen[s, :h2stor_min_level] * eTotH2StoCap[s] <= vH2SOCLast[s, w])
