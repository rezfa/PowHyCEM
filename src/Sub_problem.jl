# ─── Sub-Problems (SP) ───────────────────────────────────────────────────────

SP_models = Vector{Model}(undef, length(W))

for w in W
    SP_models[w] = Model(Gurobi.Optimizer)
    set_optimizer_attribute(SP_models[w], "Threads",      1)
    set_optimizer_attribute(SP_models[w], "Method",       2)
    set_optimizer_attribute(SP_models[w], "Crossover",    1)
    set_optimizer_attribute(SP_models[w], "OutputFlag",   0)
    set_optimizer_attribute(SP_models[w], "LogToConsole", 0)
    set_optimizer_attribute(SP_models[w], "MIPGap",       1e-3)
    set_optimizer_attribute(SP_models[w], "BarConvTol",   1e-3)
    set_optimizer_attribute(SP_models[w], "OptimalityTol",1e-3)

    # ---- SP Variables --------------------------------------------------------
    # Power generation
    @variable(SP_models[w], vPowGen[g in G, t in T] >= 0)
    @variable(SP_models[w], vPowGenFirst[g in G]    >= 0)
    @variable(SP_models[w], vPowGenOnline[g in G_ther, t in T]      >= 0)
    @variable(SP_models[w], vPowGenOnlineFirst[g in G_ther]          >= 0)
    @variable(SP_models[w], vPowGenStart[g in G_ther, t in T]        >= 0)
    @variable(SP_models[w], vPowGenShut[g in G_ther, t in T]         >= 0)
    @variable(SP_models[w], vPowResUp[g in G_ther, t in T]           >= 0)
    @variable(SP_models[w], vPowResDn[g in G_ther, t in T]           >= 0)
    # Power storage
    @variable(SP_models[w], vPowStoCha[s in S, t in T]  >= 0)
    @variable(SP_models[w], vPowStoDis[s in S, t in T]  >= 0)
    @variable(SP_models[w], vPowSOC[s in S, t in T]     >= 0)
    @variable(SP_models[w], vPowSOCFirst[s in S]         >= 0)
    # Power transmission
    @variable(SP_models[w], vPowFlow[l in L, t in T])
    # H2 generation
    @variable(SP_models[w], vH2Gen[h in H, t in T]           >= 0)
    @variable(SP_models[w], vH2GenFirst[h in H]               >= 0)
    @variable(SP_models[w], vH2GenStart[h in H_ther, t in T]  >= 0)
    @variable(SP_models[w], vH2GenShut[h in H_ther, t in T]   >= 0)
    @variable(SP_models[w], vH2GenOnline[h in H_ther, t in T]       >= 0)
    @variable(SP_models[w], vH2GenOnlineFirst[h in H_ther]            >= 0)
    # H2 storage
    @variable(SP_models[w], vH2StoCha[s in Q, t in T]   >= 0)
    @variable(SP_models[w], vH2StoDis[s in Q, t in T]   >= 0)
    @variable(SP_models[w], vH2StoSOC[s in Q, t in T]   >= 0)
    @variable(SP_models[w], vH2StoSOCFirst[s in Q]       >= 0)
    # H2 transmission
    @variable(SP_models[w], vH2FlowPos[i in I, t in T] >= 0)
    @variable(SP_models[w], vH2FlowNeg[i in I, t in T] >= 0)
    # Policy / slack variables
    @variable(SP_models[w], vH2NSD[z in Z, t in T]   >= 0)
    @variable(SP_models[w], vPowNSD[z in Z, t in T]  >= 0)
    @variable(SP_models[w], vPowCrt[z in Z, t in T]  >= 0)
    @variable(SP_models[w], vH2Crt[z in Z, t in T]   >= 0)
    @variable(SP_models[w], vExtraEmission             >= 0)
    # H2 storage SOC coupling variables
    @variable(SP_models[w], eH2SOCLast[s in Q]      >= 0)
    @variable(SP_models[w], eH2SOCFirst[s in Q]     >= 0)
    @variable(SP_models[w], slack_end[s in Q]        >= 0)
    @variable(SP_models[w], slack_start[s in Q]      >= 0)
    @variable(SP_models[w], slack_start_neg[s in Q]  >= 0)
    # Reserve slack variables
    @variable(SP_models[w], vSlackResUp[z in Z, t in T] >= 0)
    @variable(SP_models[w], vSlackResDn[z in Z, t in T] >= 0)
    # Capacity availability parameters (fixed by MP solution via coupling constraints)
    @variable(SP_models[w], eAvailPowGenCap[g in G]         >= 0)
    @variable(SP_models[w], eAvailPowGenUnit[g in G_ther]   >= 0)
    @variable(SP_models[w], eAvailPowStoCap[s in S]         >= 0)
    @variable(SP_models[w], eAvailPowTraCap[l in L]         >= 0)
    @variable(SP_models[w], eAvailH2GenCap[h in H]          >= 0)
    @variable(SP_models[w], eAvailH2GenUnit[h in H_ther]    >= 0)
    @variable(SP_models[w], eAvailH2StoCap[s in Q]          >= 0)
    @variable(SP_models[w], eAvailH2Pipe[i in I]            >= 0)
    @variable(SP_models[w], eAvailH2PipeCompCap[i in I]     >= 0)
    @variable(SP_models[w], eAvailH2StoCompCap[s in Q]      >= 0)
    @variable(SP_models[w], eMaxEmissionByWeek               >= 0)

    # ---- SP Expressions ------------------------------------------------------
    # Power generation
    @expression(SP_models[w], ePowGenByZone[z in Z, t in T],
        sum(vPowGen[g, t] * (pow_gen[g, :zone] == z ? 1 : 0) for g in G))
    @expression(SP_models[w], eTotPowGenCapByZone[z in Z],
        sum(eAvailPowGenCap[g] * (pow_gen[g, :zone] == z ? 1 : 0) for g in G_ther))
    # Power storage
    @expression(SP_models[w], ePowStoChaByZone[z in Z, t in T],
        sum(vPowStoCha[s, t] * (pow_gen[s, :zone] == z ? 1 : 0) for s in S))
    @expression(SP_models[w], ePowStoDisByZone[z in Z, t in T],
        sum(vPowStoDis[s, t] * (pow_gen[s, :zone] == z ? 1 : 0) for s in S))
    # Power transmission
    @expression(SP_models[w], eNet_Pow_Flow[z in Z, t in T],
        sum(Pow_Network[l, z] * vPowFlow[l, t] for l in L))
    @expression(SP_models[w], ePow_Loss_By_Zone[z in Z, t in T],
        sum(abs(Pow_Network[l, z]) * 0.5 * vPowFlow[l, t] * pow_lines[l, :line_loss_percentage]
            for l in L))
    # H2 generation
    @expression(SP_models[w], eH2FlowNet[i in I, t in T],
        vH2FlowPos[i, t] - vH2FlowNeg[i, t])
    @expression(SP_models[w], eH2GenEvap[h in H, t in T],
        hsc_gen[h, :boil_off] * vH2Gen[h, t])
    @expression(SP_models[w], eH2GenByZone[z in Z, t in T],
        sum((vH2Gen[h, t] - eH2GenEvap[h, t]) * (hsc_gen[h, :zone] == z ? 1 : 0) for h in H))
    # H2 storage
    @expression(SP_models[w], eH2StoChaByZone[z in Z, t in T],
        sum(vH2StoCha[s, t] * (hsc_gen[s, :zone] == z ? 1 : 0) for s in Q))
    @expression(SP_models[w], eH2StoDisByZone[z in Z, t in T],
        sum(vH2StoDis[s, t] * (hsc_gen[s, :zone] == z ? 1 : 0) for s in Q))
    # H2 transmission
    @expression(SP_models[w], eH2FlowByZone[i in I, z in Z, t in T],
        H2_Network[i, z] * eH2FlowNet[i, t])
    @expression(SP_models[w], eNet_H2_Flow[z in Z, t in T],
        sum(H2_Network[i, z] * eH2FlowNet[i, t] for i in I))
    @expression(SP_models[w], eH2_Loss_By_Zone[z in Z, t in T],
        sum(abs(H2_Network[i, z]) * 0.5 * eH2FlowNet[i, t] * hsc_pipelines[i, :pipe_loss_coeff]
            for i in I))
    # Sector-coupling demands
    @expression(SP_models[w], ePowDemandHSC[z in Z, t in T],
        sum(hsc_gen[h, :pow_demand_mwh_p_tonne] * vH2Gen[h, t] *
            (hsc_gen[h, :zone] == z ? 1 : 0) for h in H) +
        sum(hsc_gen[s, :h2charge_mwh_p_tonne] * vH2StoCha[s, t] *
            (hsc_gen[s, :zone] == z ? 1 : 0) for s in Q))
    @expression(SP_models[w], eH2DemandPow[z in Z, t in T],
        sum(pow_gen[g, :h2_demand_tonne_p_mwh] * vPowGen[g, t] *
            (pow_gen[g, :zone] == z ? 1 : 0) for g in G))
    @expression(SP_models[w], pow_D[t in T, z in Z],
        pow_demand[w][t, z] .+ ePowDemandHSC[z, t])
    @expression(SP_models[w], H2_D[t in T, z in Z],
        h2_demand[w][t, z] .+ eH2DemandPow[z, t])
    # Reserve requirements
    @expression(SP_models[w], ePowResReqUp[z in Z, t in T], 0.10 * pow_demand[w][t, z])
    @expression(SP_models[w], ePowResReqDn[z in Z, t in T], 0.05 * pow_demand[w][t, z])
    # Emissions
    @expression(SP_models[w], eEmissionByWeek,
        sum(vPowGen[g, t]      * pow_gen[g, :heat_rate_mmbtu_per_mwh]   * CO2_content[pow_gen[g, :fuel]]
            for g in G, t in T) +
        sum(vH2Gen[h, t]       * hsc_gen[h, :heat_rate_mmbtu_p_tonne]   * CO2_content[hsc_gen[h, :fuel]]
            for h in H, t in T) +
        sum(vPowGenStart[g, t] * pow_gen[g, :heat_rate_mmbtu_per_mwh]   * CO2_content[pow_gen[g, :fuel]]
            for g in G_ther, t in T) +
        sum(vH2GenStart[h, t]  * hsc_gen[h, :heat_rate_mmbtu_p_tonne]   * CO2_content[hsc_gen[h, :fuel]]
            for h in H_ther, t in T))
    @expression(SP_models[w], eEmissionByWeekZone[z in Z],
        sum(vPowGen[g, t]      * pow_gen[g, :heat_rate_mmbtu_per_mwh]   * CO2_content[pow_gen[g, :fuel]]   * (pow_gen[g, :zone] == z ? 1 : 0) for g in G, t in T) +
        sum(vH2Gen[h, t]       * hsc_gen[h, :heat_rate_mmbtu_p_tonne]   * CO2_content[hsc_gen[h, :fuel]]   * (hsc_gen[h, :zone] == z ? 1 : 0) for h in H, t in T) +
        sum(vPowGenStart[g, t] * pow_gen[g, :heat_rate_mmbtu_per_mwh]   * CO2_content[pow_gen[g, :fuel]]   * (pow_gen[g, :zone] == z ? 1 : 0) for g in G_ther, t in T) +
        sum(vH2GenStart[h, t]  * hsc_gen[h, :heat_rate_mmbtu_p_tonne]   * CO2_content[hsc_gen[h, :fuel]]   * (hsc_gen[h, :zone] == z ? 1 : 0) for h in H_ther, t in T))
    # Operating costs
    @expression(SP_models[w], eCostPowGenVar,
        sum((pow_gen[g, :vom_cost_mwh] +
             pow_gen[g, :heat_rate_mmbtu_per_mwh] .* fuel_costs[pow_gen[g, :fuel]][t]) .*
            vPowGen[g, t] for g in G, t in T))
    @expression(SP_models[w], eCostPowGenStart,
        sum(pow_gen[g, :start_cost_per_mw] .* pow_gen[g, :rep_capacity] .* vPowGenStart[g, t]
            for g in G_ther, t in T))
    @expression(SP_models[w], eCostPowStoVar,
        sum(vPowStoCha[s, t] .* pow_gen[s, :vom_cost_mwh_charge] for s in S, t in T))
    @expression(SP_models[w], eCostH2GenVar,
        sum((hsc_gen[h, :vom_cost_p_tonne] +
             hsc_gen[h, :heat_rate_mmbtu_p_tonne] .* fuel_costs[lowercase(hsc_gen[h, :fuel])][t]) .*
            vH2Gen[h, t] for h in H, t in T))
    @expression(SP_models[w], eCostH2StoVar,
        sum(hsc_gen[s, :var_om_cost_charge_p_tonne] * vH2StoCha[s, t] for s in Q, t in T))
    @expression(SP_models[w], eCostH2GenStart,
        sum(hsc_gen[h, :startup_cost_p_tonne_hr] .* hsc_gen[h, :rep_capacity] .* vH2GenStart[h, t]
            for h in H_ther, t in T))
    @expression(SP_models[w], eCostH2TraVar,
        sum(hsc_pipelines[i, :vom_per_tonne] * (vH2FlowPos[i, t] + vH2FlowNeg[i, t])
            for i in I, t in T))
    @expression(SP_models[w], eCostH2NSD,
        sum(vH2NSD[z, t] .* zones[z, :voll_hsc] for z in Z, t in T))
    @expression(SP_models[w], eCostPowNSD,
        sum(vPowNSD[z, t] .* zones[z, :voll_pow] for z in Z, t in T))
    @expression(SP_models[w], eCostPowCrt,
        sum(vPowCrt[z, t] .* zones[z, :pow_curtail_cost] for z in Z, t in T))
    @expression(SP_models[w], eCostH2Crt,
        sum(vH2Crt[z, t] .* zones[z, :h2_curtail_cost] for z in Z, t in T))
    @expression(SP_models[w], eEmissionCost,
        vExtraEmission * zones[1, :emission_cost])
    @expression(SP_models[w], SlackCost,
        sum((slack_end[s] + slack_start[s] + slack_start_neg[s]) * 1_000_000 for s in Q))
    @expression(SP_models[w], eCostResSlack,
        sum(vSlackResUp[z, t] * zones[z, :voll_pow] + vSlackResDn[z, t] * zones[z, :voll_pow]
            for z in Z, t in T))

    # ---- SP Objective --------------------------------------------------------
    @objective(SP_models[w], Min,
        eCostPowGenVar  .+ eCostPowStoVar  .+ eCostPowGenStart .+
        eCostH2GenVar   .+ eCostH2StoVar   .+ eCostH2GenStart  .+
        eCostPowNSD     .+ eCostH2NSD      .+ eCostPowCrt      .+
        eCostH2Crt      .+ eEmissionCost   .+ SlackCost        .+ eCostResSlack)

    # ---- SP Constraints ------------------------------------------------------
    # Energy balance
    @constraint(SP_models[w], cPowerBalance[z in Z, t in T],
        ePowGenByZone[z, t] .+ eNet_Pow_Flow[z, t] .- ePow_Loss_By_Zone[z, t] .+
        ePowStoDisByZone[z, t] .- ePowStoChaByZone[z, t] .+ vPowNSD[z, t] .- vPowCrt[z, t] ==
        pow_D[t, z])
    @constraint(SP_models[w], cH2Balance[z in Z, t in T],
        eH2GenByZone[z, t] .+ eNet_H2_Flow[z, t] .- eH2_Loss_By_Zone[z, t] .+
        eH2StoDisByZone[z, t] .- eH2StoChaByZone[z, t] .+ vH2NSD[z, t] .- vH2Crt[z, t] ==
        H2_D[t, z])

    # Power generation limits
    @constraint(SP_models[w], cMaxPowGen[g in G_ren, t in T],
        vPowGen[g, t] - eAvailPowGenCap[g] * pow_gen_var[H_w[w][t], pow_gen[g, :resource]] <= 0)
    @constraint(SP_models[w], cPowOnlineUnits[g in G_ther, t in T],
        vPowGenOnline[g, t] .- eAvailPowGenUnit[g] <= 0)
    @constraint(SP_models[w], cPowStartLimits[g in G_ther, t in T],
        vPowGenStart[g, t] .- eAvailPowGenUnit[g] .+ vPowGenOnline[g, t] <= 0)
    @constraint(SP_models[w], cMaxPowGenTher[g in G_ther, t in T],
        vPowGen[g, t] .- (pow_gen[g, :rep_capacity] .* pow_gen[g, :max_op_level] .* vPowGenOnline[g, t]) <= 0)
    @constraint(SP_models[w], cMinPowGenTher[g in G_ther, t in T],
        (pow_gen[g, :rep_capacity] * pow_gen[g, :min_op_level] * vPowGenOnline[g, t]) .- vPowGen[g, t] <= 0)
    @constraint(SP_models[w], cPowShutLimits[g in G_ther, t in T],
        vPowGenShut[g, t] .- vPowGenOnline[g, t] <= 0)

    # Spinning reserve
    @constraint(SP_models[w], cPowResUpMax[g in G_ther, t in T],
        vPowResUp[g, t] .+ vPowGen[g, t] .-
        pow_gen[g, :rep_capacity] * pow_gen[g, :max_op_level] * vPowGenOnline[g, t] <= 0)
    @constraint(SP_models[w], cPowResDnMax[g in G_ther, t in T],
        vPowResDn[g, t] .+ pow_gen[g, :rep_capacity] * pow_gen[g, :min_op_level] *
        vPowGenOnline[g, t] .- vPowGen[g, t] <= 0)
    @constraint(SP_models[w], cPowResUP[g in G_ther, t in T],
        vPowResUp[g, t] .- eAvailPowGenCap[g] * pow_gen[g, :ramp_up] <= 0)
    @constraint(SP_models[w], cPowResDn[g in G_ther, t in T],
        vPowResDn[g, t] .- eAvailPowGenCap[g] * pow_gen[g, :ramp_dn] <= 0)
    @constraint(SP_models[w], cPowResReqUp[z in Z, t in T],
        ePowResReqUp[z, t] .-
        sum(vPowResUp[g, t] * (pow_gen[g, :zone] == z ? 1 : 0) for g in G_ther) <= vSlackResUp[z, t])
    @constraint(SP_models[w], cPowResReqDn[z in Z, t in T],
        ePowResReqDn[z, t] .-
        sum(vPowResDn[g, t] * (pow_gen[g, :zone] == z ? 1 : 0) for g in G_ther) <= vSlackResDn[z, t])

    # Thermal unit commitment (cyclic)
    @constraint(SP_models[w], cPowUnitOnlineCon[g in G_ther, t in 2:length(T)],
        vPowGenOnline[g, t] - vPowGenOnline[g, t-1] == vPowGenStart[g, t] - vPowGenShut[g, t])
    @constraint(SP_models[w], cPowUnitOnlineFirst[g in G_ther],
        vPowGenOnline[g, 1] - vPowGenOnlineFirst[g] == vPowGenStart[g, 1] - vPowGenShut[g, 1])
    @constraint(SP_models[w], cPowUnitOnlineCycle[g in G_ther],
        vPowGenOnlineFirst[g] == vPowGenOnline[g, 168])

    # Ramp constraints — renewable units
    @constraint(SP_models[w], cPowGenRampUp[g in G_ren, t in 2:length(T)],
        vPowGen[g, t] - vPowGen[g, t-1] .- pow_gen[g, :ramp_up] * eAvailPowGenCap[g] <= 0)
    @constraint(SP_models[w], cPowGenRampUpFirst[g in G_ren],
        vPowGen[g, 1] - vPowGenFirst[g] .- pow_gen[g, :ramp_up] * eAvailPowGenCap[g] <= 0)
    @constraint(SP_models[w], cPowGenRampDn[g in G_ren, t in 2:length(T)],
        vPowGen[g, t-1] - vPowGen[g, t] .- pow_gen[g, :ramp_dn] * eAvailPowGenCap[g] <= 0)
    @constraint(SP_models[w], cPowGenRampDnFirst[g in G_ren],
        vPowGenFirst[g] - vPowGen[g, 1] .- pow_gen[g, :ramp_dn] * eAvailPowGenCap[g] <= 0)
    @constraint(SP_models[w], cPowGenCycle[g in G_ren],
        vPowGenFirst[g] == vPowGen[g, 168])

    # Ramp constraints — thermal units
    @constraint(SP_models[w], cTherPowGenRampDn[g in G_ther, t in 2:length(T)],
        vPowGen[g, t-1] .- vPowGen[g, t] .-
        pow_gen[g, :rep_capacity] * pow_gen[g, :ramp_dn] * (vPowGenOnline[g, t] .- vPowGenStart[g, t]) .+
        pow_gen[g, :rep_capacity] * vPowGenStart[g, t]   * pow_gen[g, :min_op_level] .-
        vPowGenShut[g, t] * pow_gen[g, :rep_capacity] *
        min(pow_gen_var[H_w[w][t], pow_gen[g, :resource]],
            max(pow_gen[g, :min_op_level], pow_gen[g, :ramp_dn])) <= 0)
    @constraint(SP_models[w], cTherPowGenRampDnFirst[g in G_ther],
        vPowGenFirst[g] .- vPowGen[g, 1] .-
        pow_gen[g, :rep_capacity] * pow_gen[g, :ramp_dn] * (vPowGenOnline[g, 1] .- vPowGenStart[g, 1]) .+
        pow_gen[g, :rep_capacity] * vPowGenStart[g, 1]   * pow_gen[g, :min_op_level] .-
        vPowGenShut[g, 1] * pow_gen[g, :rep_capacity] *
        min(pow_gen_var[H_w[w][1], pow_gen[g, :resource]],
            max(pow_gen[g, :min_op_level], pow_gen[g, :ramp_dn])) <= 0)
    @constraint(SP_models[w], cTherPowGenRampUp[g in G_ther, t in 2:length(T)],
        vPowGen[g, t] .- vPowGen[g, t-1] .-
        pow_gen[g, :rep_capacity] * pow_gen[g, :ramp_up] * (vPowGenOnline[g, t] .- vPowGenStart[g, t]) .+
        pow_gen[g, :rep_capacity] * vPowGenShut[g, t]    * pow_gen[g, :min_op_level] .-
        vPowGenStart[g, t] * pow_gen[g, :rep_capacity] *
        min(pow_gen_var[H_w[w][t], pow_gen[g, :resource]],
            max(pow_gen[g, :min_op_level], pow_gen[g, :ramp_up])) <= 0)
    @constraint(SP_models[w], cTherPowGenRampUpFirst[g in G_ther],
        vPowGen[g, 1] .- vPowGenFirst[g] .-
        pow_gen[g, :rep_capacity] * pow_gen[g, :ramp_up] * (vPowGenOnline[g, 1] .- vPowGenStart[g, 1]) .-
        pow_gen[g, :rep_capacity] * vPowGenShut[g, 1]    * pow_gen[g, :min_op_level] .+
        vPowGenStart[g, 1] * pow_gen[g, :rep_capacity] *
        min(pow_gen_var[H_w[w][1], pow_gen[g, :resource]],
            max(pow_gen[g, :min_op_level], pow_gen[g, :ramp_up])) <= 0)

    # Minimum up/down time — thermal power
    @constraint(SP_models[w], cMinUpTimePowGen[g in G_ther, t in T],
        sum(vPowGenStart[g, tt] for tt in intersect(T, (t - pow_gen[g, :up_time]):t)) .-
        vPowGenOnline[g, t] <= 0)
    @constraint(SP_models[w], cMinDnTimePowGen[g in G_ther, t in T],
        sum(vPowGenShut[g, tt] for tt in intersect(T, (t - pow_gen[g, :down_time]):t)) .-
        eAvailPowGenUnit[g] .+ vPowGenOnline[g, t] <= 0)

    # Power storage
    @constraint(SP_models[w], cPowStoBalance[s in S, t in 2:length(T)],
        vPowSOC[s, t] == (1 - pow_gen[s, :etta_self_dis]) * vPowSOC[s, t-1] +
        pow_gen[s, :etta_cha] * vPowStoCha[s, t] - (1/pow_gen[s, :etta_dis]) * vPowStoDis[s, t])
    @constraint(SP_models[w], cPowStoBalanceFirst[s in S],
        vPowSOC[s, 1] == (1 - pow_gen[s, :etta_self_dis]) * vPowSOCFirst[s] +
        pow_gen[s, :etta_cha] * vPowStoCha[s, 1] .- (1/pow_gen[s, :etta_dis]) * vPowStoDis[s, 1])
    @constraint(SP_models[w], cPowStoMaxDis[s in S, t in 2:length(T)],
        vPowStoDis[s, t] .- pow_gen[s, :etta_dis] * vPowSOC[s, t-1] <= 0)
    @constraint(SP_models[w], cPowStoMaxDisFirst[s in S],
        vPowStoDis[s, 1] .- pow_gen[s, :etta_dis] * vPowSOCFirst[s] <= 0)
    @constraint(SP_models[w], cPowSOCCycle[s in S],
        vPowSOCFirst[s] == vPowSOC[s, 168])
    @constraint(SP_models[w], cPowStoMaxCha[s in S, t in T],
        vPowStoCha[s, t] .- eAvailPowStoCap[s] <= 0)
    @constraint(SP_models[w], cPowStoSOCMax[s in S, t in T],
        vPowSOC[s, t] .- eAvailPowStoCap[s] * pow_gen[s, :max_op_level] <= 0)
    @constraint(SP_models[w], cPowStoSOCMin[s in S, t in T],
        eAvailPowStoCap[s] * pow_gen[s, :min_op_level] .- vPowSOC[s, t] <= 0)

    # Power transmission
    @constraints(SP_models[w], begin
        cMaxPowFlowOut[l in L, t in T],  vPowFlow[l, t] .- eAvailPowTraCap[l] <= 0
        cMaxPowFlowIn[l in L, t in T],  -eAvailPowTraCap[l] .- vPowFlow[l, t] <= 0
    end)

    # H2 generation limits — dispatchable
    @constraint(SP_models[w], cMaxH2GenVar[h in H_dis, t in T],
        vH2Gen[h, t] .- eAvailH2GenCap[h] * hsc_gen[h, :max_op_level] *
        hsc_gen_var[H_w[w][t], hsc_gen[h, :resource]] <= 0)
    @constraint(SP_models[w], cMinH2GenVar[h in H_dis, t in T],
        eAvailH2GenCap[h] * hsc_gen[h, :min_op_level] *
        hsc_gen_var[H_w[w][t], hsc_gen[h, :resource]] .- vH2Gen[h, t] <= 0)

    # H2 thermal unit commitment
    @constraint(SP_models[w], cH2OnlineUnits[h in H_ther, t in T],
        vH2GenOnline[h, t] .- eAvailH2GenUnit[h] <= 0)
    @constraint(SP_models[w], cH2StartLimits[h in H_ther, t in T],
        vH2GenStart[h, t] .- eAvailH2GenUnit[h] .+ vH2GenOnline[h, t] <= 0)
    @constraint(SP_models[w], cMaxH2GenTher[h in H_ther, t in T],
        vH2Gen[h, t] .- hsc_gen[h, :max_op_level] * hsc_gen[h, :rep_capacity] *
        vH2GenOnline[h, t] * hsc_gen_var[H_w[w][t], hsc_gen[h, :resource]] <= 0)
    @constraint(SP_models[w], cMinH2GenTher[h in H_ther, t in T],
        hsc_gen[h, :min_op_level] * hsc_gen[h, :rep_capacity] *
        vH2GenOnline[h, t] * hsc_gen_var[H_w[w][t], hsc_gen[h, :resource]] .- vH2Gen[h, t] <= 0)
    @constraint(SP_models[w], cH2ShutLimits[h in H_ther, t in T],
        vH2GenShut[h, t] .- vH2GenOnline[h, t] <= 0)

    # H2 thermal unit commitment (cyclic)
    @constraint(SP_models[w], cH2UnitOnlineCon[h in H_ther, t in 2:length(T)],
        vH2GenOnline[h, t] .- vH2GenOnline[h, t-1] == vH2GenStart[h, t] .- vH2GenShut[h, t])
    @constraint(SP_models[w], cH2UnitOnlineConFirst[h in H_ther],
        vH2GenOnline[h, 1] .- vH2GenOnlineFirst[h] == vH2GenStart[h, 1] .- vH2GenShut[h, 1])
    @constraint(SP_models[w], cH2GenOnlineConCycle[h in H_ther],
        vH2GenOnlineFirst[h] == vH2GenOnline[h, 168])

    # Minimum up/down time — H2 thermal
    @constraint(SP_models[w], cMinUpTimeH2Gen[h in H_ther, t in T],
        sum(vH2GenStart[h, tt] for tt in intersect(T, (t - hsc_gen[h, :up_time]):t)) .-
        vH2GenOnline[h, t] <= 0)
    @constraint(SP_models[w], cMinDnTimeH2Gen[h in H_ther, t in T],
        sum(vH2GenShut[h, tt] for tt in intersect(T, (t - hsc_gen[h, :down_time]):t)) .-
        eAvailH2GenUnit[h] .+ vH2GenOnline[h, t] <= 0)

    # Ramp constraints — H2 dispatchable
    @constraint(SP_models[w], cH2GenCycle[g in H_dis],
        vH2GenFirst[g] == vH2Gen[g, 168])
    @constraint(SP_models[w], cH2GenRampUp[g in H_dis, t in 2:length(T)],
        vH2Gen[g, t] - vH2Gen[g, t-1] .- hsc_gen[g, :ramp_up_percentage] * eAvailH2GenCap[g] <= 0)
    @constraint(SP_models[w], cH2GenRampUpFirst[g in H_dis],
        vH2Gen[g, 1] - vH2GenFirst[g] .- hsc_gen[g, :ramp_up_percentage] * eAvailH2GenCap[g] <= 0)
    @constraint(SP_models[w], cH2GenRampDn[g in H_dis, t in 2:length(T)],
        vH2Gen[g, t-1] - vH2Gen[g, t] .- hsc_gen[g, :ramp_down_percentage] * eAvailH2GenCap[g] <= 0)
    @constraint(SP_models[w], cH2GenRampDnFirst[g in H_dis],
        vH2GenFirst[g] - vH2Gen[g, 1] .- hsc_gen[g, :ramp_down_percentage] * eAvailH2GenCap[g] <= 0)

    # Ramp constraints — H2 thermal
    @constraint(SP_models[w], cTherH2GenRampDn[h in H_ther, t in 2:length(T)],
        vH2Gen[h, t-1] .- vH2Gen[h, t] .-
        hsc_gen[h, :rep_capacity] * hsc_gen[h, :ramp_down_percentage] * (vH2GenOnline[h, t] .- vH2GenStart[h, t]) .+
        hsc_gen[h, :rep_capacity] * vH2GenStart[h, t] * hsc_gen[h, :min_op_level] .-
        vH2GenShut[h, t] * hsc_gen[h, :rep_capacity] *
        min(hsc_gen_var[H_w[w][t], hsc_gen[h, :resource]],
            max(hsc_gen[h, :min_op_level], hsc_gen[h, :ramp_down_percentage])) <= 0)
    @constraint(SP_models[w], cTherH2GenRampDnFirst[h in H_ther],
        vH2GenFirst[h] .- vH2Gen[h, 1] .-
        hsc_gen[h, :rep_capacity] * hsc_gen[h, :ramp_down_percentage] * (vH2GenOnline[h, 1] .- vH2GenStart[h, 1]) .+
        hsc_gen[h, :rep_capacity] * vH2GenStart[h, 1] * hsc_gen[h, :min_op_level] .-
        vH2GenShut[h, 1] * hsc_gen[h, :rep_capacity] *
        min(hsc_gen_var[H_w[w][1], hsc_gen[h, :resource]],
            max(hsc_gen[h, :min_op_level], hsc_gen[h, :ramp_down_percentage])) <= 0)
    @constraint(SP_models[w], cTherH2GenRampUp[h in H_ther, t in 2:length(T)],
        vH2Gen[h, t] .- vH2Gen[h, t-1] .-
        hsc_gen[h, :rep_capacity] * hsc_gen[h, :ramp_up_percentage] * (vH2GenOnline[h, t] .- vH2GenStart[h, t]) .+
        hsc_gen[h, :rep_capacity] * vH2GenShut[h, t]  * hsc_gen[h, :min_op_level] .-
        vH2GenStart[h, t] * hsc_gen[h, :rep_capacity] *
        min(hsc_gen_var[H_w[w][t], hsc_gen[h, :resource]],
            max(hsc_gen[h, :min_op_level], hsc_gen[h, :ramp_up_percentage])) <= 0)
    @constraint(SP_models[w], cTherH2GenRampUpFirst[g in H_ther],
        vH2Gen[g, 1] .- vH2GenFirst[g] .-
        hsc_gen[g, :rep_capacity] * hsc_gen[g, :ramp_up_percentage] * (vH2GenOnline[g, 1] .- vH2GenStart[g, 1]) .-
        hsc_gen[g, :rep_capacity] * vH2GenShut[g, 1]  * hsc_gen[g, :min_op_level] .+
        vH2GenStart[g, 1] * hsc_gen[g, :rep_capacity] *
        min(hsc_gen_var[H_w[w][1], hsc_gen[g, :resource]],
            max(hsc_gen[g, :min_op_level], hsc_gen[g, :ramp_up_percentage])) <= 0)

    # H2 storage
    @constraint(SP_models[w], cH2StoBalance[s in Q, t in 2:length(T)],
        vH2StoSOC[s, t] == (1 - hsc_gen[s, :etta_self_dis]) * vH2StoSOC[s, t-1] +
        vH2StoCha[s, t] * hsc_gen[s, :etta_cha] - (1/hsc_gen[s, :etta_dis]) * vH2StoDis[s, t])
    @constraint(SP_models[w], cH2StoBalanceFirst[s in Q],
        vH2StoSOC[s, 1] == (1 - hsc_gen[s, :etta_self_dis]) * eH2SOCFirst[s] +
        vH2StoCha[s, 1] * hsc_gen[s, :etta_cha] - (1/hsc_gen[s, :etta_dis]) * vH2StoDis[s, 1] +
        slack_start[s] - slack_start_neg[s])
    @constraint(SP_models[w], cMaxH2StoSOC[s in Q, t in T],
        vH2StoSOC[s, t] .- hsc_gen[s, :h2stor_max_level] * eAvailH2StoCap[s] <= 0)
    @constraint(SP_models[w], cMinH2StoSOC[s in Q, t in T],
        hsc_gen[s, :h2stor_min_level] * eAvailH2StoCap[s] .- vH2StoSOC[s, t] <= 0)
    @constraint(SP_models[w], cMaxH2StoChar[s in Q, t in T],
        vH2StoCha[s, t] .- eAvailH2StoCompCap[s] <= 0)
    @constraint(SP_models[w], cMaxH2StoDis[s in Q, t in 2:length(T)],
        vH2StoDis[s, t] .- hsc_gen[s, :etta_dis] * vH2StoSOC[s, t-1] <= 0)
    @constraint(SP_models[w], cH2StoSOCLast[s in Q],
        eH2SOCLast[s] .+ slack_end[s] == vH2StoSOC[s, 168])

    # H2 transmission
    @constraints(SP_models[w], begin
        cMaxH2PipeFlowOut[i in I, t in T],
            vH2FlowPos[i, t] .- eAvailH2Pipe[i] * hsc_pipelines[i, :max_op_level] *
            hsc_pipelines[i, :max_pipe_cap_tonne] <= 0
        cMaxH2PipeFlowIn[i in I, t in T],
            vH2FlowNeg[i, t] .- eAvailH2Pipe[i] * hsc_pipelines[i, :max_op_level] *
            hsc_pipelines[i, :max_pipe_cap_tonne] <= 0
    end)
    @constraints(SP_models[w], begin
        cMaxH2PipeFlowOutComp[i in I, t in T],
            vH2FlowPos[i, t] .- eAvailH2PipeCompCap[i] <= 0
        cMaxH2PipeFlowInComp[i in I, t in T],
            vH2FlowNeg[i, t] .- eAvailH2PipeCompCap[i] <= 0
    end)

    # Curtailment / emission policy
    @constraint(SP_models[w], cPowCrt[z in Z, t in T],
        vPowCrt[z, t] .-
        sum(vPowGen[g, t] * (pow_gen[g, :zone] == z ? 1 : 0) for g in G) <= 0)
    @constraint(SP_models[w], cEmissionCapByWeek,
        eEmissionByWeek .- vExtraEmission .- eMaxEmissionByWeek <= 0)
end
