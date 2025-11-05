function write_costs_files(datadir,
    SP_models,
    pow_gen, hsc_gen,
    G, S, H, Q,
    PowGen_vals, PowStoCha_vals,
    H2Gen_vals, H2StoCha_vals,
    W, T
)
    results_dir = joinpath(datadir, "Results")
    # Total timesteps
    timesteps = [(w,t) for w in W for t in T]

    # ─────────────────────────────────────────────────────────────────────
    # 1) POWER COSTS
    # ─────────────────────────────────────────────────────────────────────
    power_cost_attrs = [
      "Inv_Cost",
      "FOM_Cost",
      "VOM_Cost",
      "Fuel_Cost",
      "Start_Cost"
    ]
    df_pc = DataFrame(Attribute = power_cost_attrs)

    for g in sort(collect(G))
        name   = pow_gen[g, :resource]
        caprep = pow_gen[g, :rep_capacity]
        invc   = pow_gen[g, :inv_cost_per_mwyr]
        fomc   = pow_gen[g, :fom_cost_per_mwyr]
        vomc   = pow_gen[g, :vom_cost_mwh]
        fuelc  = pow_gen[g, :heat_rate_mmbtu_per_mwh]*fuel_costs[pow_gen[g, :fuel]][1]
        startc = pow_gen[g, :start_cost_per_mw]

        # new & total capacities
        new_units = value(vNewPowGenCap[g])
        ret_units = value(vRetPowGenCap[g])
        existing  = pow_gen[g, :existing_cap] / caprep
        tot_units = existing + new_units - ret_units
        tot_cap   = tot_units * caprep
        new_cap   = new_units * caprep

        # 1) Inv_Cost = inv_cost * new_cap
        row1 = invc * new_cap

        # 2) FOM_Cost = fom_cost * tot_cap
        row2 = fomc * tot_cap

        # 3) VOM_Cost = vom_cost * (sum over PowGen_vals[g, w, t])
        gen_sum = 0.0
        for (w,t) in timesteps
            gen_sum += PowGen_vals[(g, w, t)]
        end
        row3 = vomc * gen_sum

        # 4) Fuel_Cost = fuel_cost * gen_sum
        row4 = fuelc * gen_sum
        
        start_sum = 0.0
        # 5) Start_Cost = placeholder 0.0 (replace with startc * num_starts if tracked)
        if g in G_ther
            for w in W, t in T
                start_sum += value(SP_models[w][:vPowGenStart][g, t])
            end
        else 
            start_sum = 0.0
        end
        row5 = startc * start_sum

        df_pc[!, Symbol(name)] = [row1, row2, row3, row4, row5]
    end

    for s in sort(collect(S))
        name   = pow_gen[s, :resource]
        caprep = pow_gen[s, :rep_capacity]
        invc   = pow_gen[s, :inv_cost_per_mwhyr]
        fomc   = pow_gen[s, :fom_cost_per_mwhyr]
        vomc   = pow_gen[s, :vom_cost_mwh_charge]
        fuelc  = 0
        startc = 0

        new_units = value(vNewPowStoCap[s])
        ret_units = value(vRetPowStoCap[s])
        existing  = pow_gen[s, :existing_cap_mwh] / caprep
        tot_units = existing + new_units - ret_units
        tot_cap   = tot_units * caprep
        new_cap   = new_units * caprep

        row1 = invc * new_cap
        row2 = fomc * tot_cap

        # VOM_Cost for storage = vom * (sum of charge)
        charge_sum = 0.0
        for (w,t) in timesteps
            charge_sum += PowStoCha_vals[(s, w, t)]
        end
        row3 = vomc * charge_sum

        # Fuel_Cost = 0 (no fuel for storage)
        row4 = 0.0

        row5 = 0.0

        df_pc[!, Symbol(name)] = [row1, row2, row3, row4, row5]
    end

    CSV.write(joinpath(results_dir,"21_Costs_Power.csv"), df_pc)

    # ─────────────────────────────────────────────────────────────────────
    # 2) H₂ COSTS
    # ─────────────────────────────────────────────────────────────────────
    h2_cost_attrs = [
      "Inv_Cost",
      "FOM_Cost",
      "VOM_Cost",
      "Fuel_Cost",
      "Start_Cost"
    ]
    df_hc = DataFrame(Attribute = h2_cost_attrs)

    for h in sort(collect(H))
        name   = hsc_gen[h, :resource]
        caprep = hsc_gen[h, :rep_capacity]
        invc   = hsc_gen[h, :inv_cost_tonne_hr_p_yr]
        fomc   = hsc_gen[h, :fom_cost_p_tonne_p_hr_yr]
        vomc   = hsc_gen[h, :vom_cost_p_tonne]
        fuelc  = hsc_gen[h, :heat_rate_mmbtu_p_tonne]*fuel_costs[hsc_gen[h, :fuel]][1]
        startc = hsc_gen[h, :startup_cost_p_tonne_hr]

        new_units = value(vNewH2GenCap[h])
        ret_units = value(vRetH2GenCap[h])
        existing  = hsc_gen[h, :existing_cap_tonne_p_hr] / caprep
        tot_units = existing + new_units - ret_units
        tot_cap   = tot_units * caprep
        new_cap   = new_units * caprep

        row1 = invc * new_cap
        row2 = fomc * tot_cap

        # VOM_Cost for H₂ gen = vom * (sum of H2Gen)
        gen_sum = 0.0
        for (w,t) in timesteps
            gen_sum += H2Gen_vals[(h, w, t)]
        end
        row3 = vomc * gen_sum

        # Fuel_Cost = fuelc * gen_sum
        row4 = fuelc * gen_sum

        row5 = startc * tot_cap

        df_hc[!, Symbol(name)] = [row1, row2, row3, row4, row5]
    end

    for q in sort(collect(Q))
        # Storage
        name   = hsc_gen[q, :resource]
        caprep = hsc_gen[q, :rep_capacity]
        invc   = hsc_gen[q, :inv_cost_tonne_p_yr]
        fomc   = hsc_gen[q, :fom_cost_p_tonne_p_yr]
        vomc   = hsc_gen[q, :vom_cost_p_tonne]
        fuelc  = 0
        startc = 0

        new_units = value(vNewH2StoCap[q])
        ret_units = value(vRetH2StoCap[q])
        existing  = hsc_gen[q, :existing_cap_tonne] / caprep
        tot_units = existing + new_units - ret_units
        tot_cap   = tot_units * caprep
        new_cap   = new_units * caprep

        row1 = invc * new_cap
        row2 = fomc * tot_cap

        # VOM_Cost for H₂ storage = vom * (sum of charge)
        charge_sum = 0.0
        for (w,t) in timesteps
            charge_sum += H2StoCha_vals[(q, w, t)]
        end
        row3 = vomc * charge_sum

        row4 = 0.0
        row5 = 0.0

        df_hc[!, Symbol(name)] = [row1, row2, row3, row4, row5]
    end

    for q in sort(collect(Q))
        # Storage compressor
        name   = string(hsc_gen[q, :resource], "_Compressor")
        invc   = hsc_gen[q, :inv_cost_comp_tonne_hr_p_yr]
        fomc   = hsc_gen[q, :fom_cost_comp_tonne_hr_p_yr]
        vomc   = 0

        new_cap = value(vNewH2StoCompCap[q])
        ret_cap = value(vRetH2StoCompCap[q])
        existing  = hsc_gen[q, :existing_cap_comp_tonne_hr]
        tot_cap = existing + new_cap - ret_cap

        row1 = invc * new_cap
        row2 = fomc * tot_cap

        # VOM_Cost for H₂ sto compressor: vom * (sum of H2StoCha)
        charge_sum = 0.0
        for (w,t) in timesteps
            charge_sum += H2StoCha_vals[(q, w, t)]
        end
        row3 = vomc * charge_sum

        row4 = 0.0
        row5 = 0.0

        df_hc[!, Symbol(name)] = [row1, row2, row3, row4, row5]
    end

    CSV.write(joinpath(results_dir, "22_Costs_H2.csv"), df_hc)
end

function write_line_costs_power(datadir,pow_lines, L)
    results_dir = joinpath(datadir, "Results")
    attrs = ["Inv_Cost"]
    df = DataFrame(Attribute = attrs)
    
    for l in sort(collect(L))
        invc   = pow_lines[l, :line_reinforcement_cost_per_mwyr]   
        new_cap = value(vNewPowTraCap[l]) 
    
        # Investment cost = inv_cost * new_cap
        df[!, Symbol(string(l))] = [invc * new_cap]
    end
    
    CSV.write(joinpath(results_dir, "23_Power_Line_Costs.csv"), df)

end

function write_h2_pipe_costs(datadir,
    hsc_pipelines, I,
    H2FlowPos_vals, H2FlowNeg_vals,
    W, T
    )
    timesteps = [(w,t) for w in W for t in T]
    results_dir = joinpath(datadir, "Results")
    attrs = [
      "Investment Cost",
      "FOM Cost",
      "VOM Cost",
      "Inv Cost Comp",
      "FOM Cost Comp",
      "VOM Cost Comp"
    ]
    df = DataFrame(Attribute = attrs)
    
    for i in sort(collect(I))
        
        inv_pipe   = hsc_pipelines[i, :investment_cost_per_length]*hsc_pipelines[i, :distance]  
        fom_pipe   = hsc_pipelines[i, :fom_per_length]*hsc_pipelines[i, :distance]
        vom_pipe   = hsc_pipelines[i, :vom_per_tonne]  
    
        new_pipe_units = value(vNewH2Pipe[i])
        ret_pipe_units = value(vRetH2Pipe[i])
        existing_pipes = hsc_pipelines[i, :existing_num_pipes]
        total_pipe_units = existing_pipes + new_pipe_units - ret_pipe_units
    
        # (1) Inv_Cost_Pipe = inv_pipe * new_pipe_units
        row1 = inv_pipe * new_pipe_units
    
        # (2) FOM_Cost_Pipe = fom_pipe * total_pipe_units
        row2 = fom_pipe * total_pipe_units
    
        # (3) VOM_Cost_Pipe = vom_pipe * sum_over_all_hours(net_flow)
        net_sum = 0.0
        for (w,t) in timesteps
            net_sum += (H2FlowPos_vals[(i,w,t)] + H2FlowNeg_vals[(i,w,t)])
        end
        row3 = vom_pipe * net_sum
    
        # COMPRESSOR side
        inv_comp   = hsc_pipelines[i, :compressor_inv_per_length]*hsc_pipelines[i, :distance]
        fom_comp   = 0
        
    
        new_comp_units = value(vNewH2PipeCompCap[i])
        ret_comp_units = value(vRetH2PipeCompCap[i])
        existing_comp  = hsc_pipelines[i, :existing_comp_cap_tonne_hr]
        total_comp_units = existing_comp + new_comp_units - ret_comp_units
    
        # (4) Inv_Cost_Comp = inv_comp * new_comp_units
        row4 = inv_comp * new_comp_units
    
        # (5) FOM_Cost_Comp = fom_comp * total_comp_units
        row5 = fom_comp * total_comp_units
    
        row6 = 0
    
        df[!, Symbol(string(i))] = [row1, row2, row3, row4, row5, row6]
    end
    
    CSV.write(joinpath(results_dir, "24_H2_Pipe_Costs.csv"), df)
    
end