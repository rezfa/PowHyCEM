function write_LCOH(datadir,
    SP_models, hsc_gen, H, Q,
    H2Gen_vals, H2StoCha_vals,
    H2FlowPos_vals, H2FlowNeg_vals,
    W, T,
    fuel_costs::Dict{<:AbstractString, <:Any}
)
    results_dir = joinpath(datadir, "Results")
    # 1) Precompute all (w,t) pairs
    timesteps = [(w, t) for w in W for t in T]

    # 2) Prepare category accumulators
    cats = [:Green, :Blue, :Grey]
    gen_tot        = Dict(c => 0.0 for c in cats)
    inv_cost_gen   = Dict(c => 0.0 for c in cats)
    fom_cost_gen   = Dict(c => 0.0 for c in cats)
    vom_cost_gen   = Dict(c => 0.0 for c in cats)
    fuel_cost_gen  = Dict(c => 0.0 for c in cats)
    start_cost_gen = Dict(c => 0.0 for c in cats)

    # 3) Classification function for H₂ generators
    classify_gen = function(resname::AbstractString)
        if occursin("SMR_CCS", resname)
            return :Blue
        elseif occursin("AEL", resname) || occursin("PEM", resname)
            return :Green
        else
            return :Grey
        end
    end

    # 4) Loop over all H₂ generators
    for h in sort(collect(H))
        resname = String(hsc_gen[h, :resource])
        cat = classify_gen(resname)

        caprep = hsc_gen[h, :rep_capacity]
        invc   = hsc_gen[h, :inv_cost_tonne_hr_p_yr]
        fomc   = hsc_gen[h, :fom_cost_p_tonne_p_hr_yr]
        vomc   = hsc_gen[h, :vom_cost_p_tonne]
        fuel_name   = String(hsc_gen[h, :fuel])
        mmbtu_price = fuel_costs[fuel_name][1]
        heat_rate   = hsc_gen[h, :heat_rate_mmbtu_p_tonne]
        fuelc  = heat_rate * mmbtu_price
        startc = hsc_gen[h, :startup_cost_p_tonne_hr]

        # New & total capacity for h (tonne/hr)
        new_units = value(vNewH2GenCap[h])
        ret_units = value(vRetH2GenCap[h])
        existing  = hsc_gen[h, :existing_cap_tonne_p_hr] / caprep
        tot_units = existing + new_units - ret_units
        tot_cap   = tot_units * caprep
        new_cap   = new_units * caprep

        # (1) Total generation by h
        gen_sum = 0.0
        for (w, t) in timesteps
            gen_sum += H2Gen_vals[(h, w, t)]
        end
        gen_tot[cat] += gen_sum

        # (2) Gen_Inv_Cost
        inv_cost_gen[cat] += invc * new_cap

        # (3) Gen_FOM_Cost
        fom_cost_gen[cat] += fomc * tot_cap

        # (4) Gen_VOM_Cost
        vom_cost_gen[cat] += vomc * gen_sum

        # (5) Gen_Fuel_Cost
        fuel_cost_gen[cat] += fuelc * gen_sum

        # (6) Gen_Start_Cost (only for non‐Green)
        start_sum = 0.0
        if cat != :Green
            for w in W
                let var_h2start = SP_models[w][:vH2GenStart]
                    ax1 = axes(var_h2start, 1)
                    ax2 = axes(var_h2start, 2)
                    if h ∈ ax1
                        for t in T
                            if t ∈ ax2
                                start_sum += value(var_h2start[h, t])
                            end
                        end
                    end
                end
            end
        end
        start_cost_gen[cat] += startc * start_sum
    end

    # 5) Compute total storage cost (goes only into "All")
    storage_cost = 0.0
    for q in sort(collect(Q))
        caprep = hsc_gen[q, :rep_capacity]
        invc   = hsc_gen[q, :inv_cost_tonne_p_yr]
        fomc   = hsc_gen[q, :fom_cost_p_tonne_p_hr_yr]
        vomc   = hsc_gen[q, :vom_cost_p_tonne]

        new_cap = value(vNewH2StoCap[q])
        ret_cap = value(vRetH2StoCap[q])
        existing  = hsc_gen[q, :existing_cap_tonne]
        tot_cap = existing + new_cap - ret_cap
        

        storage_cost += invc * new_cap
        storage_cost += fomc * tot_cap

        charge_sum = 0.0
        for (w, t) in timesteps
            charge_sum += H2StoCha_vals[(q, w, t)]
        end
        storage_cost += vomc * charge_sum
    end

    # 6) Compute pipeline + compressor cost (goes only into "All")
    pipeline_cost = 0.0
    pipe_ids = unique(first.(keys(H2FlowPos_vals)))
    for i in sort(pipe_ids)
        # PIPE side: use per‐length cost × distance
        inv_pipe    = hsc_pipelines[i, :investment_cost_per_length] * hsc_pipelines[i, :distance]
        vom_pipe    = hsc_pipelines[i, :vom_per_tonne]

        new_pipe_units = value(vNewH2Pipe[i])
        ret_pipe_units = value(vRetH2Pipe[i])
        existing_pipes = hsc_pipelines[i, :existing_num_pipes]
        total_pipe_units = existing_pipes + new_pipe_units - ret_pipe_units

        pipeline_cost += inv_pipe * new_pipe_units

        net_sum = 0.0
        for (w, t) in timesteps
            net_sum += (H2FlowPos_vals[(i, w, t)] - H2FlowNeg_vals[(i, w, t)])
        end
        pipeline_cost += vom_pipe * net_sum

        # COMPRESSOR side: cost per length × distance
        inv_comp = hsc_pipelines[i, :compressor_inv_per_length] * hsc_pipelines[i, :distance]
        fom_comp = hsc_pipelines[i, :fom_comp_p_tonne_hr]          * hsc_pipelines[i, :distance]

        new_comp_units = value(vNewH2PipeCompCap[i])
        ret_comp_units = value(vRetH2PipeCompCap[i])
        existing_comp  = hsc_pipelines[i, :existing_comp_cap_tonne_hr]
        total_comp_units = existing_comp + new_comp_units - ret_comp_units

        pipeline_cost += inv_comp * new_comp_units
        pipeline_cost += fom_comp * total_comp_units

    end

    # 7) Build DataFrame with 8 cost rows
    attrs = [
      "Total_Generation",
      "Gen_Inv_Cost",
      "Gen_FOM_Cost",
      "Gen_VOM_Cost",
      "Gen_Fuel_Cost",
      "Gen_Start_Cost",
      "Storage_Cost",
      "Pipeline_Cost",
      "Emission_Cost"
    ]

    df = DataFrame(
      Attribute = attrs,
      Green_H2  = zeros(Float64, length(attrs)),
      Blue_H2   = zeros(Float64, length(attrs)),
      Grey_H2   = zeros(Float64, length(attrs)),
      All       = zeros(Float64, length(attrs))
    )

    df[1, :Green_H2] = gen_tot[:Green]
    df[1, :Blue_H2]  = gen_tot[:Blue]
    df[1, :Grey_H2]  = gen_tot[:Grey]
    df[1, :All]      = sum(values(gen_tot))

    df[2, :Green_H2] = inv_cost_gen[:Green]
    df[2, :Blue_H2]  = inv_cost_gen[:Blue]
    df[2, :Grey_H2]  = inv_cost_gen[:Grey]
    df[2, :All]      = sum(values(inv_cost_gen))

    df[3, :Green_H2] = fom_cost_gen[:Green]
    df[3, :Blue_H2]  = fom_cost_gen[:Blue]
    df[3, :Grey_H2]  = fom_cost_gen[:Grey]
    df[3, :All]      = sum(values(fom_cost_gen))

    df[4, :Green_H2] = vom_cost_gen[:Green]
    df[4, :Blue_H2]  = vom_cost_gen[:Blue]
    df[4, :Grey_H2]  = vom_cost_gen[:Grey]
    df[4, :All]      = sum(values(vom_cost_gen))

    df[5, :Green_H2] = fuel_cost_gen[:Green]
    df[5, :Blue_H2]  = fuel_cost_gen[:Blue]
    df[5, :Grey_H2]  = fuel_cost_gen[:Grey]
    df[5, :All]      = sum(values(fuel_cost_gen))

    df[6, :Green_H2] = start_cost_gen[:Green]
    df[6, :Blue_H2]  = start_cost_gen[:Blue]
    df[6, :Grey_H2]  = start_cost_gen[:Grey]
    df[6, :All]      = sum(values(start_cost_gen))

    df[7, :All] = storage_cost
    df[8, :All] = pipeline_cost

    em_cost_weeks = Float64[]

    for w in W
        # 1) H₂ from gen (all numeric)
        h2_gen = sum(
          H2Gen_vals[(h,w,t)] *
            hsc_gen[h, :heat_rate_mmbtu_p_tonne] *
            CO2_content[hsc_gen[h, :fuel]][1]
          for h in H, t in T
        )
    
        # 2) H₂ from starts (value(...) turns vH2GenStart into Float64)
        h2_start = sum(
          value(SP_models[w][:vH2GenStart][h,t]) *
            hsc_gen[h, :heat_rate_mmbtu_p_tonne] *
            CO2_content[hsc_gen[h, :fuel]][1]
          for h in H_ther, t in T
        )
    
        h2_emi = h2_gen + h2_start         # Float64
    
        # 3) total CO₂ that week (value(...) → Float64)
        tot_emi = value(SP_models[w][:eEmissionByWeek])
    
        # 4) fraction due to H₂
        frac_H2 = tot_emi > 0 ? (h2_emi / tot_emi) : 0.0
    

    
        # 6) CO₂ price (if it’s a JuMP var or expr, call value; 
        #    if it’s a pure constant you can skip value)
        price   = value(SP_models[w][:eEmissionCost])
    
        # 7) allocate and push a pure Float64
        push!(em_cost_weeks, frac_H2 * price)
        
    end
    

    emission_cost = sum(em_cost_weeks)
    
    df[9, :All] = emission_cost    # Emission_Cost row

    # Allow missing values in the three category columns plus "All" for the LCOH row
    for col in [:Green_H2, :Blue_H2, :Grey_H2, :All]
        allowmissing!(df, col)
    end

    # 8) Compute LCOH = (sum of all cost rows) / (total_generation * 1000)
    total_cost = sum(df[2:9, :All])
    total_gen  = df[1, :All]  # in tonnes
    lcoh       = total_cost / (total_gen * 1000)  # €/kg

    # Append LCOH row (missing for category columns, actual value in "All")
    push!(df, (
        "LCOH_EUR_per_kg",
        missing,
        missing,
        missing,
        lcoh
    ))

    # 9) Write to CSV
    CSV.write(joinpath(results_dir,"25_LCOH_Detail.csv"), df)

end