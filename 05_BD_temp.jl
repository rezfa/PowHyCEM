using JuMP, Gurobi, DataFrames, CSV, Plots


#function run_benders(datadir = joinpath("/Users/rez/Documents/Engineering/Coding/Julia/PowHyCEM/Input_Data"))
    
    
    datadir = joinpath("/Users/rez/Documents/Engineering/Coding/Julia/PowHyCEM/Input_Data")
 # ---------- 1.  raw tables ------------------------------------------------
    pow_gen       = CSV.read(joinpath(datadir, "Powe_Gen_Data.csv"),  DataFrame)
    pow_gen_var   = CSV.read(joinpath(datadir, "Pow_Gen_Var.csv"),    DataFrame)
    pow_load      = CSV.read(joinpath(datadir, "Pow_Load.csv"),       DataFrame)
    pow_lines     = CSV.read(joinpath(datadir, "Pow_Network.csv"),    DataFrame)

    hsc_gen       = CSV.read(joinpath(datadir, "HSC_Gen_Data.csv"),   DataFrame)
    hsc_gen_var   = CSV.read(joinpath(datadir, "HSC_Gen_Var.csv"),    DataFrame)
    hsc_load      = CSV.read(joinpath(datadir, "HSC_load.csv"),       DataFrame)
    hsc_pipelines = CSV.read(joinpath(datadir, "HSC_Pipelines.csv"),  DataFrame)

    fuel          = CSV.read(joinpath(datadir, "Fuels_data.csv"),     DataFrame)
    zones         = CSV.read(joinpath(datadir, "Zone_Char.csv"),      DataFrame)

    # ---------- 2.  make headers lower-case ----------------------------------
    for f in (pow_gen, pow_load, pow_lines,
              hsc_gen, hsc_load, hsc_pipelines,
              fuel,    zones)
        rename!(f, lowercase.(names(f)))
    end

    # ---------- 3.  convenience data frames ----------------------------------
    dfGen    = pow_gen[pow_gen.gen_type .>= 1, :]
    dfStor   = pow_gen[pow_gen.stor_type .> 0, :]
    dfH2Gen  = hsc_gen[hsc_gen.h_gen_type .>= 1, :]
    dfH2Stor = hsc_gen[hsc_gen.h_stor    .>= 1, :]

    # ---------- 4.  core sets -------------------------------------------------
    G      = dfGen.r_id
    S      = dfStor.r_id
    Z      = zones.zones
    W      = 1:52
    T      = 1:168
    hours_per_week = length(T)
    H_w    = [((w-1)*hours_per_week+1):(w*hours_per_week) for w in W]

    G_ther = dfGen[dfGen.gen_type .== 1, :r_id]
    G_ren  = dfGen[dfGen.gen_type .>  1, :r_id]
    V      = dfGen[dfGen.gen_type .== 3, :r_id]
    L      = pow_lines.network_lines

    H      = dfH2Gen[dfH2Gen.h_gen_type .>= 1, :r_id]
    Q      = dfH2Stor.r_id
    I      = hsc_pipelines.hsc_pipelines
    H_ther = hsc_gen[hsc_gen.h_gen_type .== 1, :r_id]
    H_dis  = hsc_gen[hsc_gen.h_gen_type .== 2, :r_id]

    # ---------- 5.  weekly demand tensors ------------------------------------
    start_p = findfirst(==("load_mw_z1"), names(pow_load))
    pow_demand = Dict{Int,Matrix}()
    for w in W
        pow_demand[w] = Matrix(pow_load[H_w[w], start_p:start_p+length(Z)-1])
    end

    start_h  = findfirst(==("load_hsc_tonne_z1"), names(hsc_load))
    h2_demand = Dict{Int,Matrix}()
    for w in W
        h2_demand[w] = Matrix(hsc_load[H_w[w], start_h:start_h+length(Z)-1])
    end

    # ---------- 6.  fuel price & emission dictionaries -----------------------
    fuels        = names(fuel)[2:end]
    costs        = Matrix(fuel[2:end, 2:end])
    CO2_content  = fuel[1, 2:end]              # ton CO₂ / MMBtu
    fuel_costs = Dict(fuels[i] => costs[:, i]       for i in eachindex(fuels))
    fuel_CO2   = Dict(fuels[i] => CO2_content[i]    for i in eachindex(fuels))

    # ---------- 7.  network incidence matrices -------------------------------
    col_p       = findfirst(==("z1"), names(pow_lines))
    Pow_Network = pow_lines[1:length(L), col_p:col_p+length(Z)-1]

    col_h       = findfirst(==("z1"), names(hsc_pipelines))
    H2_Network  = hsc_pipelines[1:length(I), col_h:col_h+length(Z)-1]

    # ---------- 8.  bundle & return ------------------------------------------
    
    iters    = Int[]
    LB_hist  = Float64[]
    UB_hist  = Float64[]
    gap_hist = Float64[]

    # Create the empty figure and the three series (legend only here):
    plt = plot(size=(1000,600),
           xlabel="Iteration",
           ylabel="Cost (×10⁶)",
           legend=:topright,
           ylim=(0,300))

    # initialize the series (empty for now)
    plot!(plt, iters, LB_hist,  label="LB",  color=:blue)
    plot!(plt, iters, UB_hist,  label="UB",  color=:red)
    #plot!(plt, iters, gap_hist, label="Gap", color=:green, linestyle=:dash)


    #Defining the Master Problem
    MP = Model(Gurobi.Optimizer)
    set_optimizer_attribute(MP, "Method", 2)      
    #set_optimizer_attribute(MP, "Crossover", 0)
    set_optimizer_attribute(MP,"MIPGap",1e-3)
    set_optimizer_attribute(MP, "OutputFlag", 0)
    set_optimizer_attribute(MP, "LogToConsole", 0)

    
    # ---- investment Variables ---- #
    @variable(MP, vNewPowGenCap[g in G]>=0, Int)
    @variable(MP, vRetPowGenCap[g in G]>=0, Int)
    @variable(MP, vNewPowStoCap[s in S]>=0, Int)
    @variable(MP, vRetPowStoCap[s in S]>=0, Int)
    @variable(MP, vNewPowTraCap[l in L]>=0, Int)
    @variable(MP, vNewH2GenCap[h in H]>=0, Int)
    @variable(MP, vRetH2GenCap[h in H]>=0, Int)
    @variable(MP, vNewH2StoCap[s in Q]>=0, Int)
    @variable(MP, vRetH2StoCap[s in Q]>=0, Int)
    @variable(MP, vNewH2Pipe[i in I]>=0, Int)
    @variable(MP, vRetH2Pipe[i in I]>=0, Int)
    @variable(MP, vNewH2PipeCompCap[i in I]>=0, Int)
    @variable(MP, vRetH2PipeCompCap[i in I]>=0, Int)
    @variable(MP, vNewH2StoCompCap[s in Q]>=0, Int)
    @variable(MP, vRetH2StoCompCap[s in Q]>=0, Int)   
    # ---- Emission Budget Variable ---- #
    @variable(MP, vMaxEmissionByWeek[w in W]>=0)
    # ---- Set of Cuts Variable ---- #
    @variable(MP, theta[w in W] >=0)

    # ---- MP Expressions ---- #
    @expression(MP, eTotPowGenCap[g in G], pow_gen[g, :existing_cap] .+ pow_gen[g, :rep_capacity]*(vNewPowGenCap[g] .- vRetPowGenCap[g]))
    @expression(MP, eTotPowStoCap[s in S], pow_gen[s, :existing_cap] .+ pow_gen[s, :rep_capacity]*(vNewPowStoCap[s] .- vRetPowStoCap[s]))
    @expression(MP, eTotH2GenCap[h in H], hsc_gen[h, :existing_cap_tonne_p_hr] .+ hsc_gen[h, :rep_capacity]*(vNewH2GenCap[h] .- vRetH2GenCap[h]))
    @expression(MP, eTotH2StoCap[s in Q], hsc_gen[s, :existing_cap_tonne] + hsc_gen[s, :rep_capacity]*(vNewH2StoCap[s] - vRetH2StoCap[s]))
    @expression(MP, eTotH2StoCompCap[s in Q], hsc_gen[s, :existing_cap_comp_tonne_hr] + vNewH2StoCompCap[s] - vRetH2StoCompCap[s])
    @expression(MP, eTotH2Pipe[i in I], hsc_pipelines[i, :existing_num_pipes] + vNewH2Pipe[i] - vRetH2Pipe[i])
    @expression(MP, ePowGenLandUse[z in Z], sum((vNewPowGenCap[g] - vRetPowGenCap[g])*pow_gen[g, :rep_capacity]*pow_gen[g, :land_use_km2_p_cap]*(pow_gen[g,:zone]==z ? 1 : 0) for g in G))
    @expression(MP, ePowStoLandUse[z in Z], sum((vNewPowStoCap[s] - vRetPowStoCap[s])*pow_gen[s, :rep_capacity]*pow_gen[s, :land_use_km2_p_cap]*(pow_gen[s,:zone]==z ? 1 : 0) for s in S))
    @expression(MP, eH2GenLandUse[z in Z], sum((vNewH2GenCap[h]-vRetH2GenCap[h])*hsc_gen[h, :rep_capacity]*hsc_gen[h, :land_use_km2_p_cap]*(hsc_gen[h,:zone]==z ? 1 : 0) for h in H))
    @expression(MP, eH2StoLandUse[z in Z], sum((vNewH2StoCap[s]-vRetH2StoCap[s])*hsc_gen[s, :rep_capacity]* hsc_gen[s, :land_use_km2_p_cap]*(hsc_gen[s,:zone]==z ? 1 : 0) for s in Q))
    @expression(MP, eH2PipeLandUse[z in Z],0.5 * sum(hsc_pipelines[i, :land_use_km2_p_length]*hsc_pipelines[i, :distance]*hsc_pipelines[i, :land_use_km2_p_length]*(vNewH2Pipe[i]-vRetH2Pipe[i])*abs(H2_Network[i,z]) for i in I))
    @expression(MP, eTotH2PipeCompCap[i in I], hsc_pipelines[i, :existing_comp_cap_tonne_hr] + vNewH2PipeCompCap[i] - vRetH2PipeCompCap[i])
    #Cost Expressions
    @expression(MP, eCostPowGenInv, sum(pow_gen[g, :inv_cost_per_mwyr] .* vNewPowGenCap[g] .* pow_gen[g, :rep_capacity] .+ pow_gen[g, :fom_cost_per_mwyr] .* eTotPowGenCap[g] for g in G))
    @expression(MP, eCostPowStoInv, 
        sum(pow_gen[s, :inv_cost_per_mwhyr] .* vNewPowStoCap[s] *pow_gen[s, :rep_capacity] .+ (pow_gen[s, :fom_cost_per_mwhyr] .* eTotPowStoCap[s]) for s in S)
    ) 
    @expression(MP, eCostPowTraInv, sum(pow_lines[l, :line_reinforcement_cost_per_mwyr] .* vNewPowTraCap[l] for l in L))
    @expression(MP, eCostH2GenInv, sum(hsc_gen[h, :inv_cost_tonne_hr_p_yr]*vNewH2GenCap[h]*hsc_gen[h, :rep_capacity] + hsc_gen[h, :fom_cost_p_tonne_p_hr_yr]*eTotH2GenCap[h] for h in H))
    @expression(MP, eCostH2StoInv, sum(hsc_gen[s, :inv_cost_tonne_p_yr]*vNewH2StoCap[s]*hsc_gen[s, :rep_capacity] + hsc_gen[s, :inv_cost_comp_tonne_hr_p_yr]*vNewH2StoCompCap[s] +
                                     hsc_gen[s, :fom_cost_p_tonne_p_yr]*eTotH2StoCap[s] + hsc_gen[s, :fom_cost_comp_tonne_hr_p_yr]*eTotH2StoCompCap[s] for s in Q)
    )
    @expression(MP, eCostH2TraInv, sum(hsc_pipelines[i, :investment_cost_per_length]*hsc_pipelines[i, :distance]*vNewH2Pipe[i] +
                                        hsc_pipelines[i, :fom_per_length]*hsc_pipelines[i, :distance]*eTotH2Pipe[i] +
                                        hsc_pipelines[i, :compressor_inv_per_length]*hsc_pipelines[i, :distance]*vNewH2PipeCompCap[i] + 
                                        hsc_pipelines[i, :fom_comp_p_tonne_hr]*eTotH2PipeCompCap[i] for i in I)
    )

    # ---- MP Objective ---- # 
    MP_obj = eCostPowGenInv .+ eCostPowStoInv .+ eCostPowTraInv .+ eCostH2GenInv .+ eCostH2StoInv .+ eCostH2TraInv 
    @objective(MP, Min, MP_obj .+ sum(theta[w] for w in W))

    # ---- MP Constraints ---- #
    @constraint(MP, cMaxPowGenRetCapTher[g in G], vRetPowGenCap[g]*pow_gen[g, :rep_capacity] <= pow_gen[g, :existing_cap])
    for g in G
      if 0 <= pow_gen[g, :max_cap_mw]
        @constraint(MP, cMaxPowGenCap[g in G], eTotPowGenCap[g] <= pow_gen[g, :max_cap_mw])
      end
    end
    @constraint(MP, cPowGenTotCapMin[g in G], pow_gen[g, :min_cap] <= eTotPowGenCap[g])
    @constraint(MP, cMaxRetPowSto[s in S], vRetPowStoCap[s]*pow_gen[s, :rep_capacity] <= pow_gen[s, :existing_cap_mwh])
    for s in S
      if 0 <= pow_gen[s, :max_cap_mwh]
        @constraint(MP, cMaxPowStoCap[s in S],eTotPowStoCap[s]<= pow_gen[s, :max_cap_mwh])
      end
    end
    @constraint(MP, cPowStoTotCapMin[s in S], pow_gen[s, :min_cap_mwh] .- eTotPowStoCap[s] <= 0)
    for l in L
      if 0 <= pow_lines[l, :line_max_reinforcement_mw]
        @constraint(MP, [l in L], vNewPowTraCap[l] <= pow_lines[l, :line_max_reinforcement_mw])
      end
    end
    @constraint(MP, cMaxRetH2Cap[h in H], vRetH2GenCap[h]*hsc_gen[h, :rep_capacity] <= hsc_gen[h, :existing_cap_tonne_p_hr])
    for h in H
      if 0 <= hsc_gen[h, :max_cap_tonne_p_hr]
        @constraint(MP, [h in H], eTotH2GenCap[h]<= hsc_gen[h, :max_cap_tonne_p_hr])
      end
    end
    @constraint(MP, cMinH2GenCap[h in H], hsc_gen[h, :min_cap_tonne_p_hr] <= eTotH2GenCap[h])
    @constraint(MP, cMaxRetH2StoCap[s in Q], vRetH2StoCap[s]*hsc_gen[s, :rep_capacity] <= hsc_gen[s, :existing_cap_tonne])
    for s in Q
      if 0 <= hsc_gen[s, :max_cap_stor_tonne]
        @constraint(MP, [s in S], eTotH2StoCap[s]<= hsc_gen[s, :max_cap_stor_tonne])
      end
    end
    @constraint(MP, cMinH2StoCap[s in Q], hsc_gen[s, :min_cap_stor_tonne] <= eTotH2StoCap[s])
    @constraint(MP, cMaxRetH2StorCompCap[s in Q], vRetH2StoCompCap[s] <= hsc_gen[s, :existing_cap_comp_tonne_hr])
    @constraint(MP, cMaxH2StorCompcCap[s in Q], eTotH2StoCompCap[s]<= eTotH2StoCap[s])
    @constraint(MP, cMinH2StorCompcCap[s in Q], 0.01*eTotH2StoCap[s] <= eTotH2StoCompCap[s])
    @constraint(MP, cMaxH2PipeNum[i in I], eTotH2Pipe[i] <= hsc_pipelines[i, :max_num_pipes])
    @constraint(MP, cMaxRetH2PipeNum[i in I], vRetH2Pipe[i] <= hsc_pipelines[i, :existing_num_pipes])
    @constraint(MP, cMaxRetH2PipeCompCap[i in I], vRetH2PipeCompCap[i]<=hsc_pipelines[i, :existing_comp_cap_tonne_hr])
    #Land Use Constraint on each zone
    @constraint(MP, cLandUse[z in Z], ePowGenLandUse[z] + ePowStoLandUse[z] + eH2GenLandUse[z] + eH2StoLandUse[z] + eH2PipeLandUse[z] <= zones[z, :available_land])
    #Policy
    @constraint(MP, cZonalEmissionCap, sum(vMaxEmissionByWeek[w] for w in W) <= 0.05*sum((h2_demand[w][t,z]*33.3) +pow_demand[w][t,z] for t in T, w in W, z in Z))
   

    # Defining the Sub Problem
    SP_models = Vector{Model}(undef, length(W))
    SP = Model(Gurobi.Optimizer)
    
    for w in W
        SP_models[w] = Model(Gurobi.Optimizer)
        set_optimizer_attribute(SP_models[w], "Method", 2)      # use barrier method
        #set_optimizer_attribute(SP_models[w], "Crossover", 0)
        set_optimizer_attribute(SP_models[w], "OutputFlag", 0)
        set_optimizer_attribute(SP_models[w], "LogToConsole", 0)
        set_optimizer_attribute(SP_models[w],"MIPGap",1e-3)
    
        # ---- SP Variables ---- #
        # Power Generation DV #
        @variable(SP_models[w], vPowGen[g in G, t in T]>=0)
        @variable(SP_models[w], vPowGenFirst[g in G]>=0)
        @variable(SP_models[w], vPowGenOnline[g in G_ther, t in T]>=0)
        @variable(SP_models[w], vPowGenOnlineFirst[g in G_ther]>=0)
        @variable(SP_models[w], vPowGenStart[g in G_ther, t in T]>=0)
        @variable(SP_models[w], vPowGenShut[g in G_ther, t in T]>=0)
        # Power Storage DV 
        @variable(SP_models[w], vPowStoCha[s in S, t in T]>=0)
        @variable(SP_models[w], vPowStoDis[s in S, t in T]>=0)
        @variable(SP_models[w], vPowSOC[s in S, t in T]>=0)
        @variable(SP_models[w], vPowSOCFirst[s in S]>=0)
        # Power Transmission DV #
        @variable(SP_models[w], vPowFlow[l in L, t in T])
        # HSC Generation DV
        @variable(SP_models[w], vH2Gen[h in H, t in T]>=0)
        @variable(SP_models[w], vH2GenFirst[h in H]>=0)
        @variable(SP_models[w], vH2GenStart[h in H_ther, t in T]>=0)
        @variable(SP_models[w], vH2GenShut[h in H_ther, t in T]>=0)
        @variable(SP_models[w], vH2GenOnline[h in H_ther, t in T]>=0)
        @variable(SP_models[w], vH2GenOnlineFirst[h in H_ther]>=0)
        # HSC Storage DV
        @variable(SP_models[w], vH2StoCha[s in Q, t in T]>=0)
        @variable(SP_models[w], vH2StoDis[s in Q, t in T]>=0)
        @variable(SP_models[w], vH2StoSOC[s in Q, t in T]>=0)
        @variable(SP_models[w], vH2StoSOCFirst[s in Q]>=0)
        # HSC Transmission DV
        @variable(SP_models[w], vH2Flow[i in I, t in T])
        # Policy Variables
        @variable(SP_models[w], vH2NSD[z in Z, t in T]>=0)
        @variable(SP_models[w], vPowNSD[z in Z, t in T]>=0)
        @variable(SP_models[w], vPowCrt[z in Z, t in T]>=0)
        @variable(SP_models[w], vH2Crt[z in Z, t in T]>=0)
        @variable(SP_models[w], vExtraEmissionByZone>=0)
        # Variables foe availability of generation, storage and transmission
        @variable(SP_models[w], eAvailPowGenCap[g in G]>=0)
        @variable(SP_models[w], eAvailPowGenUnit[g in G_ther]>=0)
        @variable(SP_models[w], eAvailPowStoCap[s in S]>=0)
        @variable(SP_models[w], eAvailPowTraCap[l in L]>=0)
        @variable(SP_models[w], eAvailH2GenCap[h in H]>=0)
        @variable(SP_models[w], eAvailH2GenUnit[h in H_ther]>=0)
        @variable(SP_models[w], eAvailH2StoCap[s in Q]>=0)
        @variable(SP_models[w], eAvailH2Pipe[i in I]>=0)
        @variable(SP_models[w], eAvailH2PipeCompCap[i in I]>=0)
        @variable(SP_models[w], eAvailH2StoCompCap[s in Q]>=0)
        @variable(SP_models[w], eMaxEmissionByWeek>=0)
        
        # ---- SP_models[w] Expressions ---- #
        # Power Generation Expressions #
        @expression(SP_models[w], ePowGenByZone[z in Z, t in T], sum(vPowGen[g,t] * (pow_gen[g, :zone] == z ? 1 : 0) for g in G))
        # Power Storage Expressions #
        @expression(SP_models[w], ePowStoChaByZone[z in Z, t in T], sum(vPowStoCha[s,t] * (pow_gen[s, :zone] == z ? 1 : 0) for s in S))
        @expression(SP_models[w], ePowStoDisByZone[z in Z, t in T], sum(vPowStoDis[s,t] * (pow_gen[s, :zone] == z ? 1 : 0) for s in S))
        # Power Transmission Expressions #
        @expression(SP_models[w], eNet_Pow_Flow[z in Z, t in T], sum(Pow_Network[l,z] * vPowFlow[l,t] for l in L))
        @expression(SP_models[w], ePow_Loss_By_Zone[z in Z, t in T], sum(abs(Pow_Network[l,z]) * (1/2) *vPowFlow[l,t] * pow_lines[l, :line_loss_percentage] for l in L))
        # HSC Generation Expressions
        @expression(SP_models[w], eH2GenEvap[h in H, t in T], hsc_gen[h, :boil_off]*vH2Gen[h,t])
        @expression(SP_models[w], eH2GenByZone[z in Z, t in T], sum((vH2Gen[h,t] - eH2GenEvap[h,t]) * (hsc_gen[h,:zone] == z ? 1 : 0) for h in H))
        @expression(SP_models[w], eH2StoChaByZone[z in Z, t in T], sum(vH2StoCha[s,t] * (hsc_gen[s, :zone] == z ? 1 : 0) for s in Q))
        @expression(SP_models[w], eH2StoDisByZone[z in Z, t in T], sum(vH2StoDis[s,t] * (hsc_gen[s, :zone] == z ? 1 : 0) for s in Q))
        @expression(SP_models[w], eH2FlowByZone[i in I, z in Z, t in T], H2_Network[i,z] * vH2Flow[i,t])  # H2 flow contribution to zone z (H2_Network is pipeline incidence matrix)
        @expression(SP_models[w], eNet_H2_Flow[z in Z, t in T], sum(H2_Network[i,z] * vH2Flow[i,t] for i in I))
        @expression(SP_models[w], eH2_Loss_By_Zone[z in Z, t in T], sum(abs(H2_Network[i,z]) * 0.5 * vH2Flow[i,t] * hsc_pipelines[i,:pipe_loss_coeff] for i in I))
        # Sector-coupling Demands
        @expression(SP_models[w], ePowDemandHSC[z in Z, t in T], 
            sum(hsc_gen[h, :pow_demand_mwh_p_tonne]*vH2Gen[h,t]*(hsc_gen[h, :zone]==z ? 1 : 0) for h in H) +
            sum(hsc_gen[s, :h2charge_mwh_p_tonne]*vH2StoCha[s,t]*(hsc_gen[s, :zone]==z ? 1 : 0) for s in Q) +
            sum(hsc_pipelines[i, :comp_pow_mwh_per_tonne]*vH2Flow[i,t]*(H2_Network[i,z]==1 ? 1 : 0) for i in I)
        )
        @expression(SP_models[w], eH2DemandPow[z in Z, t in T], sum(pow_gen[g, :h2_demand_tonne_p_mwh]*vPowGen[g,t]*(pow_gen[g, :zone]==z ? 1 : 0) for g in G))
        @expression(SP_models[w], pow_D[t in T, z in Z], pow_demand[w][t,z] .+ ePowDemandHSC[z,t])
        @expression(SP_models[w], H2_D[t in T, z in Z], h2_demand[w][t,z] .+ eH2DemandPow[z,t])
        # Emission Expressions
        @expression(SP_models[w], eEmissionByWeek, 
            sum(vPowGen[g,t]*pow_gen[g, :heat_rate_mmbtu_per_yr]*CO2_content[pow_gen[g, :fuel]] for g in G, t in T) + 
            sum(vH2Gen[h,t]*hsc_gen[h, :heat_rate_mmbtu_p_tonne]*CO2_content[hsc_gen[h, :fuel]] for h in H, t in T)
        )
        
        # Power Cost Expressions
        @expression(SP_models[w], eCostPowGenVar, sum((pow_gen[g, :vom_cost_mwh] + pow_gen[g, :heat_rate_mmbtu_per_yr] .* fuel_costs[pow_gen[g, :fuel]][t]) .* vPowGen[g,t] for g in G, t in T))
        @expression(SP_models[w], eCostPowGenStart, sum(pow_gen[g, :start_cost_per_mw] .* pow_gen[g, :rep_capacity] .* vPowGenStart[g,t] for g in G_ther, t in T)) 
        @expression(SP_models[w], eCostPowStoVar, sum(vPowStoCha[s,t] .* pow_gen[s, :vom_cost_mwh_charge] for s in S, t in T))
        # HSC Cost Expressions
        @expression(SP_models[w], eCostH2GenVar, sum((hsc_gen[h, :vom_cost_p_tonne] + hsc_gen[h, :heat_rate_mmbtu_p_tonne] .* fuel_costs[lowercase(hsc_gen[h, :fuel])][t]) .* vH2Gen[h,t] for h in H, t in T)    )
        @expression(SP_models[w], eCostH2StoVar, sum(hsc_gen[s, :var_om_cost_charge_p_tonne]*vH2StoCha[s,t] for s in Q, t in T))
        @expression(SP_models[w], eCostH2GenStart, sum(hsc_gen[h, :startup_cost_p_tonne_hr] .* hsc_gen[h, :rep_capacity] .* vH2GenStart[h,t] for h in H_ther, t in T))
        #NSD and Curtailment Cost Expressions 
        @expression(SP_models[w], eCostH2NSD, sum(vH2NSD[z,t] .* zones[z, :voll_hsc] for z in Z, t in T))
        @expression(SP_models[w], eCostPowNSD, sum(vPowNSD[z,t] .* zones[z, :voll_pow] for z in Z, t in T)) 
        @expression(SP_models[w], eCostPowCrt, sum(vPowCrt[z,t] .* zones[z, :pow_curtail_cost] for z in Z, t in T))
        @expression(SP_models[w], eCostH2Crt, sum(vH2Crt[z,t] .* zones[z, :h2_curtail_cost] for z in Z, t in T))
        #Emission Cost Expression
        @expression(SP_models[w], eEmissionCost, vExtraEmissionByZone * zones[1, :emission_cost] )
        
        # ---- SP_models[w] Objective ---- #  
        
        @objective(SP_models[w], Min, eCostPowGenVar .+ eCostPowStoVar .+ eCostPowGenStart .+ eCostH2GenVar .+ eCostH2StoVar .+ eCostH2GenStart .+ eCostPowNSD .+ eCostH2NSD .+ eCostPowCrt .+ eCostH2Crt .+ eEmissionCost)

        # ---- SP Constraints ---- #
        # Power Balance #
        @constraint(SP_models[w], cPowerBalance[z in Z, t in T],
                ePowGenByZone[z,t] .- eNet_Pow_Flow[z,t] .- ePow_Loss_By_Zone[z,t] .+ ePowStoDisByZone[z,t] .- ePowStoChaByZone[z,t] .+ vPowNSD[z,t] .- vPowCrt[z,t] == pow_D[t,z]
        )
        # H2 Balance #
        @constraint(SP_models[w], cH2Balance[z in Z, t in T],
            eH2GenByZone[z,t] .- eNet_H2_Flow[z,t ] .- eH2_Loss_By_Zone[z,t] .+ eH2StoDisByZone[z,t] .- eH2StoChaByZone[z,t] .+ vH2NSD[z,t] .- vH2Crt[z,t] == H2_D[t,z]
        )

        # Power Generation #
        @constraint(SP_models[w], cMaxPowGen[g in G_ren, t in T], vPowGen[g,t] - eAvailPowGenCap[g] * pow_gen_var[H_w[w][t], pow_gen[g, :resource]] <= 0)
        @constraint(SP_models[w], cPowOnlineUnits[g in G_ther, t in T], vPowGenOnline[g,t] .- eAvailPowGenUnit[g]<= 0)
        @constraint(SP_models[w], cPowStartLimits[g in G_ther,t in T], vPowGenStart[g,t] .- eAvailPowGenUnit[g] .+ vPowGenOnline[g,t]<= 0)
        @constraint(SP_models[w], cMaxPowGenTher[g in G_ther, t in T], vPowGen[g,t] .- (pow_gen[g, :rep_capacity].*pow_gen[g,:max_op_level].*vPowGenOnline[g,t]) <= 0)
        @constraint(SP_models[w], cMinPowGenTher[g in G_ther, t in T], (pow_gen[g, :rep_capacity]*pow_gen[g,:min_op_level]*vPowGenOnline[g,t]) .- vPowGen[g,t] <= 0)
        @constraint(SP_models[w], cPowShutLimits[g in G_ther, t in T], vPowGenShut[g,t] .- vPowGenOnline[g,t] <= 0)
        #Cyclic constraint on Thermal power units
        @constraint(SP_models[w], cPowUnitOnlineCon[g in G_ther, t in 2:length(T)], vPowGenOnline[g,t] - vPowGenOnline[g,t-1] == vPowGenStart[g,t]-vPowGenShut[g,t])
        @constraint(SP_models[w], cPowUnitOnlineFirst[g in G_ther], vPowGenOnline[g,1] - vPowGenOnlineFirst[g] == vPowGenStart[g,1]-vPowGenShut[g,1])
        @constraint(SP_models[w], cPowUnitOnlineCycle[g in G_ther], vPowGenOnlineFirst[g] == vPowGenOnline[g,168])
        #Ramping of non-thermal units 
        @constraint(SP_models[w], cPowGenRampUp[g in G_ren, t in 2:length(T)], vPowGen[g,t]-vPowGen[g,t-1] .- pow_gen[g, :ramp_up] * eAvailPowGenCap[g]<= 0) 
        @constraint(SP_models[w], cPowGenRampUpFirst[g in G_ren], vPowGen[g,1]-vPowGenFirst[g] .- pow_gen[g, :ramp_up]*eAvailPowGenCap[g]<=0)
        @constraint(SP_models[w], cPowGenRampDn[g in G_ren, t in 2:length(T)], vPowGen[g,t-1]-vPowGen[g,t] .- pow_gen[g, :ramp_dn] * eAvailPowGenCap[g]<= 0)
        @constraint(SP_models[w], cPowGenRampDnFirst[g in G_ren], vPowGenFirst[g]-vPowGen[g,1] .- pow_gen[g, :ramp_dn] * eAvailPowGenCap[g]<= 0)
        #Ramping of Thermal units
        @constraint(SP_models[w], cTherPowGenRampDn[g in G_ther,t in 2:length(T)], 
            vPowGen[g,t-1] .- vPowGen[g,t] .- pow_gen[g, :rep_capacity]*pow_gen[g,:ramp_dn] *(vPowGenOnline[g,t].-vPowGenStart[g,t])
            .+ pow_gen[g, :rep_capacity]*vPowGenStart[g,t]*pow_gen[g,:min_op_level] .- vPowGenShut[g,t]*pow_gen[g, :rep_capacity]*min(pow_gen_var[H_w[w][t], pow_gen[g, :resource]],max(pow_gen[g,:min_op_level],pow_gen[g,:ramp_dn]))<=0
        )
        @constraint(SP_models[w], cTherPowGenRampDnFirst[g in G_ther],    
        vPowGenFirst[g] .- vPowGen[g,1] .- pow_gen[g, :rep_capacity]*pow_gen[g,:ramp_dn] *(vPowGenOnline[g,1].-vPowGenStart[g,1])
        .+ pow_gen[g, :rep_capacity]*vPowGenStart[g,1]*pow_gen[g,:min_op_level] .- vPowGenShut[g,1]*pow_gen[g, :rep_capacity]*min(pow_gen_var[H_w[w][w], pow_gen[g, :resource]],max(pow_gen[g,:min_op_level],pow_gen[g,:ramp_dn]))<=0
        )

        @constraint(SP_models[w], cTherPowGenRampUp[g in G_ther, t in 2:length(T)], 
        vPowGen[g,t] .- vPowGen[g,t-1] .- pow_gen[g, :rep_capacity]*pow_gen[g,:ramp_up] *(vPowGenOnline[g,t].-vPowGenStart[g,t])
        .+ pow_gen[g, :rep_capacity]*vPowGenShut[g,t]*pow_gen[g,:min_op_level] .- vPowGenStart[g,t]*pow_gen[g, :rep_capacity]*min(pow_gen_var[H_w[w][t], pow_gen[g, :resource]],max(pow_gen[g,:min_op_level],pow_gen[g,:ramp_up]))<=0
        )   

        @constraint(SP_models[w], cTherPowGenRampUpFirst[g in G_ther], 
        vPowGen[g,1] .- vPowGenFirst[g] .- pow_gen[g, :rep_capacity]*pow_gen[g,:ramp_up] *(vPowGenOnline[g,1].-vPowGenStart[g,1])
        .- pow_gen[g, :rep_capacity]*vPowGenShut[g,1]*pow_gen[g,:min_op_level] .+ vPowGenStart[g,1]*pow_gen[g, :rep_capacity]*min(pow_gen_var[H_w[w][w], pow_gen[g, :resource]],max(pow_gen[g,:min_op_level],pow_gen[g,:ramp_up]))<=0
        )
        #Minimum up/donw time for thermal generators
        @constraint(SP_models[w], cMinUpTimePowGen[g in G_ther,t in T], sum(vPowGenStart[g,tt] for tt in intersect(T, (t - pow_gen[g, :up_time]):t)) .- vPowGenOnline[g,t]<= 0)
        @constraint(SP_models[w], cMinDnTimePowGen[g in G_ther,t in T], sum(vPowGenShut[g,tt] for tt in intersect(T, (t - pow_gen[g, :down_time]):t)) .- eAvailPowGenUnit[g] .+vPowGenOnline[g,t]<= 0)
        #Power storage constraints
        @constraint(SP_models[w], cPowStoBalance[s in S, t in 2:length(T)], vPowSOC[s,t] == (1-pow_gen[s,:etta_self_dis])*vPowSOC[s,t-1] + pow_gen[s,:etta_cha]*vPowStoCha[s,t] - (1/pow_gen[s,:etta_dis])*vPowStoDis[s,t])
        @constraint(SP_models[w], cPowStoBalanceFirst[s in S], vPowSOC[s,1] == (1-pow_gen[s,:etta_self_dis])*vPowSOCFirst[s] + pow_gen[s,:etta_cha]*vPowStoCha[s,1] .- (1/pow_gen[s,:etta_dis])*vPowStoDis[s,1])
        @constraint(SP_models[w], cPowStoMaxDis[s in S, t in 2:length(T)], vPowStoDis[s,t] .- pow_gen[s,:etta_dis]*vPowSOC[s,t-1]<= 0)
        @constraint(SP_models[w], cPowStoMaxDisFirst[s in S], vPowStoDis[s,1] .- pow_gen[s,:etta_dis]*vPowSOCFirst[s]<= 0)
        @constraint(SP_models[w], cPowSOCCycle[s in S], vPowSOCFirst[s] == vPowSOC[s,168])
        @constraint(SP_models[w], cPowStoMaxCha[s in S, t in T], vPowStoCha[s,t] .- eAvailPowStoCap[s]<=0)
        @constraint(SP_models[w], cPowStoSOCMax[s in S, t in T], vPowSOC[s,t] .- eAvailPowStoCap[s]*pow_gen[s, :max_op_level]<=0)
        @constraint(SP_models[w], cPowStoSOCMin[s in S, t in T], eAvailPowStoCap[s]*pow_gen[s, :min_op_level] .- vPowSOC[s,t]<=0)
        # Power Transmission #
        @constraints(SP_models[w], begin
        cMaxPowFlowOut[l in L, t in T],  vPowFlow[l,t] .- eAvailPowTraCap[l]<=0 
        cMaxPowFlowIn[l in L, t in T],  -eAvailPowTraCap[l] .- vPowFlow[l,t] <=0 
        end)
        
        @constraint(SP_models[w], cMaxH2GenVar[h in H_dis, t in T], vH2Gen[h,t] .- eAvailH2GenCap[h]*hsc_gen[h, :max_op_level]*hsc_gen_var[H_w[w][t], hsc_gen[h, :resource]] <= 0)
        @constraint(SP_models[w], cMinH2GenVar[h in H_dis, t in T], eAvailH2GenCap[h]*hsc_gen[h,:min_op_level]*hsc_gen_var[H_w[w][t], hsc_gen[h, :resource]]  .- vH2Gen[h,t] <= 0)
        # H2 Thermal units constraints
        @constraint(SP_models[w], cH2OnlineUnits[h in H_ther, t in T], vH2GenOnline[h,t] .- eAvailH2GenUnit[h]<= 0)
        @constraint(SP_models[w], cH2StartLimits[h in H_ther, t in T], vH2GenStart[h,t] .- eAvailH2GenUnit[h] .+ vH2GenOnline[h,t]<= 0 ) 
        @constraint(SP_models[w], cMaxH2GenTher[h in H_ther, t in T], vH2Gen[h,t] .- hsc_gen[h, :max_op_level] * hsc_gen[h, :rep_capacity] * vH2GenOnline[h,t]*hsc_gen_var[H_w[w][t], hsc_gen[h, :resource]]<= 0)
        @constraint(SP_models[w], cMinH2GenTher[h in H_ther, t in T], hsc_gen[h, :min_op_level] * hsc_gen[h, :rep_capacity] * vH2GenOnline[h,t]*hsc_gen_var[H_w[w][t], hsc_gen[h, :resource]] .- vH2Gen[h,t]<= 0)
        @constraint(SP_models[w], cH2ShutLimits[h in H_ther, t in T], vH2GenShut[h,t] .- vH2GenOnline[h,t]<= 0)
        @constraint(SP_models[w], cPowGenCycle[g in G_ren], vPowGenFirst[g] == vPowGen[g,168])
        # H2 thermal units cyclic constraints
        @constraint(SP_models[w], cH2UnitOnlineCon[h in H_ther, t in 2:length(T)], vH2GenOnline[h,t] .- vH2GenOnline[h,t-1] == vH2GenStart[h,t].- vH2GenShut[h,t])
        @constraint(SP_models[w], cH2UnitOnlineConFirst[h in H_ther], vH2GenOnline[h,1] .- vH2GenOnlineFirst[h] == vH2GenStart[h,1] .- vH2GenShut[h,1])
        @constraint(SP_models[w], cH2GenOnlineConCycle[h in H_ther], vH2GenOnlineFirst[h] == vH2GenOnline[h,168])
        # Min Up and Down time for Thermal H2 generators
        @constraint(SP_models[w], cMinUpTimeH2Gen[h in H_ther, t in T], sum(vH2GenStart[h,tt] for tt in intersect(T, (t - hsc_gen[h, :up_time]):t)) .- vH2GenOnline[h,t]<= 0)
        @constraint(SP_models[w], cMinDnTimeH2Gen[h in H_ther, t in T], sum(vH2GenShut[h,tt] for tt in intersect(T, (t - hsc_gen[h, :down_time]):t)) .- eAvailH2GenUnit[h] .+ vH2GenOnline[h,t]<= 0)

        @constraint(SP_models[w], cH2GenCycle[g in H_dis], vH2GenFirst[g] == vH2Gen[g,168])
        # Ramp constraints for diSP_models[w]atachable units
        @constraint(SP_models[w], cH2GenRampUp[g in H_dis, t in 2:length(T)], vH2Gen[g,t]-vH2Gen[g,t-1] .- hsc_gen[g, :ramp_up_percentage] * eAvailH2GenCap[g]<= 0) 
        @constraint(SP_models[w], cH2GenRampUpFirst[g in H_dis], vH2Gen[g,1]-vH2GenFirst[g] .- hsc_gen[g, :ramp_up_percentage]*eAvailH2GenCap[g]<=0)
        @constraint(SP_models[w], cH2GenRampDn[g in H_dis, t in 2:length(T)], vH2Gen[g,t-1]-vH2Gen[g,t] .- hsc_gen[g, :ramp_down_percentage] * eAvailH2GenCap[g]<= 0)
        @constraint(SP_models[w], cH2GenRampDnFirst[g in H_dis], vH2GenFirst[g]-vH2Gen[g,1] .- hsc_gen[g, :ramp_down_percentage] * eAvailH2GenCap[g]<= 0)
        #Ramp Constraint for Thermal units
        @constraint(SP_models[w], cTherH2GenRampDn[h in H_ther, t in 2:length(T)], 
        vH2Gen[h,t-1] .- vH2Gen[h,t] .- hsc_gen[h, :rep_capacity]*hsc_gen[h,:ramp_down_percentage] *(vH2GenOnline[h,t].-vH2GenStart[h,t])
        .+ hsc_gen[h, :rep_capacity]*vH2GenStart[h,t]*hsc_gen[h,:min_op_level] .- vH2GenShut[h,t]*hsc_gen[h, :rep_capacity]*min(hsc_gen_var[H_w[w][t], hsc_gen[h, :resource]],max(hsc_gen[h,:min_op_level],hsc_gen[h,:ramp_down_percentage]))<=0
        )
        @constraint(SP_models[w], cTherH2GenRampDnFirst[h in H_ther], 
            vH2GenFirst[h] .- vH2Gen[h,1] .- hsc_gen[h, :rep_capacity]*hsc_gen[h,:ramp_down_percentage] *(vH2GenOnline[h,1].-vH2GenStart[h,1])
            .+ hsc_gen[h, :rep_capacity]*vH2GenStart[h,1]*hsc_gen[h,:min_op_level] .- vH2GenShut[h,1]*hsc_gen[h, :rep_capacity]*min(hsc_gen_var[H_w[w][w], hsc_gen[h, :resource]],max(hsc_gen[h,:min_op_level],hsc_gen[h,:ramp_down_percentage]))<=0
        )

        @constraint(SP_models[w], cTherH2GenRampUp[h in H_ther, t in 2:length(T)], 
            vH2Gen[h,t] .- vH2Gen[h,t-1] .- hsc_gen[h, :rep_capacity]*hsc_gen[h,:ramp_up_percentage] *(vH2GenOnline[h,t].-vH2GenStart[h,t])
            .+ hsc_gen[h, :rep_capacity]*vH2GenShut[h,t]*hsc_gen[h,:min_op_level] .- vH2GenStart[h,t]*hsc_gen[h, :rep_capacity]*min(hsc_gen_var[H_w[w][t], hsc_gen[h, :resource]],max(hsc_gen[h,:min_op_level],hsc_gen[h,:ramp_up_percentage]))<=0
        )

        @constraint(SP_models[w], cTherH2GenRampUpFirst[g in H_ther], 
            vH2Gen[g,1] .- vH2GenFirst[g] .- hsc_gen[g, :rep_capacity]*hsc_gen[g,:ramp_up_percentage] *(vH2GenOnline[g,1].-vH2GenStart[g,1])
            .- hsc_gen[g, :rep_capacity]*vH2GenShut[g,1]*hsc_gen[g,:min_op_level] .+ vH2GenStart[g,1]*hsc_gen[g, :rep_capacity]*min(hsc_gen_var[H_w[w][w], hsc_gen[g, :resource]],max(hsc_gen[g,:min_op_level],hsc_gen[g,:ramp_up_percentage]))<=0
        )
        #H2 Storage constraints
        @constraint(SP_models[w], cH2StoBalance[s in Q, t in 2:length(T)], vH2StoSOC[s,t] == (1-hsc_gen[s, :etta_self_dis])*vH2StoSOC[s,t-1] + vH2StoCha[s,t]*hsc_gen[s,:etta_cha] - (1/hsc_gen[s,:etta_dis])*vH2StoDis[s,t])
        @constraint(SP_models[w], cH2StoBalanceFirst[s in Q], vH2StoSOC[s,1] == (1-hsc_gen[s, :etta_self_dis])*vH2StoSOCFirst[s] + vH2StoCha[s,1]*hsc_gen[s,:etta_cha] - (1/hsc_gen[s,:etta_dis])*vH2StoDis[s,1])
        @constraint(SP_models[w], cMaxH2StoSOC[s in Q, t in T], vH2StoSOC[s,t] .- hsc_gen[s,:h2stor_max_level]*eAvailH2StoCap[s]<= 0)
        @constraint(SP_models[w], cMinH2StoSOC[s in Q, t in T], hsc_gen[s,:h2stor_min_level]*eAvailH2StoCap[s] .- vH2StoSOC[s,t]<= 0 )
        @constraint(SP_models[w], cMaxH2StoChar[s in Q,t in T], vH2StoCha[s,t] .- eAvailH2StoCompCap[s] <= 0)
        
        # H2 Transmission constraints
        @constraints(SP_models[w], begin 
                        cMaxH2PipeFlowOut[i in I, t in T], vH2Flow[i,t] .- eAvailH2Pipe[i]*hsc_pipelines[i, :max_op_level] <= 0
                        cMaxH2PipeFlowIn[i in I, t in T], - vH2Flow[i,t] .- eAvailH2Pipe[i]*hsc_pipelines[i, :max_op_level] <= 0
        end)    
        @constraints(SP_models[w], begin
                        cMaxH2PipeFlowOutComp[i in I,t in T], vH2Flow[i,t] .- eAvailH2PipeCompCap[i] <= 0
                        cMaxH2PipeFlowInComp[i in I, t in T], -vH2Flow[i,t] .- eAvailH2PipeCompCap[i] <= 0
        end) 
        # Policy constraints
        #@constraint(SP_models[w], cPowNSD[z in Z, t in T], vPowNSD[z,t] .- zones[z, :pow_nsd_share]*pow_D[w,t,z] <= 0 )
        #@constraint(SP_models[w], cH2NSD[z in Z, t in T], vH2NSD[z,t] .- zones[z, :hsc_nsd_share]*H2_D[w, t,z] <= 0)
        #Emission constraint
        @constraint(SP_models[w], cEmissionCapByWeek, eEmissionByWeek .- vExtraEmissionByZone .- eMaxEmissionByWeek<= 0)
    end

    local LB, UB
    LB = -Inf
    UB = Inf
    k = 1
    max_iter = 500
    tolerence = 1e-2

    
    # Solve initial Master Problem to get a starting investment plan
    optimize!(MP)
    @assert termination_status(MP) == MOI.OPTIMAL
    LB = objective_value(MP)  # initial lower bound
    println("Initial MP objective (LB) = ", round(LB, digits=2))
    
    coupling = Dict{Int, Dict{Symbol, Any}}()
    for w in W
        coupling[w] = Dict(
            :gencap   => ConstraintRef[],
            :genunit  => ConstraintRef[],
            :stocap   => ConstraintRef[],
            :tracap   => ConstraintRef[],
            :h2gen    => ConstraintRef[],
            :h2genunit=> ConstraintRef[],
            :h2stocap => ConstraintRef[],
            :h2pipe   => ConstraintRef[],
            :h2pipecomp => ConstraintRef[],
            :emission => ConstraintRef[]
        )
    end


    for k in 1:max_iter
        
        println("────────────────────────────────────────")
        println(" BENDERS ITERATION $k")
        println("────────────────────────────────────────")

        vNewPowGenCap_val        = value.(vNewPowGenCap)
        vRetPowGenCap_val        = value.(vRetPowGenCap)
        vNewPowStoCap_val        = value.(vNewPowStoCap)
        vRetPowStoCap_val        = value.(vRetPowStoCap)
        vNewPowTraCap_val        = value.(vNewPowTraCap)
        vNewH2GenCap_val         = value.(vNewH2GenCap)
        vRetH2GenCap_val         = value.(vRetH2GenCap)
        vNewH2StoCap_val         = value.(vNewH2StoCap)
        vRetH2StoCap_val         = value.(vRetH2StoCap)
        vNewH2Pipe_val           = value.(vNewH2Pipe)
        vRetH2Pipe_val           = value.(vRetH2Pipe)
        vNewH2PipeCompCap_val    = value.(vNewH2PipeCompCap)
        vRetH2PipeCompCap_val    = value.(vRetH2PipeCompCap)
        vNewH2StoCompCap_val     = value.(vNewH2StoCompCap)
        vRetH2StoCompCap_val     = value.(vRetH2StoCompCap)
        vMaxEmissionByWeek_val   = value.(vMaxEmissionByWeek)
      
        for w in W
            cc = coupling[w]

            # Fix available capacities for week w to the values from MP's current solution
            availPowGenCap   = SP_models[w][:eAvailPowGenCap] 
            availPowGenUnit = SP_models[w][:eAvailPowGenUnit]
            availPowStoCap   = SP_models[w][:eAvailPowStoCap]
            availPowTraCap   = SP_models[w][:eAvailPowTraCap]
            availH2GenCap    = SP_models[w][:eAvailH2GenCap]
            availH2GenUnit  = SP_models[w][:eAvailH2GenUnit]
            availH2StoCap    = SP_models[w][:eAvailH2StoCap]
            availH2Pipe      = SP_models[w][:eAvailH2Pipe]
            availH2PipeCompCap = SP_models[w][:eAvailH2PipeCompCap]
            availH2StoCompCap = SP_models[w][:eAvailH2StoCompCap]
            availeEmissionMaxByWeek = SP_models[w][:eMaxEmissionByWeek]
            
            cAvailPowGenCap = @constraint(SP_models[w], [g in G], availPowGenCap[g] == value(eTotPowGenCap[g]))
            cAvailPowGenUnit = @constraint(SP_models[w], [g in G_ther], availPowGenUnit[g] == pow_gen[g,:num_units] + (value(vNewPowGenCap[g]) - value(vRetPowGenCap[g])))
            cAvailPowStoCap = @constraint(SP_models[w], [s in S], availPowStoCap[s] == value(eTotPowStoCap[s]))
            cAvailPowTraCap = @constraint(SP_models[w], [l in L], availPowTraCap[l] == vNewPowTraCap_val[l] + pow_lines[l, :existing_transmission_cap_mw] )
            cAvailH2GenCap = @constraint(SP_models[w], [h in H_dis], availH2GenCap[h] == value(eTotH2GenCap[h]))
            cAvailH2GenUnit = @constraint(SP_models[w], [h in H_ther], availH2GenUnit[h] == hsc_gen[h,:num_units] + (value(vNewH2GenCap[h]) - value(vRetH2GenCap[h])))
            cAvailH2StoCap = @constraint(SP_models[w], [s in Q], availH2StoCap[s] == value(eTotH2StoCap[s]))
            cAvailH2Pipe = @constraint(SP_models[w], [i in I], availH2Pipe[i] == value(eTotH2Pipe[i]))
            cAvailH2PipeCompCap = @constraint(SP_models[w], [i in I], availH2PipeCompCap[i] == value(eTotH2PipeCompCap[i]))
            cAvailH2StoCompCap = @constraint(SP_models[w], [s in Q], availH2StoCompCap[s] == value(eTotH2StoCompCap[s]))
            cEmission = @constraint(SP_models[w], availeEmissionMaxByWeek == value(vMaxEmissionByWeek[w]))

            cc[:gencap]    = Dict(g => cAvailPowGenCap[g]    for g in G)
            cc[:genunit]   = Dict(g => cAvailPowGenUnit[g]   for g in G_ther)
            cc[:stocap]    = Dict(s => cAvailPowStoCap[s]    for s in S)
            cc[:tracap]    = Dict(l => cAvailPowTraCap[l]    for l in L)
            cc[:h2gen]     = Dict(h => cAvailH2GenCap[h]     for h in H_dis)
            cc[:h2genunit] = Dict(h => cAvailH2GenUnit[h]   for h in H_ther)
            cc[:h2stocap]  = Dict(s => cAvailH2StoCap[s]    for s in Q)
            cc[:h2pipe]    = Dict(i => cAvailH2Pipe[i]      for i in I)
            cc[:h2pipecomp]= Dict(i => cAvailH2PipeCompCap[i] for i in I)
            cc[:h2stocomp] = Dict(s => cAvailH2StoCompCap[s] for s in Q)
            cc[:emission]  = cEmission

        end
        
        Threads.@threads for w in W
            optimize!(SP_models[w])
        end

        for w in W
            @assert termination_status(SP_models[w]) == MOI.OPTIMAL
        end
        
        total_sp_cost = sum(objective_value(SP_models[w]) for w in W) 
        invest_cost  = value(MP_obj)
        UB_candidate = invest_cost + total_sp_cost
        UB = min(UB, UB_candidate)
        println(" → Total SP cost = ", round(total_sp_cost,  digits=2))
        println(" → Candidate UB   = ", round(UB_candidate, digits=2), " (Best UB = ",    round(UB,           digits=2), ")")

        # Optimality cuts
        for w in W
            @constraint(MP,
                theta[w] >= objective_value(SP_models[w]) +
                sum(dual(coupling[w][:gencap][g])*pow_gen[g, :rep_capacity]*(vNewPowGenCap[g] - vRetPowGenCap[g] - vNewPowGenCap_val[g] + vRetPowGenCap_val[g]) for g in G_ren) +
                sum(dual(coupling[w][:genunit][g])*(vNewPowGenCap[g] - vRetPowGenCap[g] - vNewPowGenCap_val[g] + vRetPowGenCap_val[g]) for g in G_ther) +
                sum(dual(coupling[w][:stocap][s])*pow_gen[s, :rep_capacity]*(vNewPowStoCap[s] - vRetPowStoCap[s] - vNewPowStoCap_val[s] + vRetPowStoCap_val[s]) for s in S) +
                sum(dual(coupling[w][:tracap][l])*(vNewPowTraCap[l] - vNewPowTraCap_val[l]) for l in L) +
                sum(dual(coupling[w][:h2gen][h])*hsc_gen[h, :rep_capacity]*(vNewH2GenCap[h]-vRetH2GenCap[h] - vNewH2GenCap_val[h] + vRetH2GenCap_val[h]) for h in H_dis) +
                sum(dual(coupling[w][:h2genunit][h])*hsc_gen[h, :rep_capacity]*(vNewH2GenCap[h]-vRetH2GenCap[h] - vNewH2GenCap_val[h] + vRetH2GenCap_val[h]) for h in H_ther) +
                sum(dual(coupling[w][:h2stocap][s])*hsc_gen[s, :rep_capacity]*(vNewH2StoCap[s]-vRetH2StoCap[s] - vNewH2StoCap_val[s] + vRetH2StoCap_val[s]) for s in Q) +
                sum(dual(coupling[w][:h2pipe][i])*(vNewH2Pipe[i]-vRetH2Pipe[i] - vNewH2Pipe_val[i] + vRetH2Pipe_val[i]) for i in I) +
                sum(dual(coupling[w][:h2pipecomp][i])*(vNewH2PipeCompCap[i]-vRetH2PipeCompCap[i] - vNewH2PipeCompCap_val[i]) + vRetH2PipeCompCap_val[i] for i in I) +
                sum(dual(coupling[w][:h2stocomp][s])*(vNewH2StoCompCap[s]-vRetH2StoCompCap[s] - vNewH2StoCompCap_val[s] + vRetH2StoCompCap_val[s]) for s in Q) +
                dual(coupling[w][:emission])*(vMaxEmissionByWeek[w] - vMaxEmissionByWeek_val[w])
            )
        end
        
        optimize!(MP)
        println("MP status after adding cuts: ", termination_status(MP))
        @assert termination_status(MP) == MOI.OPTIMAL
        LB = objective_value(MP)
        println(" → Updated MP objective (LB) = ", round(LB, digits=2))
        println(" → Gap = ", round((UB - LB)*100/abs(LB), digits=2),"%")

        push!(iters, k)
        push!(LB_hist, LB  / 1e6)   
        push!(UB_hist, UB  / 1e6)
        #push!(gap_hist,(UB - LB) / 1e6)

        # update each curve *without* re-adding legend entries
        if k % 5 == 0 || (UB - LB)/abs(LB + eps()) <= tolerence
            # update curves without legends
            plot!(plt, iters, LB_hist, color=:blue, label=false)
            plot!(plt, iters, UB_hist, color=:red,  label=false)
            #plot!(plt, iters, gap_hist,color=:black,label="UB–LB")
    
            display(plt)  # redraw large plot
    
            println("Iter $k ▶  LB = $(round(LB/1e6, digits=3))e6, UB = $(round(UB/1e6, digits=3))e6, ",
                    "rel gap = $(round((UB-LB)/abs(LB + eps())*100, digits=3))%")
        end
        
        if (UB - LB)/abs(LB) <= tolerence
            println("Converged (gap = ", UB - LB, "). Optimal investment plan found.")
            break
        end
        
        for w in W
            m = SP_models[w]
            # collect *all* the ConstraintRefs you stored in coupling[w]
            all_crefs = ConstraintRef[]
            for entry in values(coupling[w])
                if entry === nothing
                    continue
                elseif entry isa ConstraintRef
                    push!(all_crefs, entry)
                elseif entry isa AbstractVector
                    append!(all_crefs, entry)
                elseif entry isa Dict
                    append!(all_crefs, collect(values(entry)))
                else
                    @warn "Unexpected coupling entry type: $(typeof(entry))"
                end
            end
        
            # delete them all from the model at once
            delete.(m, all_crefs)
        
            # now clear out your coupling[w] so next iteration starts fresh
            for k in keys(coupling[w])
                coupling[w][k] = nothing
            end
        end

        k += 1

    end
    
    println("Done: Iterations = $k, final gap = ", round(UB - LB, digits=4))
    if UB - LB <= tolerence
        println("Optimal objective = ", round(UB, digits=2))
    else
        println("Stopped without full convergence. Current best solution cost = ", round(UB, digits=2))
    end
    @assert termination_status(MP) == MOI.OPTIMAL

    println("Land‐use slack by zone:")
    for z in Z
        used  = value(
            ePowGenLandUse[z] +
            ePowStoLandUse[z] +
            eH2GenLandUse[z] +
            eH2StoLandUse[z] +
            eH2PipeLandUse[z]
        )
        avail = zones[z, :available_land]
        slack = avail - used

        println(" Zone $z slack = ",
                round(slack, digits=2),
                " (used=", round(used, digits=2),
                ", avail=", round(avail, digits=2), ")")
    end
    


#end

run_benders()