using JuMP, Gurobi, DataFrames, CSV


datadir = joinpath("/Users/rez/Documents/Engineering/Coding/Julia/PowHyCEM/Input_Data")


###Loading Data from CSV files###

#Loading power sector data
pow_gen = CSV.read(joinpath(datadir, "Powe_Gen_Data.csv"), DataFrame)
pow_gen_var = CSV.read(joinpath(datadir, "Pow_Gen_Var.csv"), DataFrame)
pow_load = CSV.read(joinpath(datadir, "Pow_Load.csv"), DataFrame)
pow_lines = CSV.read(joinpath(datadir, "Pow_Network.csv"), DataFrame)

#Loading HSC data
hsc_gen = CSV.read(joinpath(datadir, "HSC_Gen_Data.csv"), DataFrame)
hsc_gen_var = CSV.read(joinpath(datadir, "HSC_Gen_Var.csv"), DataFrame)
hsc_load = CSV.read(joinpath(datadir, "HSC_load.csv"), DataFrame)
hsc_pipelines = CSV.read(joinpath(datadir, "HSC_Pipelines.csv"), DataFrame)

#Loading Fuel Data 
fuel = CSV.read(joinpath(datadir, "Fuels_data.csv"), DataFrame)

#Loading Zone Data
zones = CSV.read(joinpath(datadir, "Zone_Char.csv"), DataFrame)

#rename all columns to lower Case
for f in [pow_gen, pow_load, pow_lines, hsc_gen, hsc_load, hsc_pipelines, fuel, zones]
    rename!(f,lowercase.(names(f)))
end

function value_to_df_2dim(var)
    solution = DataFrame(var.data, :auto)
    ax1 = var.axes[1]
    ax2 = var.axes[2]
    cols = names(solution)
    insertcols!(solution, 1, :r_id => ax1)
    solution = stack(solution, Not(:r_id), variable_name=:hour)
    solution.hour = foldl(replace, [cols[i] => ax2[i] for i in 1:length(ax2)], init=solution.hour)
    rename!(solution, :value => :gen)
    solution.hour = convert.(Int64,solution.hour)
    return solution
end


###Constructing the  DataFrame###
dfGen = pow_gen[pow_gen[!, :gen_type].>= 1, :] # set of all generators
dfStor = pow_gen[pow_gen[!, :stor_type].> 0, :] #set of all power storages
dfH2Gen = hsc_gen[hsc_gen[!, :h_gen_type].>= 1, :] #set of all H2 Generators
dfH2Stor = hsc_gen[hsc_gen[!, :h_stor].>= 1, :]

# Defining Sets
G = dfGen.r_id
S = pow_gen[pow_gen[!, :stor_type].> 0, :r_id] #set of all power storages
Z = zones.zones
#T = pow_load.time_index
T = pow_load[1:24, :time_index]
K = dfGen[dfGen[!, :gen_type] .==1, :r_id]
R = dfGen[dfGen[!, :gen_type] .>1, :r_id]
L = pow_lines.network_lines
H = dfH2Gen[dfH2Gen[!, :h_gen_type].>= 1, :r_id]
W = hsc_gen[hsc_gen[!, :h_stor].>= 1, :r_id]
I = hsc_pipelines.hsc_pipelines
J = hsc_gen[hsc_gen[!, :h_gen_type].== 1, :r_id]
E = hsc_gen[hsc_gen[!, :h_gen_type].==2, :r_id]

start_p = findall(s -> s == "load_mw_z1", names(pow_load))[1]
pow_demand = Matrix(pow_load[1:length(T), start_p:start_p + length(Z)-1])

start_h = findall(s -> s == "load_hsc_tonne_z1", names(hsc_load))[1]
h2_demand = Matrix(hsc_load[1:length(T), start_h:start_h + length(Z)-1])

fuels = names(fuel)[2:end]
costs = Matrix(fuel[2:end, 2:end])
CO2_content = fuel[1, 2:end] # tons CO2/MMBtu
fuel_costs = Dict{AbstractString, Array{Float64}}()
fuel_CO2 = Dict{AbstractString, Float64}()
fuel_type = dfGen[!,:fuel]
CO2_content
for i = 1:length(fuels)
    fuel_costs[fuels[i]] = costs[:,i] 
    # fuel_CO2 is kton/MMBTU with scaling, or ton/MMBTU without scaling.
    fuel_CO2[fuels[i]] = CO2_content[i] 
end

col_p = findall(s -> s == "z1", names(pow_lines))[1]
Pow_Network = pow_lines[1:length(L), col_p:col_p+length(Z)-1]

col_h = findall(s -> s == "z1", names(hsc_pipelines))[1]
H2_Network = hsc_pipelines[1:length(I), col_h:col_h+length(Z)-1]

#Defining the Model
CEM = Model(Gurobi.Optimizer)


####################################
#### Defining Decision Variables ###
####################################

  #################################
## Power Sector Decision Variables ##
  #################################

# Power Generation DV #
@variable(CEM, vNewPowGenCap[g in G]>=0)
@variable(CEM, vRetPowGenCap[g in G]>=0)
@variable(CEM, vPowGen[g in G, t in T]>=0)
@variable(CEM, vPowGenCommit[g in K, t in T], Bin)
@variable(CEM, vPowGenStart[g in K, t in T], Bin)
@variable(CEM, vPowGenShut[g in K, t in T], Bin)
@variable(CEM, vPowResUp[g in K, t in T]>=0)
@variable(CEM, vPowResDn[g in K, t in T]>=0)

# Power Storage DV #
@variable(CEM, vNewPowStoCap[s in S]>=0)
@variable(CEM, vRetPowStoCap[s in S]>=0)
@variable(CEM, vPowStoCha[s in S, t in T]>=0)
@variable(CEM, vPowStoDis[s in S, t in T]>=0)
@variable(CEM, vPowSOC[s in S, t in 0:length(T)]>=0)

# Power Transmission DV #
@variable(CEM, vNewPowTraCap[l in L]>=0)
@variable(CEM, vPowFlow[l in L, t in T])

# Non-Served Power demand
@variable(CEM, vPowNSD[z in Z, t in T]>=0)


  ########################
## HSC Decision Variables ##
  ########################

# HSC Generation DV
@variable(CEM, vNewH2GenCap[h in H]>=0)
@variable(CEM, vRetH2GenCap[h in H]>=0)
@variable(CEM, vH2Gen[h in H, t in T]>=0)
@variable(CEM, vH2GenStart[h in J, t in T], Bin)
@variable(CEM, vH2GenShut[h in J, t in T], Bin)
@variable(CEM, vH2GenCommit[h in J, t in T], Bin)

# HSC Storage DV
@variable(CEM, vNewH2StoCap[w in W]>=0)
@variable(CEM, vRetH2StoCap[w in W]>=0)
@variable(CEM, vNewH2StoCompCap[w in W]>=0)
@variable(CEM, vRetH2StoCompCap[w in W]>=0)
@variable(CEM, vH2StoCha[w in W, t in T]>=0)
@variable(CEM, vH2StoDis[w in W, t in T]>=0)
@variable(CEM, vH2StoSOC[w in W, t in 0:length(T)]>=0)

# HSC Transmission DV
@variable(CEM, vNewH2Pipe[i in I], Int)
@variable(CEM, vRetH2Pipe[i in I], Int)
@variable(CEM, vNewH2PipeCompCap[i in I]>=0)
@variable(CEM, vRetH2PipeCompCap[i in I]>=0)
@variable(CEM, vH2Flow[i in I, t in T])

# HSC NSD #
@variable(CEM, vH2NSD[z in Z, t in T]>=0)

###################
### Expressions ###
###################

  ########################
## Power Sector Expression ##
  ########################

# Power Generation Expressions #
@expression(CEM, eTotPowGenCap[g in G], dfGen[dfGen[!, :r_id] .==g, :existing_cap] .+ vNewPowGenCap[g] .- vRetPowGenCap[g])
@expression(CEM, ePowGenByZone[z in Z, t in T], sum(vPowGen[g, t] * (pow_gen[g, :zone] == z ? 1 : 0) for g in G))
@expression(CEM, eAuxPowGen_var[g in K, t in T], eTotPowGenCap[g]*vPowGenCommit[g,t])
@expression(CEM, eAuxPowGen[g in K, t in T], vPowGen[g,t] .-eAuxPowGen_var[g,t] * pow_gen[g,:min_op_level]) #Non-linear
@expression(CEM, eTotPowGenCapByZone[z in Z], sum(eTotPowGenCap[g]*(pow_gen[g,:zone] == z ? 1 : 0) for g in K))
@expression(CEM, ePowResReqUp[z in Z, t in T], 0.1 * eTotPowGenCapByZone[z] .+ 0.05 * pow_demand[t,z])
@expression(CEM, ePowResReqDn[z in Z, t in T], 0.05 * pow_demand[t,z])
@expression(CEM, ePowGenEmiByZone[z in Z], sum(CO2_content[pow_gen[g, :fuel]] * pow_gen[g,:heat_rate_mmbtu_per_yr]*vPowGen[g,t]*(pow_gen[g, :zone]==z ? 1 : 0) for g in K, t in T))
@expression(CEM, ePowGenLandUse[z in Z], sum((vNewPowGenCap[g] - vRetPowGenCap[g])*pow_gen[g, :land_use_km2_p_cap]*(pow_gen[g,:zone]==z ? 1 : 0) for g in G))

# Power Storage Expressions #
@expression(CEM, eTotPowStoCap[s in S], dfStor[dfStor.r_id .==s, :existing_cap] .+ vNewPowStoCap[s] .- vRetPowStoCap[s])
@expression(CEM, ePowStoChaByZone[z in Z, t in T], sum(vPowStoCha[s, t] * (pow_gen[s, :zone] == z ? 1 : 0) for s in S))
@expression(CEM, ePowStoDisByZone[z in Z, t in T], sum(vPowStoDis[s, t] * (pow_gen[s, :zone] == z ? 1 : 0) for s in S))
@expression(CEM, ePowStoLandUse[z in Z], sum((vNewPowStoCap[s] - vRetPowStoCap[s])*pow_gen[s, :land_use_km2_p_cap]*(pow_gen[s,:zone]==z ? 1 : 0) for s in S))

# Power Transmission Expressions #
@expression(CEM, eTotPowTraCap[l in L], pow_lines[l, :existing_transmission_cap_mw] .+ vNewPowTraCap[l]) #No retired cap considered for transmission lines + Not var cost for power flows
@expression(CEM, eNet_Export_Pow_Flow[z in Z,t in T], sum(Pow_Network[l,z] * vPowFlow[l,t] for l in L))
@expression(CEM, ePow_Loss_By_Zone[z in Z,t in T], sum(abs(Pow_Network[l,z]) * (1/2) *vPowFlow[l,t] * pow_lines[l, :line_loss_percentage] for l in L))

# Cost Expressions
@expression(CEM, eCostPowGenInv, sum(pow_gen[g, :inv_cost_per_mwyr] .* vNewPowGenCap[g] .+ dfGen[dfGen[!, :r_id] .==g, :fom_cost_per_mwyr] .* eTotPowGenCap[g] for g in G))
@expression(CEM, eCostPowGenVar, 
    sum((pow_gen[g, :vom_cost_mwh] + pow_gen[g, :heat_rate_mmbtu_per_yr] .* fuel_costs[pow_gen[g, :fuel]][t]) .* vPowGen[g,t] for g in G, t in T)    
)
@expression(CEM, eCostPowGenStart, sum(pow_gen[g, :start_cost_per_mw] .* eTotPowGenCap[g] .* vPowGenStart[g, t] for g in K, t in T))  #We have non-linearity here - Either change the equation of use piecewise-linear approximation - if you get a non-linear problem look here!
#For cost of Non-served demand we only consider $/MWh for each zone and will not consider demand segments
@expression(CEM, eCostPowNSD, sum(vPowNSD[z,t] .* zones[z, :voll_pow] for z in Z, t in T))
@expression(CEM, eCostPowStoInv, 
    sum(pow_gen[s, :inv_cost_per_mwhyr] .* vNewPowStoCap[s] .+ ((pow_gen[s, :fom_cost_per_mwhyr] + pow_gen[s, :fom_cost_charge_per_mwyr]) .* eTotPowStoCap[s]) for s in S)
)
@expression(CEM, eCostPowStoVar, sum((vPowStoCha[s,t] + vPowStoDis[s,t]) .* pow_gen[s, :vom_cost_mwh_charge] for s in S, t in T))
@expression(CEM, eCostPowTraInv, sum(pow_lines[l, :line_reinforcement_cost_per_mwyr] .* vNewPowTraCap[l] for l in L))


  ################
## HSC Expression ##
  ################

# HSC Generation Expressions
@expression(CEM, eTotH2GenCap[h in H], hsc_gen[h, :existing_cap_tonne_p_hr] .+ vNewH2GenCap[h] .- vRetH2GenCap[h])  
@expression(CEM, eH2GenEvap[h in H, t in T], hsc_gen[h, :boil_off]*vH2Gen[h,t])
@expression(CEM, eAuxH2Gen[h in J, t in T], vH2Gen[h,t] - eTotH2GenCap[h]*hsc_gen[h, :min_op_level]*vH2GenCommit[h,t])
@expression(CEM, eH2GenEmiByZone[z in Z], sum(CO2_content[hsc_gen[h, :fuel]] * hsc_gen[h,:heat_rate_mmbtu_p_tonne]*vH2Gen[h,t]*(1-hsc_gen[h,:ccs_rate])*(hsc_gen[h, :zone]==z ? 1 : 0) for h in J, t in T))
@expression(CEM, eH2GenLandUse[z in Z], sum((vNewH2GenCap[h]-vRetH2GenCap[h])*hsc_gen[h, :land_use_km2_p_cap]*(hsc_gen[h,:zone]==z ? 1 : 0) for h in H))

# HSC Storage Expression
@expression(CEM, eTotH2StoCap[w in W], hsc_gen[w, :existing_cap_tonne] + vNewH2StoCap[w] - vRetH2StoCap[w])
@expression(CEM, eTotH2StoCompCap[w in W], hsc_gen[w, :existing_cap_comp_tonne_hr] + vNewH2StoCompCap[w] - vRetH2StoCompCap[w])
@expression(CEM, eH2StoLandUse[z in Z], sum((vNewH2StoCap[w]-vRetH2StoCap[w])* hsc_gen[w, :land_use_km2_p_cap]*(hsc_gen[w,:zone]==z ? 1 : 0) for w in W))

# HSC Tramsmission Expression
@expression(CEM, eTotH2Pipe[i in I], hsc_pipelines[i, :existing_num_pipes] + vNewH2Pipe[i] - vRetH2Pipe[i])
@expression(CEM, eTotH2PipeCompCap[i in I], hsc_pipelines[i, :existing_comp_cap_tonne_hr] + vNewH2PipeCompCap[i] - vRetH2PipeCompCap[i])
@expression(CEM, eH2PipeLandUse[z in Z],0.5 * sum(hsc_pipelines[i, :land_use_km2_p_length]*hsc_pipelines[i, :distance]*(vNewH2Pipe[i]-vRetH2Pipe[i])*abs(H2_Network[i,z]) for i in I))

# HSC Cost Expressions
@expression(CEM, eCostH2GenInv, sum(hsc_gen[h, :inv_cost_tonne_hr_p_yr]*vNewH2GenCap[h] + hsc_gen[h, :fom_cost_p_tonne_p_hr_yr]*eTotH2GenCap[h] for h in H))

@expression(CEM, eCostH2StoInv, sum(hsc_gen[w, :inv_cost_tonne_p_yr]*vNewH2StoCap[w] + hsc_gen[w, :inv_cost_comp_tonne_hr_p_yr]*vNewH2StoCompCap[w]
                                  + hsc_gen[w, :fom_cost_p_tonne_p_yr]*eTotH2StoCap[w] + hsc_gen[w, :fom_cost_comp_tonne_hr_p_yr]*eTotH2StoCompCap[w] for w in W)
)
@expression(CEM, eCostH2StoCompInv, sum((hsc_gen[w, :inv_cost_comp_tonne_hr_p_yr]*vNewH2StoCompCap[w]) + (hsc_gen[w, :fom_cost_comp_tonne_hr_p_yr]*eTotH2StoCompCap[w]) for w in W))

@expression(CEM, eCostH2TraInv, sum(hsc_pipelines[l, :investment_cost_per_length]*hsc_pipelines[l, :distance]*vNewH2Pipe[l] +
                                    hsc_pipelines[l, :fom_per_length]*hsc_pipelines[l, :distance]*eTotH2Pipe[l] +
                                    hsc_pipelines[l, :compressor_inv_per_length]*hsc_pipelines[l, :distance]*vNewH2PipeCompCap[l] + 
                                    hsc_pipelines[l, :fom_comp_p_tonne_hr]*eTotH2PipeCompCap[l] for l in L)
)

@expression(CEM, eCostH2GenVar, 
    sum((hsc_gen[h, :vom_cost_p_tonne] + hsc_gen[h, :heat_rate_mmbtu_p_tonne] .* fuel_costs[lowercase(hsc_gen[h, :fuel])][t]) .* vH2Gen[h,t] for h in H, t in T)    
)

@expression(CEM, eCostH2GenStart, sum(hsc_gen[h, :startup_cost_p_tonne_hr] .* eTotH2GenCap[h] .* vH2GenStart[h, t] for h in J, t in T)) 

@expression(CEM, eCostH2NSD, sum(vH2NSD[z,t] .* zones[z, :voll_hsc] for z in Z, t in T))

@expression(CEM, ePowDemandHSC[t in T, z in Z], sum(hsc_gen[h, :pow_demand_mwh_p_tonne]*vH2Gen[h,t]*(hsc_gen[h, :zone]==z ? 1 : 0) for h in H) + 
                                                sum(hsc_gen[w, :h2charge_mwh_p_tonne]*vH2StoCha[w,t]*(hsc_gen[w, :zone]==z ? 1 : 0) for w in W) +
                                                sum(hsc_pipelines[i, :comp_pow_mwh_per_tonne]*vH2Flow[i,t]*(H2_Network[i,z]==1 ? 1 : 0) for i in I) 
)

@expression(CEM, eH2DemandPow[z in Z,t in T], sum(pow_gen[g, :h2_demand_tonne_p_mwh]*vPowGen[g,t]*(pow_gen[g, :zone]==z ? 1 : 0) for g in G))
@expression(CEM, pow_D[t in T, z in Z], pow_demand[t,z] .+ ePowDemandHSC[t,z])
@expression(CEM, H2_D[t in T, z in Z], h2_demand[t,z] .+ eH2DemandPow[z,t])

#Defining Objective Function

obj = (eCostPowGenInv .+ eCostPowStoInv .+ eCostPowTraInv .+ eCostPowGenVar .+ eCostPowStoVar .+ eCostPowNSD .+ eCostPowGenStart) .+ (eCostH2GenInv .+ eCostH2StoInv .+ eCostH2TraInv .+ eCostH2GenVar .+ eCostH2NSD .+ eCostH2StoCompInv .+ eCostH2GenStart)

@objective(CEM, Min, obj)


###################
### Constraints ###
###################

  ##########################
## Power Sector Constraints ##
  #########################

# Power Balance #
@constraint(CEM, cPowerBalance[z in Z, t in T],
    ePowGenByZone[z,t] .- eNet_Export_Pow_Flow[z,t] .- ePow_Loss_By_Zone[z,t] .+ ePowStoDisByZone[z,t] .- ePowStoChaByZone[z,t] .+ vPowNSD[z,t] == pow_D[t,z]
)

# Power Generation #
@constraint(CEM, cMaxPowGenRetCap[g in G], vRetPowGenCap[g] <= pow_gen[g, :existing_cap])

for g in G
  if 0 <= pow_gen[g, :max_cap_mw]
    @constraint(CEM, cMaxPowGenCap[g in G], eTotPowGenCap[g] <= pow_gen[g, :max_cap_mw])
  end
end

@constraint(CEM, cPowGenTotCapMin[g in G], pow_gen[g, :min_cap] <= eTotPowGenCap[g])
@constraint(CEM, cMaxPowGen[g in R, t in T], vPowGen[g,t] .- eTotPowGenCap[g] .* pow_gen_var[t, pow_gen[g, :resource]] <= 0)
@constraint(CEM, cMaxPowGenTher[g in K, t in T], vPowGen[g,t] .- (eTotPowGenCap[g]*pow_gen[g,:max_op_level]*vPowGenCommit[g,t]) <= 0)
@constraint(CEM, cMinPowGenTher[g in K, t in T], (eTotPowGenCap[g]*pow_gen[g,:min_op_level]*vPowGenCommit[g,t]) .- vPowGen[g,t] <= 0) ## Non-linear
@constraint(CEM, cMinUpTimePowGen[g in K, t in T], sum(vPowGenStart[g,tt] for tt in intersect(T, (t - pow_gen[g, :up_time]):t)) <= vPowGenCommit[g,t])
@constraint(CEM, cMinDnTimePowGen[g in K, t in T], sum(vPowGenShut[g,tt] for tt in intersect(T, (t - pow_gen[g, :down_time]):t)) <= 1 - vPowGenCommit[g,t])
@constraint(CEM, cPowGenCommitCon[g in K, t in 1:length(T)-1], vPowGenCommit[g,t+1]-vPowGenCommit[g,t] == vPowGenStart[g,t+1]-vPowGenShut[g,t+1])
#Ramp and auxilary power generation constraints
@constraint(CEM, [g in K, t in T], eAuxPowGen_var[g,t]<= eTotPowGenCap[g])
for g in K
  if pow_gen[g, :max_cap_mw] != -1
    @constraint(CEM, [g in K, t in T], eAuxPowGen_var[g,t]<= vPowGenCommit[g,t]*pow_gen[g, :max_cap_mw])
  else
    @constraint(CEM, [g in K, t in T], eAuxPowGen_var[g,t]<= 10*eTotPowGenCap[g])
  end
end

for g in K
  if pow_gen[g, :max_cap_mw] != -1
    @constraint(CEM, [g in K, t in T], eAuxPowGen_var[g,t]>= eTotPowGenCap[g] - pow_gen[g, :max_cap_mw]*(1-vPowGenCommit[g,t]))
  else
    @constraint(CEM, [g in K, t in T], eAuxPowGen_var[g,t] >= eTotPowGenCap[g] - 10*eTotPowGenCap[g]*(1-vPowGenCommit[g,t]))
  end
end

@constraint(CEM, cAuxPowGenRampUp[g in K,t in 1:length(T)-1], eAuxPowGen[g,t+1] .- eAuxPowGen[g,t] <= eTotPowGenCap[g]*pow_gen[g,:ramp_up]) #For Thermal units
@constraint(CEM, cAuxPowGenRampDn[g in K,t in 1:length(T)-1], eAuxPowGen[g,t] .- eAuxPowGen[g,t+1]  <= eTotPowGenCap[g]*pow_gen[g,:ramp_dn])
@constraint(CEM, cPowGenRampUp[g in R, t in 1:length(T)-1], vPowGen[g,t+1] .- vPowGen[g,t] .- eTotPowGenCap[g]*pow_gen[g,:ramp_up]<= 0)  # For Non-Thermal units
@constraint(CEM, cPowGenRampDn[g in R, t in 1:length(T)-1], vPowGen[g,t] .- vPowGen[g,t+1] .- eTotPowGenCap[g]*pow_gen[g,:ramp_dn]<= 0)
#Spinning Reserve Constraints
@constraint(CEM, cPowResUpMax[g in K, t in T], vPowResUp[g,t] .+ vPowGen[g,t] .-  eTotPowGenCap[g]*pow_gen[g,:max_op_level]*vPowGenCommit[g,t] <=0) #Non-linear
@constraint(CEM, cPowResDnMax[g in K, t in T], vPowResDn[g,t] .+ eTotPowGenCap[g]*pow_gen[g,:min_op_level]*vPowGenCommit[g,t] .- vPowGen[g,t] <= 0) #Non-linear
@constraint(CEM, cPowResUP[g in K, t in T], vPowResUp[g,t] .- eTotPowGenCap[g]*pow_gen[g,:ramp_up] <=0)
@constraint(CEM, cPowResDn[g in K, t in T], vPowResDn[g,t] .- eTotPowGenCap[g]*pow_gen[g,:ramp_dn] <=0)
@constraint(CEM, cPowResReqUp[z in Z, t in T], ePowResReqUp[z,t] .- sum(vPowResUp[g,t] * (pow_gen[g,:zone]==z ? 1 : 0) for g in K)<= 0)
@constraint(CEM, cPowResReqDn[z in Z, t in T], ePowResReqDn[z,t] .- sum(vPowResDn[g,t] * (pow_gen[g,:zone]==z ? 1 : 0) for g in K)<= 0)

# Power Storage #
@constraint(CEM, cMaxRetPowSto[s in S], vRetPowStoCap[s] <= pow_gen[s, :existing_cap_mwh])
for s in S
  if 0 <= pow_gen[s, :max_cap_mwh]
    @constraint(CEM, cMaxPowStoCap[s in S],eTotPowStoCap[s]<= pow_gen[s, :max_cap_mwh])
  end
end
@constraint(CEM, cPowStoTotCapMin[s in S], pow_gen[s, :min_cap_mwh] <= eTotPowStoCap[s])
@constraint(CEM, cPowStoBalanceInit[s in S, t=[0]], vPowSOC[s,t]==0)
@constraint(CEM, cPowStoBalance[s in S, t in 1:length(T)], vPowSOC[s,t] == (1-pow_gen[s,:etta_self_dis])*vPowSOC[s,t-1] + pow_gen[s,:etta_cha]*vPowStoCha[s,t] - (1/pow_gen[s,:etta_dis])*vPowStoDis[s,t])
@constraint(CEM, cPowStoMaxDis[s in S, t in T], vPowStoDis[s,t] <= pow_gen[s,:etta_dis]*vPowSOC[s,t-1])

# Power Transmission #
@constraints(CEM, begin
                cMaxPowFlowOut[l in L, t in T],  vPowFlow[l,t] <= eTotPowTraCap[l]
                cMaxPowFlowIn[l in L, t in T], vPowFlow[l,t] >= -eTotPowTraCap[l]
end)

for l in L
  if 0 <= pow_lines[l, :line_max_reinforcement_mw]
    @constraint(CEM, [l in L], eTotPowTraCap[l] <= pow_lines[l, :line_max_reinforcement_mw])
  end
end

# Power Non-served Demand #
@constraint(CEM, cPowNSD[z in Z, t in T], vPowNSD[z,t] <= zones[z, :pow_nsd_share]*pow_D[t,z] )


  #################
## HSC Constraints ##
  #################

# H2 Balance constraint
@constraint(CEM, cH2Balance[z in Z,t in T],
  sum((vH2Gen[h,t] .- eH2GenEvap[h,t])*(hsc_gen[h, :zone]==z ? 1 : 0) for h in H) .- sum(H2_Network[i,z]*vH2Flow[i,t] for i in I) .- 0.5*sum(H2_Network[i,z]*vH2Flow[i,t]*hsc_pipelines[i, :pipe_loss_coeff] for i in I) .+
  sum((vH2StoDis[w,t] - vH2StoCha[w,t])*(hsc_gen[w, :zone]==z ? 1 : 0) for w in W) + vH2NSD[z,t] == h2_demand[t,z]
)


# H2 Generation constraint
@constraint(CEM, cMaxRetH2Cap[h in H], vRetH2GenCap[h] <= hsc_gen[h, :existing_cap_tonne_p_hr])
for h in H
  if 0 <= hsc_gen[h, :max_cap_tonne_p_hr]
    @constraint(CEM, [h in H], eTotH2GenCap[h]<= hsc_gen[h, :max_cap_tonne_p_hr])
  end
end
@constraint(CEM, cMinH2GenCap[h in H], hsc_gen[h, :min_cap_tonne_p_hr] <= eTotH2GenCap[h])
@constraint(CEM, cMaxH2GenVar[h in E, t in T], vH2Gen[h,t] <= eTotH2GenCap[h]*hsc_gen[h, :max_op_level]*hsc_gen_var[t, hsc_gen[h,:resource]])
@constraint(CEM, cMinH2GenVar[h in E, t in T], eTotH2GenCap[h]*hsc_gen[h,:min_op_level]  <= vH2Gen[h,t])
@constraint(CEM, cMaxH2GenTher[h in J, t in T], vH2Gen[h,t] <= hsc_gen[h, :max_op_level]*eTotH2GenCap[h]*vH2GenCommit[h,t])
@constraint(CEM, cMinH2GenTher[h in J, t in T], hsc_gen[h, :min_op_level]*eTotH2GenCap[h]*vH2GenCommit[h,t] <= vH2Gen[h,t])
@constraint(CEM, cMinUpTimeH2Gen[h in J, t in T], sum(vH2GenStart[h,tt] for tt in intersect(T, (t - hsc_gen[h, :up_time]):t)) <= vH2GenCommit[h,t])
@constraint(CEM, cMinDnTimeH2Gen[h in J, t in T], sum(vH2GenShut[h,tt] for tt in intersect(T, (t - hsc_gen[h, :down_time]):t)) <= 1 - vH2GenCommit[h,t])
@constraint(CEM, cH2GenCommitCon[h in J, t in 1:length(T)-1], vH2GenCommit[h,t+1]-vH2GenCommit[h,t] == vH2GenStart[h,t+1]-vH2GenShut[h,t+1])
#Ramp and Auxilary H2 Gen Constraints
@constraint(CEM, cAuxH2GenRampUp[h in J, t in 1:length(T)-1], eAuxH2Gen[h,t+1]-eAuxH2Gen[h,t] <= hsc_gen[h, :ramp_up_percentage]*eTotH2GenCap[h])
@constraint(CEM, cAuxH2GenRampDn[h in J, t in 1:length(T)-1], eAuxH2Gen[h,t] - eAuxH2Gen[h,t+1] <= hsc_gen[h, :ramp_down_percentage]*eTotH2GenCap[h])
@constraint(CEM, cH2GenRampUp[h in J, t in 1:length(T)-1], vH2Gen[h,t+1] - vH2Gen[h,t] <= hsc_gen[h, :ramp_up_percentage]*eTotH2GenCap[h])
@constraint(CEM, cH2GenRampDn[h in J, t in 1:length(T)-1], vH2Gen[h,t] - vH2Gen[h,t+1] <= hsc_gen[h, :ramp_down_percentage]*eTotH2GenCap[h])

# H2 Storage constraint
@constraint(CEM, cMaxRetH2StoCap[w in W], vRetH2StoCap[w] <= hsc_gen[w, :existing_cap_tonne])
for w in W
  if 0 <= hsc_gen[w, :max_cap_stor_tonne]
    @constraint(CEM, [w in W], eTotH2StoCap[w]<= hsc_gen[w, :max_cap_stor_tonne])
  end
end
@constraint(CEM, cMinH2StoCap[w in W], hsc_gen[w, :min_cap_stor_tonne] <= eTotH2StoCap[w])
@constraint(CEM, cH2StoBalance[w in W, t in T], vH2StoSOC[w,t] == (1-hsc_gen[w, :etta_self_dis])*vH2StoSOC[w,t-1] + vH2StoCha[w,t]*hsc_gen[w,:etta_cha] - (1/hsc_gen[w,:etta_dis])*vH2StoDis[w,t])
@constraint(CEM, cH2StoSOCInit[w in W, t=[0]], vH2StoSOC[w,t]==0)
@constraint(CEM, cMaxH2StoSOC[w in W, t in T], vH2StoSOC[w,t]<= hsc_gen[w,:h2stor_max_level]*eTotH2StoCap[w])
@constraint(CEM, cMinH2StoSOC[w in W, t in T], hsc_gen[w,:h2stor_min_level]*eTotH2StoCap[w]<= vH2StoSOC[w,t])
@constraint(CEM, cMaxRetH2StorCompCap[w in W], vRetH2StoCompCap[w] <= hsc_gen[w, :existing_cap_comp_tonne_hr])
@constraint(CEM, cMaxH2StorCompcCap[w in W], eTotH2StoCompCap[w]<= eTotH2StoCap[w])
@constraint(CEM, cMinH2StorCompcCap[w in W], 0.1*eTotH2StoCap[w] <= eTotH2StoCompCap[w])
@constraint(CEM, cMaxH2StoChar[w in W,t in T], vH2StoCha[w,t] <= eTotH2StoCompCap[w])

# H2 Transmission constraints
@constraint(CEM, cMaxH2PipeNum[i in I], eTotH2Pipe[i] <= hsc_pipelines[i, :max_num_pipes])
@constraint(CEM, cMaxRetH2PipeNum[i in I], vRetH2Pipe[i] <= hsc_pipelines[i, :existing_num_pipes])
#@constraints(CEM, begin
#                  cMinH2PipeFlowOut[i in I, t in T], hsc_pipelines[i, :max_flow_tonne_per_h_per_pipe]*hsc_pipelines[i, :min_op_level] <= vH2Flow[i,t]
#                  cMinH2PipeFlowiIn[i in I, t in T], - hsc_pipelines[i, :max_flow_tonne_per_h_per_pipe]*hsc_pipelines[i, :min_op_level] >= vH2Flow[i,t]
#end)

@constraints(CEM, begin 
                  cMaxH2PipeFlowOut[i in I, t in T], vH2Flow[i,t] <= hsc_pipelines[i, :max_flow_tonne_per_h_per_pipe]*hsc_pipelines[i, :max_op_level]*eTotH2Pipe[i]
                  cMaxH2PipeFlowIn[i in I, t in T], vH2Flow[i,t] >= -hsc_pipelines[i, :max_flow_tonne_per_h_per_pipe]*hsc_pipelines[i, :max_op_level]*eTotH2Pipe[i]
end)

@constraint(CEM, cMaxRetH2PipeCompCap[i in I], vRetH2PipeCompCap[i]<=hsc_pipelines[i, :existing_comp_cap_tonne_hr])
@constraints(CEM, begin
                  cMaxH2PipeFlowOutComp[i in I,t in T], vH2Flow[i,t] <= eTotH2PipeCompCap[i]
                  cMaxH2PipeFlowInComp[i in I, t in T], vH2Flow[i,t] >= -eTotH2PipeCompCap[i]
end)

# H2 NSD Constraints
@constraint(CEM, cH2NSD[z in Z, t in T], vH2NSD[z,t] <= zones[z, :hsc_nsd_share]*H2_D[t,z])


# System Emission Constraint by zone
# H2 Emissions # Method 1  = Emission should not exceed specified emission cap
@constraint(CEM, cEmission[z in Z], ePowGenEmiByZone[z] + eH2GenEmiByZone[z] <= zones[z, :emission_cap_tonne_co2])

#method 2 = We can produced as much emission as we want, however we will be fined if we cross the cap
#=
if zones[z, :emission_cap_tonne_co2] <= ePowGenEmiByZone[z] + eH2GenEmiByZone[z]
  CEM(:eObj) += sum((ePowGenEmiByZone[z] + eH2GenEmiByZone[z] - zones[z, :emission_cap_tonne_co2]) * zones[z, :emission_cost] for z in Z)
end
=#

# Land Use Constraint on each zone
@constraint(CEM, cLandUse[z in Z], ePowGenLandUse[z] + ePowStoLandUse[z] + eH2GenLandUse[z] + eH2StoLandUse[z] + eH2PipeLandUse[z] <= zones[z, :available_land])

#set_optimizer_attribute(CEM, "DualReductions", 0)
try
  optimize!(CEM)
  if termination_status(CEM) == MOI.OPTIMAL
      println("Objective value: ", objective_value(CEM))
  elseif termination_status(CEM) == MOI.INFEASIBLE
    println("Model is infeasible. Computing IIS...")

    set_optimizer_attribute(CEM, "IISMethod", 1)  # Enable IIS computation
    compute_conflict!(CEM)  # Compute conflicting constraints
    println("Conflicting constraints:")
    for (c, status) in list_conflicting_constraints(CEM)
        println(c)
    end
  else
      println("Model status: ", termination_status(CEM))
  end
catch e
  println("Optimization failed: ", e)
end