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
T = pow_load[1:8760, :time_index]
G_ther = dfGen[dfGen[!, :gen_type] .==1, :r_id] #set of all thermal power units
G_ren = dfGen[dfGen[!, :gen_type] .>1, :r_id] # set of all renewable dispatchable power units
V = dfGen[dfGen[!, :gen_type].==3, :r_id]  #set of all fuel cells
L = pow_lines.network_lines

H = dfH2Gen[dfH2Gen[!, :h_gen_type].>= 1, :r_id] #set of all hydrogen generators
Q = hsc_gen[hsc_gen[!, :h_stor].>= 1, :r_id] # set of hydrogen storages
I = hsc_pipelines.hsc_pipelines
H_ther = hsc_gen[hsc_gen[!, :h_gen_type].== 1, :r_id]  #set of all thermal hydrogen generators
H_dis = hsc_gen[hsc_gen[!, :h_gen_type].==2, :r_id]   #set of all renewable hydrogen generators

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

  #--------------------------------#
## Power Sector Decision Variables ##
  #--------------------------------#

# Power Generation DV #
@variable(CEM, vNewPowGenCap[g in G], Int)
@variable(CEM, vRetPowGenCap[g in G], Int)
@variable(CEM, vPowGen[g in G, t in T]>=0)
@variable(CEM, vPowGenOnline[g in G_ther, t in T], Int)
@variable(CEM, vPowGenStart[g in G_ther, t in T], Int)
@variable(CEM, vPowGenShut[g in G_ther, t in T], Int)
@variable(CEM, vPowResUp[g in G_ther, t in T]>=0)
@variable(CEM, vPowResDn[g in G_ther, t in T]>=0)

# Power Storage DV #
@variable(CEM, vNewPowStoCap[s in S], Int)
@variable(CEM, vRetPowStoCap[s in S], Int)
@variable(CEM, vPowStoCha[s in S, t in T]>=0)
@variable(CEM, vPowStoDis[s in S, t in T]>=0)
@variable(CEM, vPowSOC[s in S, t in 0:length(T)]>=0)

# Power Transmission DV #
@variable(CEM, vNewPowTraCap[l in L], Int)
@variable(CEM, vPowFlow[l in L, t in T])

# Non-Served Power demand
@variable(CEM, vPowNSD[z in Z, t in T]>=0)


  ########################
## HSC Decision Variables ##
  ########################

# HSC Generation DV
@variable(CEM, vNewH2GenCap[h in H], Int)
@variable(CEM, vRetH2GenCap[h in H], Int)
@variable(CEM, vH2Gen[h in H, t in T]>=0)
@variable(CEM, vH2GenStart[h in H_ther, t in T], Int)
@variable(CEM, vH2GenShut[h in H_ther, t in T], Int)
@variable(CEM, vH2GenOnline[h in H_ther, t in T], Int)


# HSC Storage DV
@variable(CEM, vNewH2StoCap[s in Q]>=0, Int)
@variable(CEM, vRetH2StoCap[s in Q]>=0, Int)
@variable(CEM, vNewH2StoCompCap[s in Q]>=0, Int)
@variable(CEM, vRetH2StoCompCap[s in Q]>=0, Int)
@variable(CEM, vH2StoCha[s in Q, t in T]>=0)
@variable(CEM, vH2StoDis[s in Q, t in T]>=0)
@variable(CEM, vH2StoSOC[s in Q, t in 0:length(T)]>=0)

# HSC Transmission DV
@variable(CEM, vNewH2Pipe[i in I]>=0, Int)
@variable(CEM, vRetH2Pipe[i in I]>=0, Int)
@variable(CEM, vNewH2PipeCompCap[i in I]>=0, Int)
@variable(CEM, vRetH2PipeCompCap[i in I]>=0, Int)
@variable(CEM, vH2Flow[i in I, t in T])

# HSC NSD #
@variable(CEM, vH2NSD[z in Z, t in T]>=0)

@variable(CEM, vExtraEmmisionByZone[z in Z,t in T]>=0)
###################
### Expressions ###
###################

  ########################
## Power Sector Expression ##
  ########################

#----------#
#important: Remeber to consider the representative capacity in the cost for thermal units. Since the investment decisions are integer now, the unit capacity for every source is considered 1, while for thermal units is representative capacity of the cluster
#The same should be applied for land use --> for thermal units: Land use = (vNewPowGenCap[g] - vRetPowGenCap[g])* eRepPowGenCap *  pow_gen[s, :land_use_km2_p_cap]
  # Power Generation Expressions #
#@expression(CEM, eRepPowGenCap[g in G_ther], pow_gen[g, :existing_cap] / dfGen[g, :num_units])
@expression(CEM, eTotPowGenCap[g in G], pow_gen[g, :existing_cap] .+ pow_gen[g, :rep_capacity]*(vNewPowGenCap[g] .- vRetPowGenCap[g]))
@expression(CEM, eTotPowGenUnit[g in G_ther], pow_gen[g, :num_units]+vNewPowGenCap[g]-vRetPowGenCap[g])
@expression(CEM, ePowGenByZone[z in Z, t in T], sum(vPowGen[g, t] * (pow_gen[g, :zone] == z ? 1 : 0) for g in G))
@expression(CEM, eTotPowGenCapByZone[z in Z], sum(eTotPowGenCap[g]*(pow_gen[g,:zone] == z ? 1 : 0) for g in G_ther))
@expression(CEM, ePowResReqUp[z in Z, t in T], 0.1 * eTotPowGenCapByZone[z] .+ 0.05 * pow_demand[t,z])
@expression(CEM, ePowResReqDn[z in Z, t in T], 0.05 * pow_demand[t,z])
@expression(CEM, ePowGenEmiByZone[z in Z], sum(CO2_content[pow_gen[g, :fuel]] * pow_gen[g,:heat_rate_mmbtu_per_yr]*vPowGen[g,t]*(pow_gen[g, :zone]==z ? 1 : 0) for g in G_ther, t in T))
@expression(CEM, ePowGenLandUse[z in Z], sum((vNewPowGenCap[g] - vRetPowGenCap[g])*pow_gen[g, :rep_capacity]*pow_gen[g, :land_use_km2_p_cap]*(pow_gen[g,:zone]==z ? 1 : 0) for g in G))

# Power Storage Expressions #
@expression(CEM, eTotPowStoCap[s in S], dfStor[dfStor.r_id .==s, :existing_cap] .+ pow_gen[s, :rep_capacity]*(vNewPowStoCap[s] .- vRetPowStoCap[s]))
@expression(CEM, ePowStoChaByZone[z in Z, t in T], sum(vPowStoCha[s, t] * (pow_gen[s, :zone] == z ? 1 : 0) for s in S))
@expression(CEM, ePowStoDisByZone[z in Z, t in T], sum(vPowStoDis[s, t] * (pow_gen[s, :zone] == z ? 1 : 0) for s in S))
@expression(CEM, ePowStoLandUse[z in Z], sum((vNewPowStoCap[s] - vRetPowStoCap[s])*pow_gen[s, :rep_capacity]*pow_gen[s, :land_use_km2_p_cap]*(pow_gen[s,:zone]==z ? 1 : 0) for s in S))

# Power Transmission Expressions #
@expression(CEM, eTotPowTraCap[l in L], pow_lines[l, :existing_transmission_cap_mw] .+ vNewPowTraCap[l]) #No retired cap considered for transmission lines + Not var cost for power flows
@expression(CEM, eNet_Pow_Flow[z in Z,t in T], sum(Pow_Network[l,z] * vPowFlow[l,t] for l in L))
@expression(CEM, ePow_Loss_By_Zone[z in Z,t in T], sum(abs(Pow_Network[l,z]) * (1/2) *vPowFlow[l,t] * pow_lines[l, :line_loss_percentage] for l in L))

# Cost Expressions
@expression(CEM, eCostPowGenInv, sum(pow_gen[g, :inv_cost_per_mwyr] .* vNewPowGenCap[g] .* pow_gen[g, :rep_capacity] .+ dfGen[dfGen[!, :r_id] .==g, :fom_cost_per_mwyr] .* eTotPowGenCap[g] for g in G))
@expression(CEM, eCostPowGenVar, 
    sum((pow_gen[g, :vom_cost_mwh] + pow_gen[g, :heat_rate_mmbtu_per_yr] .* fuel_costs[pow_gen[g, :fuel]][t]) .* vPowGen[g,t] for g in G, t in T)    
)
@expression(CEM, eCostPowGenStart, sum(pow_gen[g, :start_cost_per_mw] .* pow_gen[g, :rep_capacity] .* vPowGenStart[g, t] for g in G_ther, t in T)) 
#For cost of Non-served demand we only consider $/MWh for each zone and will not consider demand segments
@expression(CEM, eCostPowNSD, sum(vPowNSD[z,t] .* zones[z, :voll_pow] for z in Z, t in T))
@expression(CEM, eCostPowStoInv, 
    sum(pow_gen[s, :inv_cost_per_mwhyr] .* vNewPowStoCap[s] *pow_gen[s, :rep_capacity] .+ ((pow_gen[s, :fom_cost_per_mwhyr] + pow_gen[s, :fom_cost_charge_per_mwyr]) .* eTotPowStoCap[s]) for s in S)
)
@expression(CEM, eCostPowStoVar, sum(vPowStoCha[s,t] .* pow_gen[s, :vom_cost_mwh_charge] for s in S, t in T))
@expression(CEM, eCostPowTraInv, sum(pow_lines[l, :line_reinforcement_cost_per_mwyr] .* vNewPowTraCap[l] for l in L))


  ################
## HSC Expression ##
  ################

# HSC Generation Expressions
@expression(CEM, eTotH2GenCap[h in H], hsc_gen[h, :existing_cap_tonne_p_hr] .+ hsc_gen[h, :rep_capacity]*(vNewH2GenCap[h] .- vRetH2GenCap[h]))
@expression(CEM, eTotH2GenUnit[h in H_ther], hsc_gen[h, :num_units]+vNewH2GenCap[h]-vRetH2GenCap[h])
@expression(CEM, eH2GenEvap[h in H, t in T], hsc_gen[h, :boil_off]*vH2Gen[h,t])
@expression(CEM, eH2GenEmiByZone[z in Z], sum(CO2_content[hsc_gen[h, :fuel]] * hsc_gen[h,:heat_rate_mmbtu_p_tonne]*vH2Gen[h,t]*(1-hsc_gen[h,:ccs_rate])*(hsc_gen[h, :zone]==z ? 1 : 0) for h in H_ther, t in T))
@expression(CEM, eH2GenLandUse[z in Z], sum((vNewH2GenCap[h]-vRetH2GenCap[h])*hsc_gen[h, :rep_capacity]*hsc_gen[h, :land_use_km2_p_cap]*(hsc_gen[h,:zone]==z ? 1 : 0) for h in H))

# HSC Storage Expression
@expression(CEM, eTotH2StoCap[s in Q], hsc_gen[s, :existing_cap_tonne] + hsc_gen[s, :rep_capacity]*(vNewH2StoCap[s] - vRetH2StoCap[s]))
@expression(CEM, eTotH2StoCompCap[s in Q], hsc_gen[s, :existing_cap_comp_tonne_hr] + vNewH2StoCompCap[s] - vRetH2StoCompCap[s])
@expression(CEM, eH2StoLandUse[z in Z], sum((vNewH2StoCap[s]-vRetH2StoCap[s])*hsc_gen[s, :rep_capacity]* hsc_gen[s, :land_use_km2_p_cap]*(hsc_gen[s,:zone]==z ? 1 : 0) for s in Q))

# HSC Tramsmission Expression
@expression(CEM, eTotH2Pipe[i in I], hsc_pipelines[i, :existing_num_pipes] + vNewH2Pipe[i] - vRetH2Pipe[i])
@expression(CEM, eNet_H2_Flow[z in Z,t in T], sum(H2_Network[i,z] * vH2Flow[i,t] for i in I))
@expression(CEM, eH2_Loss_By_Zone[z in Z,t in T], sum(abs(H2_Network[i,z]) * (1/2) *vH2Flow[i,t] * hsc_pipelines[i, :pipe_loss_coeff] for i in I ))
@expression(CEM, eTotH2PipeCompCap[i in I], hsc_pipelines[i, :existing_comp_cap_tonne_hr] + vNewH2PipeCompCap[i] - vRetH2PipeCompCap[i])
@expression(CEM, eH2PipeLandUse[z in Z],0.5 * sum(hsc_pipelines[i, :land_use_km2_p_length]*hsc_pipelines[i, :distance]*hsc_pipelines[i, :land_use_km2_p_length]*(vNewH2Pipe[i]-vRetH2Pipe[i])*abs(H2_Network[i,z]) for i in I))

# HSC Cost Expressions
@expression(CEM, eCostH2GenInv, sum(hsc_gen[h, :inv_cost_tonne_hr_p_yr]*vNewH2GenCap[h]*hsc_gen[h, :rep_capacity] + hsc_gen[h, :fom_cost_p_tonne_p_hr_yr]*eTotH2GenCap[h] for h in H))

@expression(CEM, eCostH2StoInv, sum(hsc_gen[s, :inv_cost_tonne_p_yr]*vNewH2StoCap[s]*hsc_gen[s, :rep_capacity] #=+ hsc_gen[s, :inv_cost_comp_tonne_hr_p_yr]*vNewH2StoCompCap[s]=#
                                  + hsc_gen[s, :fom_cost_p_tonne_p_yr]*eTotH2StoCap[s] #=+ hsc_gen[s, :fom_cost_comp_tonne_hr_p_yr]*eTotH2StoCompCap[w]=# for s in Q)
)
@expression(CEM, eCostH2StoCompInv, sum((hsc_gen[s, :inv_cost_comp_tonne_hr_p_yr]*vNewH2StoCompCap[s]) + (hsc_gen[s, :fom_cost_comp_tonne_hr_p_yr]*eTotH2StoCompCap[s]) for s in Q))

@expression(CEM, eCostH2TraInv, sum(hsc_pipelines[i, :investment_cost_per_length]*hsc_pipelines[i, :distance]*vNewH2Pipe[i] +
                                    hsc_pipelines[i, :fom_per_length]*hsc_pipelines[i, :distance]*eTotH2Pipe[i] +
                                    hsc_pipelines[i, :compressor_inv_per_length]*hsc_pipelines[i, :distance]*vNewH2PipeCompCap[i] + 
                                    hsc_pipelines[i, :fom_comp_p_tonne_hr]*eTotH2PipeCompCap[i] for i in I)
)

@expression(CEM, eCostH2GenVar, 
    sum((hsc_gen[h, :vom_cost_p_tonne] + hsc_gen[h, :heat_rate_mmbtu_p_tonne] .* fuel_costs[lowercase(hsc_gen[h, :fuel])][t]) .* vH2Gen[h,t] for h in H, t in T)    
)

@expression(CEM, eCostH2GenStart, sum(hsc_gen[h, :startup_cost_p_tonne_hr] .* hsc_gen[h, :rep_capacity] .* vH2GenStart[h, t] for h in H_ther, t in T)) 

@expression(CEM, eCostH2NSD, sum(vH2NSD[z,t] .* zones[z, :voll_hsc] for z in Z, t in T))

@expression(CEM, ePowDemandHSC[t in T, z in Z], sum(hsc_gen[h, :pow_demand_mwh_p_tonne]*vH2Gen[h,t]*(hsc_gen[h, :zone]==z ? 1 : 0) for h in H) + 
                                                sum(hsc_gen[s, :h2charge_mwh_p_tonne]*vH2StoCha[s,t]*(hsc_gen[s, :zone]==z ? 1 : 0) for s in Q) +
                                                sum(hsc_pipelines[i, :comp_pow_mwh_per_tonne]*vH2Flow[i,t]*(H2_Network[i,z]==1 ? 1 : 0) for i in I)
)
@expression(CEM, eH2DemandPow[z in Z,t in T], sum(pow_gen[g, :h2_demand_tonne_p_mwh]*vPowGen[g,t]*(pow_gen[g, :zone]==z ? 1 : 0) for g in G))
@expression(CEM, pow_D[t in T, z in Z], pow_demand[t,z] .+ ePowDemandHSC[t,z])
@expression(CEM, H2_D[t in T, z in Z], h2_demand[t,z] .+ eH2DemandPow[z,t])

@expression(CEM, eZonalEmissionCap[z in Z], sum(vPowGen[g,t]*pow_gen[g, :heat_rate_mmbtu_per_yr]*CO2_content[pow_gen[g, :fuel]]*(pow_gen[g, :zone] == z ? 1 : 0) for g in G, t in T) + 
                                          sum(vH2Gen[h,t]*hsc_gen[h, :heat_rate_mmbtu_p_tonne]*CO2_content[hsc_gen[h, :fuel]] *(hsc_gen[h, :zone] == z ? 1 : 0) for h in H, t in T)
)
@expression(CEM, eEmissionCost, sum(vExtraEmmisionByZone[z,t]*zones[z, :emission_cost] for z in Z, t in T))
#Defining Objective Function
obj = (eCostPowGenInv .+ eCostPowStoInv .+ eCostPowTraInv .+ eCostPowGenVar .+ eCostPowStoVar .+ eCostPowNSD .+ eCostPowGenStart).+ (eCostH2GenInv .+ eCostH2StoInv .+ eCostH2TraInv .+ eCostH2GenVar .+ eCostH2NSD .+ eCostH2StoCompInv .+ eCostH2GenStart) .+ eEmissionCost

@objective(CEM, Min, obj) # The objective is linear exept the start cost of thermal generators which is integer


###################
### Constraints ###
###################

  #------------------------ #
## Power Sector Constraints ##
  #-------------------------#

# Power Balance #
@constraint(CEM, cPowerBalance[z in Z, t in T],
    ePowGenByZone[z,t] .- eNet_Pow_Flow[z,t] .- ePow_Loss_By_Zone[z,t] .+ ePowStoDisByZone[z,t] .- ePowStoChaByZone[z,t] .+ vPowNSD[z,t] == pow_D[t,z]
)

# Power Generation #
@constraint(CEM, cMaxPowGenRetCapTher[g in G_ther], vRetPowGenCap[g]*pow_gen[g, :rep_capacity] <= pow_gen[g, :existing_cap])

for g in G
  if 0 <= pow_gen[g, :max_cap_mw]
    @constraint(CEM, cMaxPowGenCap[g in G], eTotPowGenCap[g] <= pow_gen[g, :max_cap_mw])
  end
end

@constraint(CEM, cPowGenTotCapMin[g in G], pow_gen[g, :min_cap] <= eTotPowGenCap[g])
@constraint(CEM, cMaxPowGen[g in G_ren, t in T], vPowGen[g,t] .- eTotPowGenCap[g] .* pow_gen_var[t, pow_gen[g, :resource]] <= 0)
@constraint(CEM, cMaxPowGenTher[g in G_ther, t in T], vPowGen[g,t] .- (pow_gen[g, :rep_capacity].*pow_gen[g,:max_op_level].*vPowGenOnline[g,t]) <= 0)
@constraint(CEM, cMinPowGenTher[g in G_ther, t in T], (pow_gen[g, :rep_capacity]*pow_gen[g,:min_op_level]*vPowGenOnline[g,t]) .- vPowGen[g,t] <= 0) 
@constraint(CEM, cPowOnlineUnits[g in G_ther, t in T], vPowGenOnline[g,t] <= eTotPowGenUnit[g])
@constraint(CEM, cPowStartLimits[g in G_ther, t in T], vPowGenStart[g,t]<= eTotPowGenUnit[g]-vPowGenOnline[g,t])
@constraint(CEM, cPowShutLimits[g in G_ther, t in T], vPowGenShut[g,t] <= vPowGenOnline[g,t])
@constraint(CEM, cPowUnitOnlineCon[g in G_ther, t in 2:length(T)], vPowGenOnline[g,t] - vPowGenOnline[g, t-1] == vPowGenStart[g,t]-vPowGenShut[g,t])
#Minimum up/donw time for thermal generators
@constraint(CEM, cMinUpTimePowGen[g in G_ther, t in T], sum(vPowGenStart[g,tt] for tt in intersect(T, (t - pow_gen[g, :up_time]):t)) <= vPowGenOnline[g,t])
@constraint(CEM, cMinDnTimePowGen[g in G_ther, t in T], sum(vPowGenShut[g,tt] for tt in intersect(T, (t - pow_gen[g, :down_time]):t)) <= eTotPowGenUnit[g] - vPowGenOnline[g,t])
#Ramp and auxilary power generation constraints
@constraint(CEM, cPowGenRampUp[g in G_ren, t in 2:length(T)], vPowGen[g,t]-vPowGen[g,t-1] .- pow_gen[g, :ramp_up] * eTotPowGenCap[g]<= 0) # For Non-Thermal units
@constraint(CEM, cPowGenRampDn[g in G_ren, t in 2:length(T)], vPowGen[g,t-1]-vPowGen[g,t] .- pow_gen[g, :ramp_dn] * eTotPowGenCap[g]<= 0)

@constraint(CEM, cTherPowGenRampDn[g in G_ther, t in 2:length(T)], # Quadratic
            vPowGen[g,t-1] .- vPowGen[g,t] .- (vPowGenOnline[g,t].-vPowGenStart[g,t])*pow_gen[g, :rep_capacity]*pow_gen[g,:ramp_dn] 
            .- pow_gen[g,:min_op_level]*pow_gen[g, :rep_capacity]*vPowGenStart[g,t] .+ max(pow_gen[g,:min_op_level],pow_gen[g,:ramp_dn])*vPowGenShut[g,t]*pow_gen[g, :rep_capacity]<=0
) #Ramp down constraint of thermal unit based on the Paper " Heterogenous Unit Clustering ... "

@constraint(CEM, cTherPowGenRampUp[g in G_ther, t in 2:length(T)],  #Quadratic
            vPowGen[g,t] .- vPowGen[g,t-1] .- (vPowGenOnline[g,t]-vPowGenStart[g,t])*pow_gen[g,:ramp_up]*pow_gen[g, :rep_capacity] .+ pow_gen[g, :min_op_level]*vPowGenShut[g,t]*pow_gen[g, :rep_capacity]
            .- max(pow_gen[g,:min_op_level], pow_gen[g,:ramp_up])*vPowGenStart[g,t]*pow_gen[g, :rep_capacity] <=0
)

#Spinning Reserve Constraints
#@constraint(CEM, cPowResUpMax[g in G_ther, t in T], vPowResUp[g,t] .+ vPowGen[g,t] .-  pow_gen[g, :rep_capacity]*pow_gen[g,:max_op_level]*vPowGenOnline[g,t] <=0)
#@constraint(CEM, cPowResDnMax[g in G_ther, t in T], vPowResDn[g,t] .+ pow_gen[g, :rep_capacity]*pow_gen[g,:min_op_level]*vPowGenOnline[g,t] .- vPowGen[g,t] <= 0) 
#@constraint(CEM, cPowResUP[g in G_ther, t in T], vPowResUp[g,t] .- eTotPowGenCap[g]*pow_gen[g,:ramp_up] <=0)
#@constraint(CEM, cPowResDn[g in G_ther, t in T], vPowResDn[g,t] .- eTotPowGenCap[g]*pow_gen[g,:ramp_dn] <=0)
#@constraint(CEM, cPowResReqUp[z in Z, t in T], ePowResReqUp[z,t] .- sum(vPowResUp[g,t] * (pow_gen[g,:zone]==z ? 1 : 0) for g in G_ther)<= 0)
#@constraint(CEM, cPowResReqDn[z in Z, t in T], ePowResReqDn[z,t] .- sum(vPowResDn[g,t] * (pow_gen[g,:zone]==z ? 1 : 0) for g in G_ther)<= 0)

# Power Storage #
@constraint(CEM, cMaxRetPowSto[s in S], vRetPowStoCap[s]*pow_gen[s, :rep_capacity] <= pow_gen[s, :existing_cap_mwh])
for s in S
  if 0 <= pow_gen[s, :max_cap_mwh]
    @constraint(CEM, cMaxPowStoCap[s in S],eTotPowStoCap[s]<= pow_gen[s, :max_cap_mwh])
  end
end
@constraint(CEM, cPowStoTotCapMin[s in S], pow_gen[s, :min_cap_mwh] .- eTotPowStoCap[s] <= 0)
@constraint(CEM, cPowStoBalanceInit[s in S, t=[0]], vPowSOC[s,t]==0)
@constraint(CEM, cPowStoBalance[s in S, t in 1:length(T)], vPowSOC[s,t] == (1-pow_gen[s,:etta_self_dis])*vPowSOC[s,t-1] + pow_gen[s,:etta_cha]*vPowStoCha[s,t] - (1/pow_gen[s,:etta_dis])*vPowStoDis[s,t])
@constraint(CEM, cPowStoMaxDis[s in S, t in T], vPowStoDis[s,t] <= pow_gen[s,:etta_dis]*vPowSOC[s,t-1])
@constraint(CEM, cPowStoMaxCha[s in S, t in T], vPowStoCha[s,t] .- eTotPowStoCap[s]<=0)

# Power Transmission #
@constraints(CEM, begin
                cMaxPowFlowOut[l in L, t in T],  vPowFlow[l,t] <= eTotPowTraCap[l]
                cMaxPowFlowIn[l in L, t in T], vPowFlow[l,t] >= -eTotPowTraCap[l]
end)

for l in L
  if 0 <= pow_lines[l, :line_max_reinforcement_mw]
    @constraint(CEM, [l in L], vNewPowTraCap[l] <= pow_lines[l, :line_max_reinforcement_mw])
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
  sum((vH2StoDis[s,t] - vH2StoCha[s,t])*(hsc_gen[s, :zone]==z ? 1 : 0) for s in Q) + vH2NSD[z,t] == h2_demand[t,z]
)


# H2 Generation constraint
@constraint(CEM, cMaxRetH2Cap[h in H], vRetH2GenCap[h]*hsc_gen[h, :rep_capacity] <= hsc_gen[h, :existing_cap_tonne_p_hr])
for h in H
  if 0 <= hsc_gen[h, :max_cap_tonne_p_hr]
    @constraint(CEM, [h in H], eTotH2GenCap[h]<= hsc_gen[h, :max_cap_tonne_p_hr])
  end
end
@constraint(CEM, cMinH2GenCap[h in H], hsc_gen[h, :min_cap_tonne_p_hr] <= eTotH2GenCap[h])
@constraint(CEM, cMaxH2GenVar[h in H_dis, t in T], vH2Gen[h,t] <= eTotH2GenCap[h]*hsc_gen[h, :max_op_level]*hsc_gen_var[t, hsc_gen[h,:resource]])
@constraint(CEM, cMinH2GenVar[h in H_dis, t in T], eTotH2GenCap[h]*hsc_gen[h,:min_op_level]  <= vH2Gen[h,t])
@constraint(CEM, cMaxH2GenTher[h in H_ther, t in T], vH2Gen[h,t] <=  hsc_gen[h, :max_op_level] * hsc_gen[h, :rep_capacity] * vH2GenOnline[h,t])
@constraint(CEM, cMinH2GenTher[h in H_ther, t in T], hsc_gen[h, :min_op_level] * hsc_gen[h, :rep_capacity] * vH2GenOnline[h,t] <= vH2Gen[h,t]) #Quadratic
@constraint(CEM, cH2OnlineUnits[h in H_ther, t in T], vH2GenOnline[h,t] <= eTotH2GenUnit[h])

@constraint(CEM, cH2StartLimits[h in H_ther, t in T], vH2GenStart[h,t]<= eTotH2GenUnit[h] - vH2GenOnline[h,t])
@constraint(CEM, cH2ShutLimits[h in H_ther, t in T], vH2GenShut[h,t] <= vH2GenOnline[h,t])
@constraint(CEM, cH2UnitOnlineCon[h in H_ther, t in 2:length(T)], vH2GenOnline[h,t] - vH2GenOnline[h, t-1] == vH2GenStart[h,t]-vH2GenShut[h,t])

@constraint(CEM, cMinUpTimeH2Gen[h in H_ther, t in T], sum(vH2GenStart[h,tt] for tt in intersect(T, (t - hsc_gen[h, :up_time]):t)) <= vH2GenOnline[h,t])
@constraint(CEM, cMinDnTimeH2Gen[h in H_ther, t in T], sum(vH2GenShut[h,tt] for tt in intersect(T, (t - hsc_gen[h, :down_time]):t)) <= eTotH2GenUnit[h] - vH2GenOnline[h,t])

#Ramp and Auxilary H2 Gen Constraints

@constraint(CEM, cH2GenRampUp[h in H_dis, t in 2:length(T)], vH2Gen[h,t] - vH2Gen[h,t-1] .- hsc_gen[h, :ramp_up_percentage] * eTotH2GenCap[h]<= 0) 
@constraint(CEM, cH2GenRampDn[h in H_dis, t in 2:length(T)], vH2Gen[h,t-1]-vH2Gen[h,t] .- hsc_gen[h, :ramp_down_percentage] * eTotH2GenCap[h]<= 0)

@constraint(CEM, cTherH2GenRampDn[h in H_ther, t in 2:length(T)],
            vH2Gen[h,t-1] .- vH2Gen[h,t] .- (vH2GenOnline[h,t].-vH2GenStart[h,t])*hsc_gen[h, :rep_capacity]*hsc_gen[h,:ramp_down_percentage] 
            .- hsc_gen[h,:min_op_level]*hsc_gen[h, :rep_capacity]*vH2GenStart[h,t] .+ max(hsc_gen[h,:min_op_level],hsc_gen[h,:ramp_down_percentage])*vH2GenShut[h,t]*hsc_gen[h, :rep_capacity] <=0
) #Ramp down constraint of thermal unit based on the Paper " Heterogenous Unit Clustering ... "

@constraint(CEM, cTherH2GenRampUp[h in H_ther, t in 2:length(T)], #Quadratic
            vH2Gen[h,t] .- vH2Gen[h,t-1] .- (vH2GenOnline[h,t]-vH2GenStart[h,t])*hsc_gen[h,:ramp_up_percentage]*hsc_gen[h, :rep_capacity] .+ hsc_gen[h, :min_op_level]*vH2GenShut[h,t]*hsc_gen[h, :rep_capacity]
            .- max(hsc_gen[h,:min_op_level], hsc_gen[h,:ramp_up_percentage])*vH2GenStart[h,t]*hsc_gen[h, :rep_capacity] <=0
)

# H2 Storage constraint
@constraint(CEM, cMaxRetH2StoCap[s in Q], vRetH2StoCap[s]*hsc_gen[s, :rep_capacity] <= hsc_gen[s, :existing_cap_tonne])
for s in Q
  if 0 <= hsc_gen[s, :max_cap_stor_tonne]
    @constraint(CEM, [s in S], eTotH2StoCap[s]<= hsc_gen[s, :max_cap_stor_tonne])
  end
end
@constraint(CEM, cMinH2StoCap[s in Q], hsc_gen[s, :min_cap_stor_tonne] <= eTotH2StoCap[s])
@constraint(CEM, cH2StoBalance[s in Q, t in T], vH2StoSOC[s,t] == (1-hsc_gen[s, :etta_self_dis])*vH2StoSOC[s,t-1] + vH2StoCha[s,t]*hsc_gen[s,:etta_cha] - (1/hsc_gen[s,:etta_dis])*vH2StoDis[s,t])
@constraint(CEM, cH2StoSOCInit[s in Q, t=[0]], vH2StoSOC[s,t]==0)
@constraint(CEM, cMaxH2StoSOC[s in Q, t in T], vH2StoSOC[s,t]<= hsc_gen[s,:h2stor_max_level]*eTotH2StoCap[s])
@constraint(CEM, cMinH2StoSOC[s in Q, t in T], hsc_gen[s,:h2stor_min_level]*eTotH2StoCap[s]<= vH2StoSOC[s,t])
@constraint(CEM, cMaxRetH2StorCompCap[s in Q], vRetH2StoCompCap[s] <= hsc_gen[s, :existing_cap_comp_tonne_hr])
@constraint(CEM, cMaxH2StoChar[s in Q,t in T], vH2StoCha[s,t] <= eTotH2StoCompCap[s])

# H2 Transmission constraints
@constraint(CEM, cMaxH2PipeNum[i in I], eTotH2Pipe[i] <= hsc_pipelines[i, :max_num_pipes])
@constraint(CEM, cMaxRetH2PipeNum[i in I], vRetH2Pipe[i] <= hsc_pipelines[i, :existing_num_pipes])
@constraints(CEM, begin 
                  cMaxH2PipeFlowOut[i in I, t in T], vH2Flow[i,t] <= eTotH2Pipe[i]*hsc_pipelines[i, :max_op_level]
                  cMaxH2PipeFlowIn[i in I, t in T], vH2Flow[i,t] >= -eTotH2Pipe[i]*hsc_pipelines[i, :max_op_level]
end)

@constraint(CEM, cMaxRetH2PipeCompCap[i in I], vRetH2PipeCompCap[i]<=hsc_pipelines[i, :existing_comp_cap_tonne_hr])
@constraints(CEM, begin
                  cMaxH2PipeFlowOutComp[i in I,t in T], vH2Flow[i,t] <= eTotH2PipeCompCap[i]
                  cMaxH2PipeFlowInComp[i in I, t in T], vH2Flow[i,t] >= -eTotH2PipeCompCap[i]
end)

# H2 NSD Constraints
@constraint(CEM, cH2NSD[z in Z, t in T], vH2NSD[z,t] <= zones[z, :hsc_nsd_share]*H2_D[t,z])

# System Emission Constraint by zone
@constraint(CEM, cZonalEmissionCap[z in Z, t in T], eZonalEmissionCap[z] - vExtraEmmisionByZone[z,t] <= zones[z, :emission_cap_tonne_co2])

# Land Use Constraint on each zone
@constraint(CEM, cLandUse[z in Z], ePowGenLandUse[z] + ePowStoLandUse[z] + eH2GenLandUse[z] + eH2StoLandUse[z] + eH2PipeLandUse[z] <= zones[z, :available_land])

#set_optimizer_attribute(CEM, "NonConvex", 2)
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