using JuMP, Gurobi, DataFrames, CSV


datadir = joinpath("/Users/rez/Documents/Engineering/Coding/Julia/MyCode")


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
for f in [pow_gen, pow_gen_var, pow_load, pow_lines, pow_network, hsc_gen, hsc_gen_var, hsc_load, hsc_pipelines, fuel, zones]
    rename!(f,lowercase.(names(f)))
end
#=
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
=#

###Constructing the  DataFrame###
dfGen = pow_gen[pow_gen[!, :gen_type].>= 1, :] # set of all generators
dfStor = pow_gen[pow_gen[!, :stor_type].> 0, :] #set of all power storages
dfH2Gen = hsc_gen[hsc_gen[!, :h_gen_type].>= 1, :] #set of all H2 Generators
dfH2Stor = hsc_gen[hsc_gen[!, :h_stor].>= 1, :]

start = findall(s -> s == "Load_MW_z1", names(pow_load))[1]
pow_D = Matrix(pow_load[1:length(T), start:start+length(Z)-1])
### Create Solving Function ###

# Defining Sets
G = dfGen.r_id
S = pow_gen[pow_gen[!, :stor_type].> 0, :r_id] #set of all power storages
Z = zones.zones
T = pow_load.time_index
K = dfGen[dfGen[!, :gen_type] .==1, :r_id]
R = dfGen[dfGen[!, :gen_type] .>1, :r_id]
L = pow_lines.network_lines
H = dfH2Gen[dfH2Gen[!, :h_gen_type].>= 1, :r_id]
W = hsc_gen[hsc_gen[!, :h_stor].>= 1, :r_id]
J = hsc_gen[hsc_gen[!, :h_gen_type].== 1, :]

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

col = findall(s -> s == "z1", names(pow_lines))[1]
Pow_Network = pow_lines[1:length(L), col:col+length(Z)-1]

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
@variable(CEM, vPowFlow[l in L, t in T]>=0)

# Non-Served Power demand
@variable(CEM, vPowNSD[z in Z, t in T]>=0)


  ########################
## HSC Decision Variables ##
  ########################

@variable(CEM, vNewHscGenCap[h in H]>=0)


###################
### Expressions ###
###################

  ########################
## Power Sector Expression ##
  ########################

# Power Generation Expressions #
@expression(CEM, eTotPowGenCap[g in G], dfGen[dfGen[!, :r_id] .==g, :existing_cap] .+ vNewPowGenCap[g] .- vRetPowGenCap[g])
@expression(CEM, eCostPowGenInv, sum(dfGen[dfGen[!, :r_id] .==g, :inv_cost_per_mwyr] .* vNewPowGenCap[g] .+ dfGen[dfGen[!, :r_id] .==g, :fom_cost_per_mwyr] .* eTotPowGenCap[g] for g in G))
@expression(CEM, eCostPowGenVar, 
    sum((pow_gen[g, :vom_cost_mwh] + pow_gen[g, :heat_rate_mmbtu_per_yr] .* fuel_costs[pow_gen[g, :fuel]][t]) .* vPowGen[g,t] for g in G, t in T)    
)
@expression(CEM, eCostPowGenStart, sum(pow_gen[g, :start_cost_per_mw] .* eTotPowGenCap[g] .* vPowGenStart[g, t] for g in K, t in T))       #We have non-linearity here - Either change the equation of use piecewise-linear approximation - if you get a non-linear problem look here!
@expression(CEM, ePowGenByZone[z in Z, t in T], sum(vPowGen[g, t] * (pow_gen[g, :zone] == z ? 1 : 0) for g in G))
@expression(CEM, eAuxPowGen[g in K, t in T], vPowGen[g,t] .- eTotPowGenCap[g]*pow_gen[g,:min_op_level]*vPowGenCommit[g,t]) #Non-linear
@expression(CEM, eTotPowGenCapByZoneTher[z in Z], sum(eTotPowGenCap[g]*(pow_gen[g,:zone] == z ? 1 : 0) for g in K))
@expression(CEM, ePowResReqUp[z in Z, t in T], 0.1 * eTotPowGenCapByZone[z] .+ 0.05 * pow_D[t,z])
@expression(CEM, ePowResReqDn[z in Z, t in T], 0.05 * pow_D[t,z])
@expression(CEM, ePowGenEmiByZone[z in Z], sum(CO2_content[pow_gen[g, :fuel]] * pow_gen[g,:heat_rate_mmbtu_per_yr]*vPowGen[g,t]*(pow_gen[g, :zone]==z ? 1 : 0) for g in K, t in T))
@expression(CEM, ePowGenLandUse[z in Z], sum((vNewPowGenCap[g] - vRetPowGenCap[g])*pow_gen[g, :land_use]*(pow_gen[g,:zone]==z ? 1 : 0) for g in G))

# Power Storage Expressions #
@expression(CEM, eTotPowStoCap[s in S], dfStor[dfStor.r_id .==s, :existing_cap] .+ vNewPowStoCap[s] .- vRetPowStoCap[s])
@expression(CEM, eCostPowStoInv, 
    sum(dfStor[dfStor.r_id .==s, :inv_cost_per_mwhyr] .* vNewPowStoCap[s] .+ ((dfStor[dfStor.r_id .==s, :fom_cost_per_mwhyr] + dfStor[dfStor.r_id .==s, :fom_cost_charge_per_mwyr]) .* eTotPowStoCap[s]) for s in S)
)
@expression(CEM, eCostPowStoVar, sum((vPowStoCha[s,t] + vPowStoDis[s,t]) .* pow_gen[s, :vom_cost_mwh_charge] for s in S, t in T))
@expression(CEM, ePowStoChaByZone[z in Z, t in T], sum(vPowStoCha[s, t] * (pow_gen[s, :zone] == z ? 1 : 0) for s in S))
@expression(CEM, ePowStoDisByZone[z in Z, t in T], sum(vPowStoDis[s, t] * (pow_gen[s, :zone] == z ? 1 : 0) for s in S))
@expression(CEM, ePowStoLandUse[z in Z], sum((vNewPowStoCap[s] - vRetPowStoCap[s])*pow_gen[s, :land_use]*(pow_gen[s,:zone]==z ? 1 : 0) for s in S))

# Power Transmission Expressions #
@expression(CEM, eTotPowTraCap[l in L], pow_lines[l, :existing_transmission_cap_mw] + vNewPowTraCap[l])
@expression(CEM, eCostPowTraInv, sum(pow_lines[l, :line_reinforcement_cost_per_mwyr] .* vNewPowTraCap[l] for l in L))
@expression(CEM, eNet_Export_Pow_Flow[z in Z,t in T], sum(Pow_Network[l,z] * vPowFlow[l,t] for l in L))
@expression(CEM, ePow_Loss_By_Zone[z in Z,t in T], sum(abs(Pow_Network[l,z]) * (1/2) *vPowFlow[l,t] * pow_lines[l, :line_loss_percentage] for l in L))

# Power NSD expression #
#For cost of Non-served demand we only consider $/MWh for each zone and will not consider demand segments
@expression(CEM, eCostPowNSD, sum(vPowNSD[z,t] .* zones[z, :voll_pow] for z in Z, t in T))

#Defining Objective Function
obj_pow = eCostPowGenInv .+ eCostPowStoInv .+ eCostPowTraInv .+ eCostPowGenVar .+ eCostPowStoVar .+ eCostPowNSD


############################
### Defining Constraints ###
############################

  ##########################
## Power Sector Constraints ##
  #########################

# Power Balance #
@constraint(CEM, cPowerBalance[z in Z, t in T],
    ePowGenByZone[z,t] - eNet_Export_Pow_Flow[z,t] - ePow_Loss_By_Zone[z,t] + ePowStoDisByZone[z,t] - ePowStoChaByZone[z,t] + vPowNSD[z,t] == pow_D[t,z]
)

# Power Generation #
@constraint(CEM, cMaxPowGenRetCap[g in G], vRetPowGenCap[g] <= pow_gen[g, :existing_cap])
@constraint(CEM, cPowGenTotCapMax[g in G], eTotPowGenCap[g] .- pow_gen[g, :max_cap_mw] <= 0)
@constraint(CEM, cPowGenTotCapMin[g in G], pow_gen[g, :min_cap] <= eTotPowGenCap[g])
@constraint(CEM, cMaxPowGen[g in G, t in T], vPowGen[g,t] .- eTotPowGenCap[g] .* pow_gen_var[t, pow_gen[g, :resource]] <= 0)

@constraint(CEM, cMaxPowGenTher[g in K, t in T], vPowGen[g,t] .- (eTotPowGenCap[g]*pow_gen[g,:max_op_level]*vPowGenCommit[g,t]) <= 0)
@constraint(CEM, cMinPowGenTher[g in K, t in T], (eTotPowGenCap[g]*pow_gen[g,:min_op_level]*vPowGenCommit[g,t]) .- vPowGen[g,t] <= 0) ## Non-linear
@constraint(CEM, cMinUpTimePowGen[g in K, t in T], sum(vPowGenStart[g,tt] for tt in intersect(T, (t - pow_gen[g, :up_time]):t)) <= vPowGenCommit[g,t])
@constraint(CEM, cMinDnTimePowGen[g in K, t in T], sum(vPowGenShut[g,tt] for tt in intersect(T, (t - pow_gen[g, :down_time]):t)) <= 1 - vPowGenCommit[g,t])
@constraint(CEM, cPowGenCommitCon[g in K, t in T[1:8759]], vPowGenCommit[g,t+1]-vPowGenCommit[g,t] == vPowGenStart[g,t+1]-vPowGenShut[g,t+1])
#Ramp and auxilary power generation constraints
@constraint(CEM, cAuxPowGenRampUp[g in K,t in T[1:8759]], eAuxPowGen[g,t+1] .- eAuxPowGen[g,t] <= eTotPowGenCap[g]*pow_gen[g,:ramp_up])
@constraint(CEM, cAuxPowGenRampDn[g in K,t in T[1:8759]], eAuxPowGen[g,t] .- eAuxPowGen[g,t+1]  <= eTotPowGenCap[g]*pow_gen[g,:ramp_dn])
@constraint(CEM, cPowGenRampUp[g in K, t in T[1:8759]], vPowGen[g,t+1] .- vPowGen[g,t] .- eTotPowGenCap[g]*pow_gen[g,:ramp_up]<= 0)
@constraint(CEM, cPowGenRampDn[g in K, t in T[1:8759]], vPowGen[g,t] .- vPowGen[g,t+1] .- eTotPowGenCap[g]*pow_gen[g,:ramp_dn]<= 0)
#Spinning Reserve Constraints
@constraint(CEM, cPowResUpMax[g in K, t in T], vPowResUp[g,t] .+ vPowGen[g,t] .-  eTotPowGenCap[g]*pow_gen[g,:max_op_level]*vPowGenCommit[g,t] <=0) #Non-linear
@constraint(CEM, cPowResDnMax[g in K, t in T], vPowResDn[g,t] .+ eTotPowGenCap[g]*pow_gen[g,:min_op_level]*vPowGenCommit[g,t] .- vPowGen[g,t] <= 0) #Non-linear
@constraint(CEM, cPowResUP[g in K, t in T], vPowResUp[g,t] .- eTotPowGenCap[g]*pow_gen[g,:ramp_up] <=0)
@constraint(CEM, cPowResDn[g in K, t in T], vPowResDn[g,t] .- eTotPowGenCap[g]*pow_gen[g,:ramp_dn] <=0)
@constraint(CEM, cPowResReqUp[z in Z, t in T], ePowResReqUp[z,t] .- sum(vPowResUp[g,t] * (pow_gen[g,:zone]==z ? 1 : 0) for g in K)<= 0)
@constraint(CEM, cPowResReqDn[z in Z, t in T], ePowResReqDn[z,t] .- sum(vPowResDn[g,t] * (pow_gen[g,:zone]==z ? 1 : 0) for g in K)<= 0)

# Power Storage #
@constraint(CEM, cMaxPowStorRetCap[s in S], vRetPowStoCap[s] <= pow_gen[s, :existing_cap])
@constraint(CEM, cPowStoTotCapMax[s in S], eTotPowStoCap[s] .- pow_gen[s, :max_cap_mwh] <= 0)
@constraint(CEM, cPowStoTotCapMin[s in S], pow_gen[s, :min_cap_mwh] <= eTotPowStoCap[s])
@constraint(CEM, cPowStoBalanceInit[s in S, t=0], vPowSOC[s,t]==0)
@constraint(CEM, cPowStoBalance[s in S, t in 1:length(T)], vPowSOC[s,t] == (1-pow_gen[s,:etta_self_dis])*vPowSOC[s,t-1] + pow_gen[s,:etta_cha]*vPowStoCha[s,t] - (1/pow_gen[s,:etta_dis])*vPowStoDis[s,t])
@constraint(CEM, cPowStoMaxDis[s in S, t in T], vPowStoDis[s,t] <= pow_gen[s,:etta_dis]*vPowSOC[s,t-1])

# Power Transmission #
@constraint(CEM, cMaxPowFlow[l in L, t in T], vPowFlow[l,t] <= eTotPowTraCap[l])

# Power Non-served Demand #
@constraint(CEM, cPowNSD[z in Z, t in T], vPowNSD[z,t] <= zones[z, :pow_nsd_share]*pow_D[t,z])

# Power Emissions # Method 1  = Emission should not exceed specified emission cap
@constraint(CEM, cPowEmission[z in Z], ePowGenEmiByZone[z] <= zones[z, :emission_cap_tonne_co2])

#method 2 = We can produced as much emission as we want, however we will be fined if we cross the cap
#obj_pow += sum((ePowGenEmiByZone[z] - zones[z, :emission_cap_tonne_co2]) * zones[z, :emission_cost] for z in Z)


# Power Sector Land constraint #
ePowGenLandUse[z] + ePowStoLandUse[z]



CEM
