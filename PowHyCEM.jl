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

###Constructing the  DataFrame###
dfGen = pow_gen[pow_gen[!, :gen_type].>= 1, :] # set of all generators
dfStor = pow_gen[pow_gen[!, :stor_type].> 0, :] #set of all power storages
dfH2Gen = hsc_gen[hsc_gen[!, :h_gen_type].>= 1, :] #set of all H2 Generators
dfH2Stor = hsc_gen[hsc_gen[!, :h_stor].>= 1, :]

# Defining Sets
G = dfGen.r_id
S = pow_gen[pow_gen[!, :stor_type].> 0, :r_id] #set of all power storages
Z = zones.zones
W = collect(1:52)
hours_per_week = 168
H_w = [((w - 1) * hours_per_week + 1):(w * hours_per_week) for w in W]
T = collect(1:168)
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
pow_demand = Dict{Int, Matrix}()
for w in W
  pow_demand[w] = Matrix(pow_load[H_w[w], start_p:start_p + length(Z)-1])
end

start_h = findall(s -> s == "load_hsc_tonne_z1", names(hsc_load))[1]
h2_demand = Dict{Int, Matrix}()
for w in W
  h2_demand[w] = Matrix(hsc_load[H_w[w], start_h:start_h + length(Z)-1])
end

fuels = names(fuel)[2:end]
costs = Matrix(fuel[2:end, 2:end])
CO2_content = fuel[1, 2:end] # tons CO2/MMBtu
fuel_costs = Dict{AbstractString, Array{Float64}}()
fuel_CO2 = Dict{AbstractString, Float64}()
fuel_type = dfGen[!,:fuel]

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
set_optimizer_attribute(CEM, "Method", 2)    
set_optimizer_attribute(CEM, "Crossover", 0)
#set_optimizer_attribute(CEM, "OutputFlag", 1)
#set_optimizer_attribute(CEM, "LogToConsole", 1)
set_optimizer_attribute(CEM,"MIPGap",1e-3)
set_optimizer_attribute(CEM, "BarConvTol", 1e-3)
set_optimizer_attribute(CEM, "OptimalityTol", 1e-3)

####################################
#### Defining Decision Variables ###
####################################

  #--------------------------------#
## Power Sector Decision Variables ##
  #--------------------------------#

# Power Generation DV #
@variable(CEM, vNewPowGenCap[g in G]>=0, Int)
@variable(CEM, vRetPowGenCap[g in G]>=0, Int)
@variable(CEM, vPowGen[g in G, w in W, t in T]>=0)
@variable(CEM, vPowGenOnline[g in G_ther, w in W, t in T]>=0)
@variable(CEM, vPowGenStart[g in G_ther, w in W, t in T]>=0)
@variable(CEM, vPowGenShut[g in G_ther, w in W, t in T]>=0)
@variable(CEM, vPowResUp[g in G_ther, w in W, t in T]>=0)
@variable(CEM, vPowResDn[g in G_ther, w in W, t in T]>=0)

# Power Storage DV 
@variable(CEM, vNewPowStoCap[s in S]>=0, Int)
@variable(CEM, vRetPowStoCap[s in S]>=0, Int)
@variable(CEM, vPowStoCha[s in S, w in W, t in T]>=0)
@variable(CEM, vPowStoDis[s in S, w in W, t in T]>=0)
@variable(CEM, vPowSOC[s in S, w in W, t in T]>=0)

# Power Transmission DV #
@variable(CEM, vNewPowTraCap[l in L]>=0)
@variable(CEM, vPowFlow[l in L, w in W, t in T])

  ########################
## HSC Decision Variables ##
  ########################

# HSC Generation DV
@variable(CEM, vNewH2GenCap[h in H]>=0, Int)
@variable(CEM, vRetH2GenCap[h in H]>=0, Int)
@variable(CEM, vH2Gen[h in H, w in W, t in T]>=0)
@variable(CEM, vH2GenStart[h in H_ther,  w in W, t in T]>=0)
@variable(CEM, vH2GenShut[h in H_ther,  w in W, t in T]>=0)
@variable(CEM, vH2GenOnline[h in H_ther,  w in W, t in T]>=0)

# HSC Storage DV
@variable(CEM, vNewH2StoCap[s in Q]>=0, Int)
@variable(CEM, vRetH2StoCap[s in Q]>=0, Int)
@variable(CEM, vNewH2StoCompCap[s in Q]>=0, Int)
@variable(CEM, vRetH2StoCompCap[s in Q]>=0, Int)
@variable(CEM, vH2StoCha[s in Q, w in W, t in T]>=0)
@variable(CEM, vH2StoDis[s in Q, w in W, t in T]>=0)
@variable(CEM, vH2StoSOC[s in Q, w in W, t in T]>=0)

# HSC Transmission DV
@variable(CEM, vNewH2Pipe[i in I]>=0, Int)
@variable(CEM, vRetH2Pipe[i in I]>=0, Int)
@variable(CEM, vNewH2PipeCompCap[i in I]>=0, Int)
@variable(CEM, vRetH2PipeCompCap[i in I]>=0, Int)
@variable(CEM, vH2FlowPos[i in I, w in W, t in T]>=0) #Positive flow
@variable(CEM, vH2FlowNeg[i in I, w in W, t in T]>=0) #Negative flow

# Policy Variables
@variable(CEM, vH2NSD[z in Z,  w in W, t in T]>=0)
@variable(CEM, vPowNSD[z in Z, w in W, t in T]>=0)
@variable(CEM, vPowCrt[z in Z, w in W, t in T]>=0)
@variable(CEM, vH2Crt[z in Z, w in W, t in T]>=0)
@variable(CEM, vExtraEmmisionByWeek[w in W]>=0)
@variable(CEM, vMaxEmissionByWeek[w in W]>=0)

#cyclic variables
#We define the Monolithic model the way the temporal BD with and without LDS is defined
@variable(CEM, vPowGenFirst[g in G, w in W]>=0)
@variable(CEM, vPowGenOnlineFirst[g in G_ther, w in W]>=0)
@variable(CEM, vPowSOCFirst[s in S, w in W]>=0)
@variable(CEM, vH2GenFirst[h in H, w in W]>=0)
@variable(CEM, vH2GenOnlineFirst[h in H_ther, w in W]>=0)
@variable(CEM, vH2StoSOCFirst[s in Q, w in W]>=0)
@variable(CEM, vH2StoSOCLast[s in Q, w in W]>=0)

###################
### Expressions ###
###################

  ########################
## Power Sector Expression ##
  ########################

#----------#
#important: Since the investment decisions are integer now, the unit capacity for every source is considered 1, while for thermal units is representative capacity of the cluster
#The same should be applied for land use --> for thermal units: Land use = (vNewPowGenCap[g] - vRetPowGenCap[g])* eRepPowGenCap *  pow_gen[s, :land_use_km2_p_cap]
  # Power Generation Expressions #
@expression(CEM, eTotPowGenCap[g in G], pow_gen[g, :existing_cap] .+ pow_gen[g, :rep_capacity]*(vNewPowGenCap[g] .- vRetPowGenCap[g]))
@expression(CEM, eTotPowGenUnit[g in G_ther], pow_gen[g, :num_units]+vNewPowGenCap[g]-vRetPowGenCap[g])
@expression(CEM, ePowGenByZone[z in Z, w in W, t in T], sum(vPowGen[g,w,t] * (pow_gen[g, :zone] == z ? 1 : 0) for g in G))
@expression(CEM, eTotPowGenCapByZone[z in Z], sum(eTotPowGenCap[g]*(pow_gen[g,:zone] == z ? 1 : 0) for g in G_ther))
#@expression(CEM, ePowGenEmiByZone[z in Z], sum(CO2_content[pow_gen[g, :fuel]] * pow_gen[g,:heat_rate_mmbtu_per_yr]*vPowGen[g,w,t]*(pow_gen[g, :zone]==z ? 1 : 0) for g in G_ther, w in W, t in T))
@expression(CEM, ePowGenLandUse[z in Z], sum((vNewPowGenCap[g] - vRetPowGenCap[g])*pow_gen[g, :rep_capacity]*pow_gen[g, :land_use_km2_p_cap]*(pow_gen[g,:zone]==z ? 1 : 0) for g in G))

# Power Storage Expressions #
@expression(CEM, eTotPowStoCap[s in S], dfStor[dfStor.r_id .==s, :existing_cap] .+ pow_gen[s, :rep_capacity]*(vNewPowStoCap[s] .- vRetPowStoCap[s]))
@expression(CEM, ePowStoChaByZone[z in Z, w in W, t in T], sum(vPowStoCha[s,w,t] * (pow_gen[s, :zone] == z ? 1 : 0) for s in S))
@expression(CEM, ePowStoDisByZone[z in Z, w in W, t in T], sum(vPowStoDis[s,w,t] * (pow_gen[s, :zone] == z ? 1 : 0) for s in S))
@expression(CEM, ePowStoLandUse[z in Z], sum((vNewPowStoCap[s] - vRetPowStoCap[s])*pow_gen[s, :rep_capacity]*pow_gen[s, :land_use_km2_p_cap]*(pow_gen[s,:zone]==z ? 1 : 0) for s in S))

# Power Transmission Expressions #
@expression(CEM, eTotPowTraCap[l in L], pow_lines[l, :existing_transmission_cap_mw] .+ vNewPowTraCap[l]) #No retired cap considered for transmission lines + Not var cost for power flows - Rep Cap is cosidered 1
@expression(CEM, eNet_Pow_Flow[z in Z, w in W, t in T], sum(Pow_Network[l,z] * vPowFlow[l,w,t] for l in L))
@expression(CEM, ePow_Loss_By_Zone[z in Z, w in W, t in T], sum(abs(Pow_Network[l,z]) * (1/2) *vPowFlow[l,w,t] * pow_lines[l, :line_loss_percentage] for l in L))


  ################
## HSC Expression ##
  ################

# HSC Generation Expressions
@expression(CEM, eTotH2GenCap[h in H], hsc_gen[h, :existing_cap_tonne_p_hr] .+ hsc_gen[h, :rep_capacity]*(vNewH2GenCap[h] .- vRetH2GenCap[h]))
@expression(CEM, eTotH2GenUnit[h in H_ther], hsc_gen[h, :num_units]+vNewH2GenCap[h]-vRetH2GenCap[h])
#@expression(CEM, eRepH2GenCap[h in H_ther], hsc_gen[h, :existing_cap_tonne_p_hr]/hsc_gen[h, :num_units])
@expression(CEM, eH2GenEvap[h in H, w in W, t in T], hsc_gen[h, :boil_off]*vH2Gen[h,w,t])
#@expression(CEM, eH2GenEmiByZone[z in Z], sum(CO2_content[hsc_gen[h, :fuel]] * hsc_gen[h,:heat_rate_mmbtu_p_tonne]*vH2Gen[h,w,t]*(1-hsc_gen[h,:ccs_rate])*(hsc_gen[h, :zone]==z ? 1 : 0) for h in H_ther, w in W,t in T))
@expression(CEM, eH2GenLandUse[z in Z], sum((vNewH2GenCap[h]-vRetH2GenCap[h])*hsc_gen[h, :rep_capacity]*hsc_gen[h, :land_use_km2_p_cap]*(hsc_gen[h,:zone]==z ? 1 : 0) for h in H))

# HSC Storage Expression
@expression(CEM, eTotH2StoCap[s in Q], hsc_gen[s, :existing_cap_tonne] + hsc_gen[s, :rep_capacity]*(vNewH2StoCap[s] - vRetH2StoCap[s]))
@expression(CEM, eTotH2StoCompCap[s in Q], hsc_gen[s, :existing_cap_comp_tonne_hr] + vNewH2StoCompCap[s] - vRetH2StoCompCap[s]) #rep cap is considered 1
@expression(CEM, eH2StoLandUse[z in Z], sum((vNewH2StoCap[s]-vRetH2StoCap[s])*hsc_gen[s, :rep_capacity]* hsc_gen[s, :land_use_km2_p_cap]*(hsc_gen[s,:zone]==z ? 1 : 0) for s in Q))

# HSC Tramsmission Expression
@expression(CEM, eH2FlowNet[i in I, w in W, t in T], vH2FlowPos[i,w,t] - vH2FlowNeg[i,w,t]) #Net flow of H2 in the network
@expression(CEM, eTotH2Pipe[i in I], hsc_pipelines[i, :existing_num_pipes] + vNewH2Pipe[i] - vRetH2Pipe[i])
@expression(CEM, eNet_H2_Flow[z in Z,w in W, t in T], sum(H2_Network[i,z] * eH2FlowNet[i,w,t] for i in I))
@expression(CEM, eH2_Loss_By_Zone[z in Z, w in W, t in T], sum(abs(H2_Network[i,z]) * (1/2) *eH2FlowNet[i,w,t] * hsc_pipelines[i, :pipe_loss_coeff] for i in I ))
@expression(CEM, eTotH2PipeCompCap[i in I], hsc_pipelines[i, :existing_comp_cap_tonne_hr] + vNewH2PipeCompCap[i] - vRetH2PipeCompCap[i])
@expression(CEM, eH2PipeLandUse[z in Z],0.5 * sum(hsc_pipelines[i, :land_use_km2_p_length]*hsc_pipelines[i, :distance]*hsc_pipelines[i, :land_use_km2_p_length]*(vNewH2Pipe[i]-vRetH2Pipe[i])*abs(H2_Network[i,z]) for i in I))

# HSC Cost Expressions
@expression(CEM, eCostH2GenInv, sum(hsc_gen[h, :inv_cost_tonne_hr_p_yr]*vNewH2GenCap[h]*hsc_gen[h, :rep_capacity] + hsc_gen[h, :fom_cost_p_tonne_p_hr_yr]*eTotH2GenCap[h] for h in H))

@expression(CEM, eCostH2StoInv, sum(hsc_gen[s, :inv_cost_tonne_p_yr]*vNewH2StoCap[s]*hsc_gen[s, :rep_capacity] + hsc_gen[s, :inv_cost_comp_tonne_hr_p_yr]*vNewH2StoCompCap[s] 
                                  + hsc_gen[s, :fom_cost_p_tonne_p_yr]*eTotH2StoCap[s] + hsc_gen[s, :fom_cost_comp_tonne_hr_p_yr]*eTotH2StoCompCap[s] for s in Q)
)

@expression(CEM, eCostH2TraInv, sum(hsc_pipelines[i, :investment_cost_per_length]*hsc_pipelines[i, :distance]*vNewH2Pipe[i] +
                                    hsc_pipelines[i, :fom_per_length]*hsc_pipelines[i, :distance]*eTotH2Pipe[i] +
                                    hsc_pipelines[i, :compressor_inv_per_length]*hsc_pipelines[i, :distance]*vNewH2PipeCompCap[i] + 
                                    hsc_pipelines[i, :fom_comp_p_tonne_hr]*eTotH2PipeCompCap[i] for i in I)
)

@expression(CEM, eCostH2GenVar, 
    sum((hsc_gen[h, :vom_cost_p_tonne] + hsc_gen[h, :heat_rate_mmbtu_p_tonne] .* fuel_costs[lowercase(hsc_gen[h, :fuel])][t]) .* vH2Gen[h,w,t] for h in H, w in W, t in T)    
)

@expression(CEM, eCostH2StoVar, sum(hsc_gen[s, :var_om_cost_charge_p_tonne]*vH2StoCha[s,w,t] for s in Q, w in W, t in T))
@expression(CEM, eCostH2GenStart, sum(hsc_gen[h, :startup_cost_p_tonne_hr] .* hsc_gen[h, :rep_capacity] .* vH2GenStart[h,w,t] for h in H_ther, w in W, t in T)) 

@expression(CEM, eCostH2TraVar, sum(hsc_pipelines[i, :vom_per_tonne]*(vH2FlowPos[i,w,t]+vH2FlowNeg[i,w,t]) for i in I, w in W, t in T))

@expression(CEM, eCostH2NSD, sum(vH2NSD[z,w,t] .* zones[z, :voll_hsc] for z in Z, w in W,t in T))
@expression(CEM, eCostH2Crt, sum(vH2Crt[z,w,t] .* zones[z, :h2_curtail_cost] for z in Z, w in W, t in T))
# Cost Expressions power sector
@expression(CEM, eCostPowGenInv, sum(pow_gen[g, :inv_cost_per_mwyr] .* vNewPowGenCap[g] .* pow_gen[g, :rep_capacity] .+ pow_gen[g, :fom_cost_per_mwyr] .* eTotPowGenCap[g] for g in G))
@expression(CEM, eCostPowGenVar, 
    sum((pow_gen[g, :vom_cost_mwh] + pow_gen[g, :heat_rate_mmbtu_per_mwh] .* fuel_costs[pow_gen[g, :fuel]][t]) .* vPowGen[g,w,t] for g in G, w in W, t in T)    
)
@expression(CEM, eCostPowGenStart, sum(pow_gen[g, :start_cost_per_mw] .* pow_gen[g, :rep_capacity] .* vPowGenStart[g,w,t] for g in G_ther, w in W, t in T)) 
#For cost of Non-served demand we only consider $/MWh for each zone and will not consider demand segments
@expression(CEM, eCostPowNSD, sum(vPowNSD[z,w,t] .* zones[z, :voll_pow] for z in Z, w in W, t in T)) 
@expression(CEM, eCostPowCrt, sum(vPowCrt[z,w,t] .* zones[z, :pow_curtail_cost] for z in Z, w in W, t in T))

@expression(CEM, eCostPowStoInv, 
    sum(pow_gen[s, :inv_cost_per_mwhyr] .* vNewPowStoCap[s] *pow_gen[s, :rep_capacity] .+ (pow_gen[s, :fom_cost_per_mwhyr] .* eTotPowStoCap[s]) for s in S)
) 
@expression(CEM, eCostPowStoVar, sum(vPowStoCha[s,w,t] .* pow_gen[s, :vom_cost_mwh_charge] for s in S, w in W, t in T))
@expression(CEM, eCostPowTraInv, sum(pow_lines[l, :line_reinforcement_cost_per_mwyr] .* vNewPowTraCap[l] for l in L))

#Power and H2 Demads
@expression(CEM, ePowDemandHSC[w in W, t in T, z in Z], 
            sum(hsc_gen[h, :pow_demand_mwh_p_tonne]*vH2Gen[h,w,t]*(hsc_gen[h, :zone]==z ? 1 : 0) for h in H) + 
            sum(hsc_gen[s, :h2charge_mwh_p_tonne]*vH2StoCha[s,w,t]*(hsc_gen[s, :zone]==z ? 1 : 0) for s in Q)
)
@expression(CEM, eH2DemandPow[w in W, t in T, z in Z], 
    sum(pow_gen[g, :h2_demand_tonne_p_mwh]*vPowGen[g,w,t]*(pow_gen[g, :zone]==z ? 1 : 0) for g in G)
)
@expression(CEM, pow_D[w in W, t in T, z in Z], pow_demand[w][t,z] .+ ePowDemandHSC[w,t,z])
@expression(CEM, H2_D[w in W, t in T, z in Z], h2_demand[w][t,z] .+ eH2DemandPow[w,t,z])

@expression(CEM, eEmissionByWeek[w in W], sum(vPowGen[g,w,t]*pow_gen[g, :heat_rate_mmbtu_per_mwh]*CO2_content[pow_gen[g, :fuel]] for g in G, t in T) + 
                                          sum(vH2Gen[h,w,t]*hsc_gen[h, :heat_rate_mmbtu_p_tonne]*CO2_content[hsc_gen[h, :fuel]] for h in H, t in T) +
                                          sum(vPowGenStart[g,w,t]*pow_gen[g, :heat_rate_mmbtu_per_mwh]*CO2_content[pow_gen[g, :fuel]] for g in G_ther, t in T) +
                                          sum(vH2GenStart[h,w,t]*hsc_gen[h, :heat_rate_mmbtu_p_tonne]*CO2_content[hsc_gen[h, :fuel]] for h in H_ther, t in T)
)

@expression(CEM, ePowResReqUp[z in Z, w in W, t in T], 0.1 * eTotPowGenCapByZone[z])
@expression(CEM, ePowResReqDn[z in Z, w in W, t in T], 0.05 * eTotPowGenCapByZone[z])

#Defining Objective Function
@expression(CEM, eEmissionCost, sum(vExtraEmmisionByWeek[w]*zones[1,:emission_cost] for w in W))

obj = (eCostPowGenInv .+ eCostPowStoInv .+ eCostPowTraInv .+ eCostPowGenVar .+ eCostPowStoVar .+ eCostPowNSD .+ eCostPowCrt .+ eCostPowGenStart).+ 
      (eCostH2GenInv .+ eCostH2StoInv .+ eCostH2TraInv .+ eCostH2GenVar .+ eCostH2StoVar.+ eCostH2TraVar .+ eCostH2NSD .+ eCostH2Crt.+ eCostH2GenStart) .+ eEmissionCost

@objective(CEM, Min, obj) 


###################
### Constraints ###
###################

  #------------------------ #
## Power Sector Constraints ##
  #-------------------------#

# Power Balance #
@constraint(CEM, cPowerBalance[z in Z, w in W,t in T],
    ePowGenByZone[z,w,t] .+ eNet_Pow_Flow[z,w,t] .- ePow_Loss_By_Zone[z,w,t] .+ ePowStoDisByZone[z,w,t] .- ePowStoChaByZone[z,w,t] .+ vPowNSD[z,w,t] .- vPowCrt[z,w,t] == pow_D[w,t,z]
)

# Power Generation #
@constraint(CEM, cMaxPowGenRetCapTher[g in G], vRetPowGenCap[g]*pow_gen[g, :rep_capacity] <= pow_gen[g, :existing_cap])

@constraint(CEM, cMaxPowGenCap[g in G], eTotPowGenCap[g] <= pow_gen[g, :max_cap_mw])

@constraint(CEM, cPowGenTotCapMin[g in G], pow_gen[g, :min_cap] <= eTotPowGenCap[g])
@constraint(CEM, cMaxPowGen[g in G_ren, w in W, t in T], vPowGen[g,w,t] .- eTotPowGenCap[g] .* pow_gen_var[H_w[w][t], pow_gen[g, :resource]] <= 0)
@constraint(CEM, cMaxPowGenTher[g in G_ther, w in W, t in T], vPowGen[g,w,t] .- (pow_gen[g, :rep_capacity].*pow_gen[g,:max_op_level].*vPowGenOnline[g,w,t]) <= 0)
@constraint(CEM, cMinPowGenTher[g in G_ther,w in W,t in T], (pow_gen[g, :rep_capacity]*pow_gen[g,:min_op_level]*vPowGenOnline[g,w,t]) .- vPowGen[g,w,t] <= 0)
@constraint(CEM, cPowOnlineUnits[g in G_ther, w in W, t in T], vPowGenOnline[g,w,t] .- eTotPowGenUnit[g] <= 0)
@constraint(CEM, cPowStartLimits[g in G_ther,w in W,t in T], vPowGenStart[g,w,t] .- eTotPowGenUnit[g] .+ vPowGenOnline[g,w,t] <= 0)
@constraint(CEM, cPowShutLimits[g in G_ther,w in W,t in T], vPowGenShut[g,w,t] .- vPowGenOnline[g,w,t]<= 0)
#Cyclic constraint on Thermal power units
@constraint(CEM, cPowUnitOnlineCon[g in G_ther, w in W, t in 2:length(T)], vPowGenOnline[g,w,t] - vPowGenOnline[g,w,t-1] == vPowGenStart[g,w,t]-vPowGenShut[g,w,t])
@constraint(CEM, cPowUnitOnlineFirst[g in G_ther, w in 2:length(W)], vPowGenOnline[g,w,1] - vPowGenOnlineFirst[g,w] == vPowGenStart[g,w,1]-vPowGenShut[g,w,1])
@constraint(CEM, cPowUnitOnlineCycle[g in G_ther, w in W], vPowGenOnlineFirst[g,w] == vPowGenOnline[g,w,168])

#Ramping of non-thermal units 
@constraint(CEM, cPowGenRampUp[g in G_ren, w in W, t in 2:length(T)], vPowGen[g,w,t]-vPowGen[g,w,t-1] .- pow_gen[g, :ramp_up] * eTotPowGenCap[g]<= 0)
@constraint(CEM, cPowGenRampUpFirst[g in G_ren, w in W], vPowGen[g,w,1]-vPowGenFirst[g,w] .- pow_gen[g, :ramp_up]*eTotPowGenCap[g]<=0)
@constraint(CEM, cPowGenRampDn[g in G_ren, w in W, t in 2:length(T)], vPowGen[g,w,t-1]-vPowGen[g,w,t] .- pow_gen[g, :ramp_dn] * eTotPowGenCap[g]<= 0)
@constraint(CEM, cPowGenRampDnFirst[g in G_ren, w in W], vPowGenFirst[g,w]-vPowGen[g,w,1] .- pow_gen[g, :ramp_dn] * eTotPowGenCap[g]<= 0)
@constraint(CEM, cPowGenCycle[g in G_ren, w in W], vPowGenFirst[g,w] == vPowGen[g,w,168])

#Ramping of Thermal units
@constraint(CEM, cTherPowGenRampDn[g in G_ther, w in W, t in 2:length(T)], 
            vPowGen[g,w,t-1] .- vPowGen[g,w,t] .- pow_gen[g, :rep_capacity]*pow_gen[g,:ramp_dn] *(vPowGenOnline[g,w,t].-vPowGenStart[g,w,t])
            .+ pow_gen[g, :rep_capacity]*vPowGenStart[g,w,t]*pow_gen[g,:min_op_level] .- vPowGenShut[g,w,t]*pow_gen[g, :rep_capacity]*min(pow_gen_var[H_w[w][t], pow_gen[g, :resource]],max(pow_gen[g,:min_op_level],pow_gen[g,:ramp_dn]))<=0
)
@constraint(CEM, cTherPowGenRampDnFirst[g in G_ther, w in 2:length(W)], 
            vPowGenFirst[g,w] .- vPowGen[g,w,1] .- pow_gen[g, :rep_capacity]*pow_gen[g,:ramp_dn] *(vPowGenOnline[g,w,1].-vPowGenStart[g,w,1])
            .+ pow_gen[g, :rep_capacity]*vPowGenStart[g,w,1]*pow_gen[g,:min_op_level] .- vPowGenShut[g,w,1]*pow_gen[g, :rep_capacity]*min(pow_gen_var[H_w[w][1], pow_gen[g, :resource]],max(pow_gen[g,:min_op_level],pow_gen[g,:ramp_dn]))<=0
)

@constraint(CEM, cTherPowGenRampUp[g in G_ther, w in W, t in 2:length(T)], 
            vPowGen[g,w,t] .- vPowGen[g,w,t-1] .- pow_gen[g, :rep_capacity]*pow_gen[g,:ramp_up] *(vPowGenOnline[g,w,t].-vPowGenStart[g,w,t])
            .+ pow_gen[g, :rep_capacity]*vPowGenShut[g,w,t]*pow_gen[g,:min_op_level] .- vPowGenStart[g,w,t]*pow_gen[g, :rep_capacity]*min(pow_gen_var[H_w[w][t], pow_gen[g, :resource]],max(pow_gen[g,:min_op_level],pow_gen[g,:ramp_up]))<=0
)

@constraint(CEM, cTherPowGenRampUpFirst[g in G_ther, w in 2:length(W)], 
            vPowGen[g,w,1] .- vPowGenFirst[g,w] .- pow_gen[g, :rep_capacity]*pow_gen[g,:ramp_up] *(vPowGenOnline[g,w,1].-vPowGenStart[g,w,1])
            .- pow_gen[g, :rep_capacity]*vPowGenShut[g,w,1]*pow_gen[g,:min_op_level] .+ vPowGenStart[g,w,1]*pow_gen[g, :rep_capacity]*min(pow_gen_var[H_w[w][1], pow_gen[g, :resource]],max(pow_gen[g,:min_op_level],pow_gen[g,:ramp_up]))<=0
)


#Minimum up/donw time for thermal generators
@constraint(CEM, cMinUpTimePowGen[g in G_ther,w in W,t in T], sum(vPowGenStart[g,w,tt] for tt in intersect(T, (t - pow_gen[g, :up_time]):t)) <= vPowGenOnline[g,w,t])
@constraint(CEM, cMinDnTimePowGen[g in G_ther,w in W,t in T], sum(vPowGenShut[g,w,tt] for tt in intersect(T, (t - pow_gen[g, :down_time]):t)) <= eTotPowGenUnit[g] - vPowGenOnline[g,w,t])
#Spinning Reserve Constraints
@constraint(CEM, cPowResUpMax[g in G_ther, w in W,t in T], vPowResUp[g,w,t] .+ vPowGen[g,w,t] .-  pow_gen[g, :rep_capacity]*pow_gen[g,:max_op_level]*vPowGenOnline[g,w,t] <=0)
@constraint(CEM, cPowResDnMax[g in G_ther, w in W,t in T], vPowResDn[g,w,t] .+ pow_gen[g, :rep_capacity]*pow_gen[g,:min_op_level]*vPowGenOnline[g,w,t] .- vPowGen[g,w,t] <= 0) 
@constraint(CEM, cPowResUP[g in G_ther, w in W, t in T], vPowResUp[g,w,t] .- eTotPowGenCap[g]*pow_gen[g,:ramp_up] <=0)
@constraint(CEM, cPowResDn[g in G_ther, w in W, t in T], vPowResDn[g,w,t] .- eTotPowGenCap[g]*pow_gen[g,:ramp_dn] <=0)
@constraint(CEM, cPowResReqUp[z in Z, w in W, t in T], ePowResReqUp[z,w,t] .- sum(vPowResUp[g,w,t] * (pow_gen[g,:zone]==z ? 1 : 0) for g in G_ther)<= 0)
@constraint(CEM, cPowResReqDn[z in Z, w in W, t in T], ePowResReqDn[z,w,t] .- sum(vPowResDn[g,w,t] * (pow_gen[g,:zone]==z ? 1 : 0) for g in G_ther)<= 0)

# Power Storage #
@constraint(CEM, cMaxRetPowSto[s in S], vRetPowStoCap[s]*pow_gen[s, :rep_capacity] <= pow_gen[s, :existing_cap_mwh])

@constraint(CEM, cMaxPowStoCap[s in S], eTotPowStoCap[s] .- pow_gen[s, :max_cap_mwh]<= 0)

@constraint(CEM, cPowStoTotCapMin[s in S], pow_gen[s, :min_cap_mwh] .-  eTotPowStoCap[s] <=0)
@constraint(CEM, cPowStoBalance[s in S, w in W, t in 2:length(T)], vPowSOC[s,w,t] == (1-pow_gen[s,:etta_self_dis])*vPowSOC[s,w,t-1] + pow_gen[s,:etta_cha]*vPowStoCha[s,w,t] - (1/pow_gen[s,:etta_dis])*vPowStoDis[s,w,t])
@constraint(CEM, cPowStoBalanceFirst[s in S, w in 2:length(W)], vPowSOC[s,w,1] == (1-pow_gen[s,:etta_self_dis])*vPowSOCFirst[s,w] + pow_gen[s,:etta_cha]*vPowStoCha[s,w,1] - (1/pow_gen[s,:etta_dis])*vPowStoDis[s,w,1])
@constraint(CEM, cPowStoMaxDis[s in S, w in W, t in 2:length(T)], vPowStoDis[s,w,t] <= pow_gen[s,:etta_dis]*vPowSOC[s,w,t-1])
@constraint(CEM, cPowStoMaxDisFirst[s in S, w in 2:length(W)], vPowStoDis[s,w,1] <= pow_gen[s,:etta_dis]*vPowSOCFirst[s,w])
@constraint(CEM, cPowSOCCycle[s in S, w in W], vPowSOCFirst[s,w] == vPowSOC[s,w,168])
@constraint(CEM, cPowStoMaxCha[s in S, w in W, t in T], vPowStoCha[s,w,t] .- eTotPowStoCap[s]<=0)
@constraint(CEM, cPowStoSOCMax[s in S, w in W, t in T], vPowSOC[s,w,t] .- eTotPowStoCap[s]*pow_gen[s, :max_op_level] <=0)
@constraint(CEM, cPowStoSOCMin[s in S, w in W, t in T], eTotPowStoCap[s]*pow_gen[s, :min_op_level] .- vPowSOC[s,w,t] >=0)
# Power Transmission #
@constraints(CEM, begin
                cMaxPowFlowOut[l in L, w in W, t in T],  vPowFlow[l,w,t] <= eTotPowTraCap[l]
                cMaxPowFlowIn[l in L, w in W, t in T], vPowFlow[l,w,t] >= -eTotPowTraCap[l]
end)

@constraint(CEM, cMaxPowTraCap[l in L], vNewPowTraCap[l] <= pow_lines[l, :line_max_reinforcement_mw])


  #################
## HSC Constraints ##
  #################

# H2 Balance constraint
@constraint(CEM, cH2Balance[z in Z, w in W, t in T],
  sum((vH2Gen[h,w,t] .- eH2GenEvap[h,w,t])*(hsc_gen[h, :zone]==z ? 1 : 0) for h in H) .+ eNet_H2_Flow[z,w,t] .- eH2_Loss_By_Zone[z,w,t] .+
  sum((vH2StoDis[s,w,t] - vH2StoCha[s,w,t])*(hsc_gen[s, :zone]==z ? 1 : 0) for s in Q) + vH2NSD[z,w,t] .- vH2Crt[z,w,t]== H2_D[w,t,z]
)

# H2 Generation constraint
@constraint(CEM, cMaxRetH2Cap[h in H], vRetH2GenCap[h]*hsc_gen[h, :rep_capacity] <= hsc_gen[h, :existing_cap_tonne_p_hr])

@constraint(CEM, cMaxH2GenCap[h in H], eTotH2GenCap[h]<= hsc_gen[h, :max_cap_tonne_p_hr])

@constraint(CEM, cMinH2GenCap[h in H], hsc_gen[h, :min_cap_tonne_p_hr] <= eTotH2GenCap[h])
@constraint(CEM, cMaxH2GenVar[h in H_dis, w in W, t in T], vH2Gen[h,w,t] <= eTotH2GenCap[h]*hsc_gen[h, :max_op_level]*hsc_gen_var[H_w[w][t], hsc_gen[h, :resource]])
@constraint(CEM, cMinH2GenVar[h in H_dis, w in W, t in T], eTotH2GenCap[h]*hsc_gen[h,:min_op_level]*hsc_gen_var[H_w[w][t], hsc_gen[h, :resource]]  <= vH2Gen[h,w,t])

@constraint(CEM, cMaxH2GenTher[h in H_ther, w in W, t in T], vH2Gen[h,w,t] <=  hsc_gen[h, :max_op_level] * hsc_gen[h, :rep_capacity] * vH2GenOnline[h,w,t]*hsc_gen_var[H_w[w][t], hsc_gen[h, :resource]])
@constraint(CEM, cMinH2GenTher[h in H_ther, w in W, t in T], hsc_gen[h, :min_op_level] * hsc_gen[h, :rep_capacity] * vH2GenOnline[h,w,t]*hsc_gen_var[H_w[w][t], hsc_gen[h, :resource]] <= vH2Gen[h,w,t])
@constraint(CEM, cH2OnlineUnits[h in H_ther, w in W, t in T], vH2GenOnline[h,w,t] <= eTotH2GenUnit[h])
@constraint(CEM, cH2StartLimits[h in H_ther, w in W, t in T], vH2GenStart[h,w,t]<= eTotH2GenUnit[h] - vH2GenOnline[h,w,t])
@constraint(CEM, cH2ShutLimits[h in H_ther, w in W, t in T], vH2GenShut[h,w,t] <= vH2GenOnline[h,w,t])

@constraint(CEM, cH2UnitOnlineCon[h in H_ther, w in W, t in 2:length(T)], vH2GenOnline[h,w,t] - vH2GenOnline[h,w,t-1] == vH2GenStart[h,w,t]-vH2GenShut[h,w,t])
@constraint(CEM, cH2UnitOnlineConFirst[h in H_ther, w in W], vH2GenOnline[h,w,1] - vH2GenOnlineFirst[h,w] == vH2GenStart[h,w,1]-vH2GenShut[h,w,1])
@constraint(CEM, cH2UnitOnlineCycle[h in H_ther, w in W], vH2GenOnlineFirst[h,w] == vH2GenOnline[h,w,168])

# Min Up and Down time for Thermal H2 generators
@constraint(CEM, cMinUpTimeH2Gen[h in H_ther, w in W, t in T], sum(vH2GenStart[h,w,tt] for tt in intersect(T, (t - hsc_gen[h, :up_time]):t)) <= vH2GenOnline[h,w,t])
@constraint(CEM, cMinDnTimeH2Gen[h in H_ther, w in W, t in T], sum(vH2GenShut[h,w,tt] for tt in intersect(T, (t - hsc_gen[h, :down_time]):t)) <= eTotH2GenUnit[h] - vH2GenOnline[h,w,t])

#Ramp Constraint for Thermal units
@constraint(CEM, cTherH2GenRampDn[h in H_ther, w in W, t in 2:length(T)], 
            vH2Gen[h,w,t-1] .- vH2Gen[h,w,t] .- hsc_gen[h, :rep_capacity]*hsc_gen[h,:ramp_down_percentage] *(vH2GenOnline[h,w,t].-vH2GenStart[h,w,t])
            .+ hsc_gen[h, :rep_capacity]*vH2GenStart[h,w,t]*hsc_gen[h,:min_op_level] .- vH2GenShut[h,w,t]*hsc_gen[h, :rep_capacity]*min(hsc_gen_var[H_w[w][t], hsc_gen[h, :resource]],max(hsc_gen[h,:min_op_level],hsc_gen[h,:ramp_down_percentage]))<=0
)
@constraint(CEM, cTherH2GenRampDnFirst[h in H_ther, w in W], 
            vH2GenFirst[h,w] .- vH2Gen[h,w,1] .- hsc_gen[h, :rep_capacity]*hsc_gen[h,:ramp_down_percentage] *(vH2GenOnline[h,w,1].-vH2GenStart[h,w,1])
            .+ hsc_gen[h, :rep_capacity]*vH2GenStart[h,w,1]*hsc_gen[h,:min_op_level] .- vH2GenShut[h,w,1]*hsc_gen[h, :rep_capacity]*min(hsc_gen_var[H_w[w][1], hsc_gen[h, :resource]],max(hsc_gen[h,:min_op_level],hsc_gen[h,:ramp_down_percentage]))<=0
)

@constraint(CEM, cTherH2GenRampUp[h in H_ther, w in W, t in 2:length(T)], 
            vH2Gen[h,w,t] .- vH2Gen[h,w,t-1] .- hsc_gen[h, :rep_capacity]*hsc_gen[h,:ramp_up_percentage] *(vH2GenOnline[h,w,t].-vH2GenStart[h,w,t])
            .+ hsc_gen[h, :rep_capacity]*vH2GenShut[h,w,t]*hsc_gen[h,:min_op_level] .- vH2GenStart[h,w,t]*hsc_gen[h, :rep_capacity]*min(hsc_gen_var[H_w[w][t], hsc_gen[h, :resource]],max(hsc_gen[h,:min_op_level],hsc_gen[h,:ramp_up_percentage]))<=0
)

@constraint(CEM, cTherH2GenRampUpFirst[g in H_ther, w in W], 
            vH2Gen[g,w,1] .- vH2GenFirst[g,w] .- hsc_gen[g, :rep_capacity]*hsc_gen[g,:ramp_up_percentage] *(vH2GenOnline[g,w,1].-vH2GenStart[g,w,1])
            .- hsc_gen[g, :rep_capacity]*vH2GenShut[g,w,1]*hsc_gen[g,:min_op_level] .+ vH2GenStart[g,w,1]*hsc_gen[g, :rep_capacity]*min(hsc_gen_var[H_w[w][1], hsc_gen[g, :resource]],max(hsc_gen[g,:min_op_level],hsc_gen[g,:ramp_up_percentage]))<=0
)

# Ramp constraints for dispatachable units
@constraint(CEM, cH2GenRampUp[g in H_dis, w in W, t in 2:length(T)], vH2Gen[g,w,t]-vH2Gen[g,w,t-1] .- hsc_gen[g, :ramp_up_percentage] * eTotH2GenCap[g]<= 0) 
@constraint(CEM, cH2GenRampUpFirst[g in H_dis, w in W], vH2Gen[g,w,1]-vH2GenFirst[g,w] .- hsc_gen[g, :ramp_up_percentage]*eTotH2GenCap[g]<=0)
@constraint(CEM, cH2GenRampDn[g in H_dis, w in W, t in 2:length(T)], vH2Gen[g,w,t-1]-vH2Gen[g,w,t] .- hsc_gen[g, :ramp_down_percentage] * eTotH2GenCap[g]<= 0)
@constraint(CEM, cH2GenRampDnFirst[g in H_dis, w in W], vH2GenFirst[g,w]-vH2Gen[g,w,1] .- hsc_gen[g, :ramp_down_percentage] * eTotH2GenCap[g]<= 0)
@constraint(CEM, cH2GenCycle[g in H_dis, w in W], vH2GenFirst[g,w] == vH2Gen[g,w,168])

# H2 Storage constraint
@constraint(CEM, cMaxRetH2StoCap[s in Q], vRetH2StoCap[s]*hsc_gen[s, :rep_capacity] <= hsc_gen[s, :existing_cap_tonne])

@constraint(CEM, cMaxH2StoCap[s in Q], eTotH2StoCap[s]<= hsc_gen[s, :max_cap_stor_tonne])

@constraint(CEM, cMinH2StoCap[s in Q], hsc_gen[s, :min_cap_stor_tonne] <= eTotH2StoCap[s])
@constraint(CEM, cH2StoBalance[s in Q, w in W, t in 2:length(T)], vH2StoSOC[s,w,t] == (1-hsc_gen[s, :etta_self_dis])*vH2StoSOC[s,w,t-1] + vH2StoCha[s,w,t]*hsc_gen[s,:etta_cha] - (1/hsc_gen[s,:etta_dis])*vH2StoDis[s,w,t])
@constraint(CEM, cH2StoBalanceFirst[s in Q, w in 2:length(W)], vH2StoSOC[s,w,1] == (1-hsc_gen[s, :etta_self_dis])*vH2StoSOC[s,w-1,168] + vH2StoCha[s,w,1]*hsc_gen[s,:etta_cha] - (1/hsc_gen[s,:etta_dis])*vH2StoDis[s,w,1])
@constraint(CEM, cMaxH2StoSOC[s in Q, w in W, t in T], vH2StoSOC[s,w,t]<= hsc_gen[s,:h2stor_max_level]*eTotH2StoCap[s])
@constraint(CEM, cMinH2StoSOC[s in Q, w in W, t in T], hsc_gen[s,:h2stor_min_level]*eTotH2StoCap[s]<= vH2StoSOC[s,w,t])
@constraint(CEM, cMaxRetH2StorCompCap[s in Q], vRetH2StoCompCap[s] <= hsc_gen[s, :existing_cap_comp_tonne_hr])
@constraint(CEM, cMaxH2StorCompcCap[s in Q], eTotH2StoCompCap[s]<= eTotH2StoCap[s])
@constraint(CEM, cMinH2StorCompcCap[s in Q], 0.01*eTotH2StoCap[s] <= eTotH2StoCompCap[s])
@constraint(CEM, cMaxH2StoChar[s in Q,w in W,t in T], vH2StoCha[s,w,t] <= eTotH2StoCompCap[s])

# H2 Transmission constraints
@constraint(CEM, cMaxH2PipeNum[i in I], eTotH2Pipe[i] <= hsc_pipelines[i, :max_num_pipes])
@constraint(CEM, cMaxRetH2PipeNum[i in I], vRetH2Pipe[i] <= hsc_pipelines[i, :existing_num_pipes])
@constraints(CEM, begin 
                  cMaxH2PipeFlowOut[i in I, w in W, t in T], vH2FlowPos[i,w,t] <= eTotH2Pipe[i]*hsc_pipelines[i, :max_op_level]
                  cMaxH2PipeFlowIn[i in I, w in W, t in T], vH2FlowNeg[i,w,t] <= eTotH2Pipe[i]*hsc_pipelines[i, :max_op_level]
end)

@constraint(CEM, cMaxRetH2PipeCompCap[i in I], vRetH2PipeCompCap[i]<=hsc_pipelines[i, :existing_comp_cap_tonne_hr])
@constraints(CEM, begin
                  cMaxH2PipeFlowOutComp[i in I,w in W,t in T], vH2FlowPos[i,w,t] <= eTotH2PipeCompCap[i]
                  cMaxH2PipeFlowInComp[i in I, w in W, t in T], vH2FlowNeg[i,w,t] <= eTotH2PipeCompCap[i]
end)

# Policy constraints

# Power Non-served Demand #
#@constraint(CEM, cPowNSD[z in Z, w in W, t in T], vPowNSD[z,w,t] <= zones[z, :pow_nsd_share]*pow_D[w,t,z] )
# H2 NSD Constraints
#@constraint(CEM, cH2NSD[z in Z, w in W, t in T], vH2NSD[z,w,t] <= zones[z, :hsc_nsd_share]*H2_D[w, t,z])

#Emission constraint
@constraint(CEM, cMaxEmission, sum(vMaxEmissionByWeek[w] for w in W) <= 0.05*sum((h2_demand[w][t,z]*33.3) +pow_demand[w][t,z] for z in Z, w in W, t in T))
@constraint(CEM, cZonalEmissionCapByWeek[w in W], eEmissionByWeek[w] - vExtraEmmisionByWeek[w] <= vMaxEmissionByWeek[w])

#Land Use Constraint on each zone
@constraint(CEM, cLandUse[z in Z], ePowGenLandUse[z] + ePowStoLandUse[z] + eH2GenLandUse[z] + eH2StoLandUse[z] + eH2PipeLandUse[z] <= zones[z, :available_land])


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