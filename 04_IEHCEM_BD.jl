using JuMP, Gurobi, DataFrames, CSV
using MathOptInterface
using JuMP.MOI

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
#T = pow_load.time_index
#T = pow_load[1:8736, :time_index]
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

big_M = 100000
LB = -Inf
UB = Inf
gap = UB - LB
iteration = 0
max_iter = 1000
tolerence = 1e-6
function farkas_dual(model, con_ref)
    try
        value = MOI.get(model, MOI.DualFarkas(), con_ref)
        println("Farkas dual for $con_ref: $value")
        return value
    catch e
        println("Error getting Farkas dual for $con_ref: $e")
        return 0.0
    end
end

#Defining the Model
MP = Model(Gurobi.Optimizer)

set_optimizer_attribute(MP, "InfUnbdInfo", 1)
set_optimizer_attribute(MP, "OutputFlag", 1)
set_optimizer_attribute(MP, "DualReductions", 0)
 
# ---- MP Variables ---- #
@variable(MP, vNewPowGenCap[g in G]>=0, Int)
@variable(MP, vRetPowGenCap[g in G]>=0, Int)
@variable(MP, vNewPowStoCap[s in S]>=0, Int)
@variable(MP, vRetPowStoCap[s in S]>=0, Int)
@variable(MP, vNewPowTraCap[l in L]>=0, Int)

@variable(MP, vNewH2GenCap[h in H]>=0, Int)
@variable(MP, vRetH2GenCap[h in H]>=0, Int)
@variable(MP, vNewH2StoCap[s in Q]>=0, Int)
@variable(MP, vRetH2StoCap[s in Q]>=0, Int)
@variable(MP, vNewH2StoCompCap[s in Q]>=0, Int)
@variable(MP, vRetH2StoCompCap[s in Q]>=0, Int)
@variable(MP, vNewH2Pipe[i in I]>=0, Int)
@variable(MP, vRetH2Pipe[i in I]>=0, Int)
@variable(MP, vNewH2PipeCompCap[i in I]>=0, Int)
@variable(MP, vRetH2PipeCompCap[i in I]>=0, Int)
@variable(MP, etta >=0)
@variable(MP, vMaxEmissionByWeek[z in Z, w in W]>=0)

# ---- MP Expressions ---- #
@expression(MP, eTotPowGenCap[g in G], pow_gen[g, :existing_cap] .+ pow_gen[g, :rep_capacity]*(vNewPowGenCap[g] .- vRetPowGenCap[g]))
@expression(MP, eTotPowStoCap[s in S], dfStor[dfStor.r_id .==s, :existing_cap] .+ pow_gen[s, :rep_capacity]*(vNewPowStoCap[s] .- vRetPowStoCap[s]))
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
@expression(MP, eCostH2StoInv, sum(hsc_gen[s, :inv_cost_tonne_p_yr]*vNewH2StoCap[s]*hsc_gen[s, :rep_capacity] + hsc_gen[s, :inv_cost_comp_tonne_hr_p_yr]*vNewH2StoCompCap[s]
                                  + hsc_gen[s, :fom_cost_p_tonne_p_yr]*eTotH2StoCap[s] + hsc_gen[s, :fom_cost_comp_tonne_hr_p_yr]*eTotH2StoCompCap[s] for s in Q)
)
@expression(MP, eCostH2TraInv, sum(hsc_pipelines[i, :investment_cost_per_length]*hsc_pipelines[i, :distance]*vNewH2Pipe[i] +
                                    hsc_pipelines[i, :fom_per_length]*hsc_pipelines[i, :distance]*eTotH2Pipe[i] +
                                    hsc_pipelines[i, :compressor_inv_per_length]*hsc_pipelines[i, :distance]*vNewH2PipeCompCap[i] + 
                                    hsc_pipelines[i, :fom_comp_p_tonne_hr]*eTotH2PipeCompCap[i] for i in I)
)

# ---- MP Objective ---- # 
MP_obj = eCostPowGenInv .+ eCostPowStoInv .+ eCostPowTraInv .+ eCostH2GenInv .+ eCostH2StoInv .+ eCostH2TraInv 
@objective(MP, Min, MP_obj .+ etta)

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
@constraint(MP, cMaxH2PipeNum[i in I], eTotH2Pipe[i] <= hsc_pipelines[i, :max_num_pipes])
@constraint(MP, cMaxRetH2PipeNum[i in I], vRetH2Pipe[i] <= hsc_pipelines[i, :existing_num_pipes])
@constraint(MP, cMaxRetH2PipeCompCap[i in I], vRetH2PipeCompCap[i]<=hsc_pipelines[i, :existing_comp_cap_tonne_hr])
#Land Use Constraint on each zone
@constraint(MP, cLandUse[z in Z], ePowGenLandUse[z] + ePowStoLandUse[z] + eH2GenLandUse[z] + eH2StoLandUse[z] + eH2PipeLandUse[z] <= zones[z, :available_land])
#Policy
@constraint(MP, cZonalEmissionCap[z in Z], sum(vMaxEmissionByWeek[z,w] for w in W) <= 0.05*sum((h2_demand[w][t,z]*33.3) +pow_demand[w][t,z] for t in T, w in W))

#while (gap > tolerence*UB) && (iteration <= max_iter)
    
    iteration += 1

    optimize!(MP)
    if termination_status(MP) != MOI.OPTIMAL
        error("Master Problem infeasible. No solution Exisits")
    end
    LB = objective_value(MP)

##--------------------##
## Sub Problems model ##
##--------------------##

    SP = Model(Gurobi.Optimizer)

    set_optimizer_attribute(SP, "InfUnbdInfo", 1)
    set_optimizer_attribute(SP, "OutputFlag", 1)
    set_optimizer_attribute(SP, "DualReductions", 0)
    set_optimizer_attribute(SP, "NumericFocus", 3)

# ---- SP Variables ---- #
# Power Generation DV #
    @variable(SP, vPowGen[g in G, w in W, t in T]>=0)
    @variable(SP, vPowGenFirst[g in G, w in W]>=0)
    @variable(SP, vPowGenOnline[g in G_ther, w in W, t in T]>=0)
    @variable(SP, vPowGenOnlineFirst[g in G_ther, w in W]>=0)
    @variable(SP, vPowGenStart[g in G_ther, w in W, t in T]>=0)
    @variable(SP, vPowGenShut[g in G_ther, w in W, t in T]>=0)
# Power Storage DV 
    @variable(SP, vPowStoCha[s in S, w in W, t in T]>=0)
    @variable(SP, vPowStoDis[s in S, w in W, t in T]>=0)
    @variable(SP, vPowSOC[s in S, w in W, t in T]>=0)
    @variable(SP, vPowSOCFirst[s in S, w in W]>=0)
# Power Transmission DV #
    @variable(SP, vPowFlow[l in L, w in W, t in T])
# HSC Generation DV
    @variable(SP, vH2Gen[h in H, w in W, t in T]>=0)
    @variable(SP, vH2GenFirst[h in H, w in W]>=0)
    @variable(SP, vH2GenStart[h in H_ther,  w in W, t in T]>=0)
    @variable(SP, vH2GenShut[h in H_ther,  w in W, t in T]>=0)
    @variable(SP, vH2GenOnline[h in H_ther,  w in W, t in T]>=0)
    @variable(SP, vH2GenOnlineFirst[h in H_ther, w in W]>=0)
    # HSC Storage DV
    @variable(SP, vH2StoCha[s in Q, w in W, t in T]>=0)
    @variable(SP, vH2StoDis[s in Q, w in W, t in T]>=0)
    @variable(SP, vH2StoSOC[s in Q, w in W, t in T]>=0)
    @variable(SP, vH2StoSOCFirst[s in Q, w in W]>=0)
    # HSC Transmission DV
    @variable(SP, vH2Flow[i in I, w in W, t in T])
    # Policy Variables
    @variable(SP, vH2NSD[z in Z,  w in W, t in T]>=0)
    @variable(SP, vPowNSD[z in Z, w in W, t in T]>=0)
    @variable(SP, vExtraEmissionByZone[z in Z, w in W]>=0)
    @variable(SP, vPowSlack[z in Z, w in W, t in T]>=0)
    @variable(SP, vH2Slack[z in Z, w in W, t in T]>=0)
    # ---- SP Expressions ---- #
    # Power Generation Expressions #
    @expression(SP, eAvailPowGenCap[g in G], pow_gen[g, :existing_cap] .+ pow_gen[g,:rep_capacity]*(value(MP[:vNewPowGenCap][g]) - value(MP[:vRetPowGenCap][g])))
    @expression(SP, eAvailPowGenUnit[g in G_ther], pow_gen[g, :num_units]+ value(MP[:vNewPowGenCap][g])- value(MP[:vRetPowGenCap][g]))
    @expression(SP, ePowGenByZone[z in Z, w in W, t in T], sum(vPowGen[g,w,t] * (pow_gen[g, :zone] == z ? 1 : 0) for g in G))

    # Power Storage Expressions #
    @expression(SP, eAvailPowStoCap[s in S], dfStor[dfStor.r_id .==s, :existing_cap] .+ pow_gen[s, :rep_capacity]*(value(MP[:vNewPowStoCap][s]) .- value(MP[:vRetPowStoCap][s])))
    @expression(SP, ePowStoChaByZone[z in Z, w in W, t in T], sum(vPowStoCha[s,w,t] * (pow_gen[s, :zone] == z ? 1 : 0) for s in S))
    @expression(SP, ePowStoDisByZone[z in Z, w in W, t in T], sum(vPowStoDis[s,w,t] * (pow_gen[s, :zone] == z ? 1 : 0) for s in S))
    # Power Transmission Expressions #
    @expression(SP, eAvailPowTraCap[l in L], pow_lines[l, :existing_transmission_cap_mw] .+ value(MP[:vNewPowTraCap][l])) #No retired cap considered for transmission lines + Not var cost for power flows - Rep Cap is cosidered 1
    @expression(SP, ePowFlowbyZone[l in L, z in Z, w in W, t in T], vPowFlow[l,w,t]*Pow_Network[l,z])
    @expression(SP, eNet_Pow_Flow[z in Z, w in W, t in T], sum(Pow_Network[l,z] * vPowFlow[l,w,t] for l in L))
    @expression(SP, ePow_Loss_By_Zone[z in Z, w in W, t in T], sum(abs(Pow_Network[l,z]) * (1/2) *vPowFlow[l,w,t] * pow_lines[l, :line_loss_percentage] for l in L))
    # Power Cost Expressions
    @expression(SP, eCostPowGenVar, sum((pow_gen[g, :vom_cost_mwh] + pow_gen[g, :heat_rate_mmbtu_per_yr] .* fuel_costs[pow_gen[g, :fuel]][t]) .* vPowGen[g,w,t] for g in G, w in W, t in T))
    @expression(SP, eCostPowGenStart, sum(pow_gen[g, :start_cost_per_mw] .* pow_gen[g, :rep_capacity] .* vPowGenStart[g,w,t] for g in G_ther, w in W, t in T)) 
    #For cost of Non-served demand we only consider $/MWh for each zone and will not consider demand segments
    @expression(SP, eCostPowNSD, sum(vPowNSD[z,w,t] .* zones[z, :voll_pow] for z in Z, w in W, t in T)) 
    @expression(SP, eCostPowSlack, sum(vPowSlack[z,w,t] .* big_M for z in Z, w in W, t in T))
    @expression(SP, eCostPowStoVar, sum(vPowStoCha[s,w,t] .* pow_gen[s, :vom_cost_mwh_charge] for s in S, w in W, t in T))

    # HSC Generation Expressions
    @expression(SP, eAvailH2GenCap[h in H], hsc_gen[h, :existing_cap_tonne_p_hr] .+ hsc_gen[h, :rep_capacity]*(value(MP[:vNewH2GenCap][h]) .- value(MP[:vRetH2GenCap][h])))
    @expression(SP, eAvailH2GenUnit[h in H_ther], hsc_gen[h, :num_units]+ value(MP[:vNewH2GenCap][h])- value(MP[:vRetH2GenCap][h]))
    @expression(SP, eH2GenEvap[h in H, w in W, t in T], hsc_gen[h, :boil_off]*vH2Gen[h,w,t])
    # HSC Storage Expression
    @expression(SP, eAvailH2StoCap[s in Q], hsc_gen[s, :existing_cap_tonne] + hsc_gen[s, :rep_capacity]*(value(MP[:vNewH2StoCap][s]) - value(MP[:vRetH2StoCap][s])))
    @expression(SP, eAvailH2StoCompCap[s in Q], hsc_gen[s, :existing_cap_comp_tonne_hr] + value(MP[:vNewH2StoCompCap][s]) - value(MP[:vRetH2StoCompCap][s])) #rep cap is considered 1
    # HSC Tramsmission Expression
    @expression(SP, eAvailH2Pipe[i in I], hsc_pipelines[i, :existing_num_pipes] + value(MP[:vNewH2Pipe][i]) - value(MP[:vRetH2Pipe][i]))
    @expression(SP, eH2FolwbyZone[i in I, z in Z, w in W, t in T], vH2Flow[i,w,t]*H2_Network[i,z])
    @expression(SP, eNet_H2_Flow[z in Z,w in W, t in T], sum(H2_Network[i,z] * vH2Flow[i,w,t] *(-1) for i in I))
    @expression(SP, eH2_Loss_By_Zone[z in Z, w in W, t in T], sum(abs(H2_Network[i,z]) * (1/2) *vH2Flow[i,w,t] * hsc_pipelines[i, :pipe_loss_coeff] for i in I ))
    @expression(SP, eAvailH2PipeCompCap[i in I], hsc_pipelines[i, :existing_comp_cap_tonne_hr] + value(MP[:vNewH2PipeCompCap][i]) - value(MP[:vRetH2PipeCompCap][i]))
    # HSC Cost Expressions
    @expression(SP, eCostH2GenVar, sum((hsc_gen[h, :vom_cost_p_tonne] + hsc_gen[h, :heat_rate_mmbtu_p_tonne] .* fuel_costs[lowercase(hsc_gen[h, :fuel])][t]) .* vH2Gen[h,w,t] for h in H, w in W, t in T)    )
    @expression(SP, eCostH2StoVar, sum(hsc_gen[s, :var_om_cost_charge_p_tonne]*vH2StoCha[s,w,t] for s in Q, w in W, t in T))
    @expression(SP, eCostH2GenStart, sum(hsc_gen[h, :startup_cost_p_tonne_hr] .* hsc_gen[h, :rep_capacity] .* vH2GenStart[h,w,t] for h in H_ther, w in W, t in T)) 

    @expression(SP, eCostH2NSD, sum(vH2NSD[z,w,t] .* zones[z, :voll_hsc] for z in Z, w in W,t in T))
    @expression(SP, eCostH2Slack, sum(vH2Slack[z,w,t] * big_M for z in Z, w in W, t in T))

    @expression(SP, ePowDemandHSC[w in W, t in T, z in Z], sum(hsc_gen[h, :pow_demand_mwh_p_tonne]*vH2Gen[h,w,t]*(hsc_gen[h, :zone]==z ? 1 : 0) for h in H) +
                                                        sum(hsc_gen[s, :h2charge_mwh_p_tonne]*vH2StoCha[s,w,t]*(hsc_gen[s, :zone]==z ? 1 : 0) for s in Q) +
                                                        sum(hsc_pipelines[i, :comp_pow_mwh_per_tonne]*vH2Flow[i,w,t]*(H2_Network[i,z]==1 ? 1 : 0) for i in I)
    )
    @expression(SP, eH2DemandPow[w in W, t in T, z in Z], sum(pow_gen[g, :h2_demand_tonne_p_mwh]*vPowGen[g,w,t]*(pow_gen[g, :zone]==z ? 1 : 0) for g in G))
    @expression(SP, pow_D[w in W, t in T, z in Z], pow_demand[w][t,z] .+ ePowDemandHSC[w,t,z])
    @expression(SP, H2_D[w in W, t in T, z in Z], h2_demand[w][t,z] .+ eH2DemandPow[w,t,z])
    @expression(SP, eZonalEmissionCapByWeek[z in Z, w in W], sum(vPowGen[g,w,t]*pow_gen[g, :heat_rate_mmbtu_per_yr]*CO2_content[pow_gen[g, :fuel]]*(pow_gen[g, :zone] == z ? 1 : 0) for g in G, t in T) + 
                                                              sum(vH2Gen[h,w,t]*hsc_gen[h, :heat_rate_mmbtu_p_tonne]*CO2_content[hsc_gen[h, :fuel]] *(hsc_gen[h, :zone] == z ? 1 : 0) for h in H, t in T)
    )
    #Emission Cost Expression
    @expression(SP, eEmissionCost, sum(sum(vExtraEmissionByZone[z,w] for w in W)*zones[z, :emission_cost] for z in Z))

    # ---- SP Objective ---- #  
    SP_obj = eCostPowGenVar .+ eCostPowStoVar .+ eCostPowNSD .+ eCostPowSlack .+ eCostPowGenStart .+ eCostH2GenVar .+ eCostH2StoVar .+ eCostH2NSD .+ eCostH2Slack .+ eCostH2GenStart .+ eEmissionCost
    @objective(SP, Min, SP_obj)
    # ---- SP Constraints ---- #
    # Power Balance #
    @constraint(SP, cPowerBalance[z in Z, w in W,t in T],
        ePowGenByZone[z,w,t] .- eNet_Pow_Flow[z,w,t] .- ePow_Loss_By_Zone[z,w,t] .+ ePowStoDisByZone[z,w,t] .- ePowStoChaByZone[z,w,t] .+ vPowNSD[z,w,t] + vPowSlack[z,w,t] >= pow_D[w,t,z]
    )
    # Power Generation #
    @constraint(SP, cMaxPowGen[g in G_ren, w in W, t in T], vPowGen[g,w,t] - eAvailPowGenCap[g] * pow_gen_var[H_w[w][t], pow_gen[g, :resource]] <= 0)
    @constraint(SP, cMaxPowGenTher[g in G_ther, w in W, t in T], vPowGen[g,w,t] .- (pow_gen[g, :rep_capacity].*pow_gen[g,:max_op_level].*vPowGenOnline[g,w,t]) <= 0)
    @constraint(SP, cMinPowGenTher[g in G_ther,w in W,t in T], (pow_gen[g, :rep_capacity]*pow_gen[g,:min_op_level]*vPowGenOnline[g,w,t]) .- vPowGen[g,w,t] <= 0)
    @constraint(SP, cPowOnlineUnits[g in G_ther, w in W, t in T], vPowGenOnline[g,w,t] .- eAvailPowGenUnit[g]<= 0)
    @constraint(SP, cPowStartLimits[g in G_ther,w in W,t in T], vPowGenStart[g,w,t] .- eAvailPowGenUnit[g] .+ vPowGenOnline[g,w,t]<= 0)
    @constraint(SP, cPowShutLimits[g in G_ther,w in W,t in T], vPowGenShut[g,w,t] .- vPowGenOnline[g,w,t] <= 0)
    #Cyclic constraint on Thermal power units
    @constraint(SP, cPowUnitOnlineCon[g in G_ther, w in W, t in 2:length(T)], vPowGenOnline[g,w,t] - vPowGenOnline[g,w,t-1] == vPowGenStart[g,w,t]-vPowGenShut[g,w,t])
    @constraint(SP, cPowUnitOnlineFirst[g in G_ther, w in W], vPowGenOnline[g,w,1] - vPowGenOnlineFirst[g,w] == vPowGenStart[g,w,1]-vPowGenShut[g,w,1])
    @constraint(SP, cPowUnitOnlineCycle[g in G_ther, w in W], vPowGenOnlineFirst[g,w] == vPowGenOnline[g,w,168])
    #Ramping of non-thermal units 
    @constraint(SP, cPowGenRampUp[g in G_ren, w in W, t in 2:length(T)], vPowGen[g,w,t]-vPowGen[g,w,t-1] .- pow_gen[g, :ramp_up] * eAvailPowGenCap[g]<= 0) 
    @constraint(SP, cPowGenRampUpFirst[g in G_ren, w in W], vPowGen[g,w,1]-vPowGenFirst[g,w] .- pow_gen[g, :ramp_up]*eAvailPowGenCap[g]<=0)
    @constraint(SP, cPowGenRampDn[g in G_ren, w in W, t in 2:length(T)], vPowGen[g,w,t-1]-vPowGen[g,w,t] .- pow_gen[g, :ramp_dn] * eAvailPowGenCap[g]<= 0)
    @constraint(SP, cPowGenRampDnFirst[g in G_ren, w in W], vPowGenFirst[g,w]-vPowGen[g,w,1] .- pow_gen[g, :ramp_dn] * eAvailPowGenCap[g]<= 0)
    @constraint(SP, cPowGenCycle[g in G_ren, w in W], vPowGenFirst[g,w] == vPowGen[g,w,168])
    #Ramping of Thermal units
    @constraint(SP, cTherPowGenRampDn[g in G_ther, w in W, t in 2:length(T)], 
                vPowGen[g,w,t-1] .- vPowGen[g,w,t] .- pow_gen[g, :rep_capacity]*pow_gen[g,:ramp_dn] *(vPowGenOnline[g,w,t].-vPowGenStart[g,w,t])
                .+ pow_gen[g, :rep_capacity]*vPowGenStart[g,w,t]*pow_gen[g,:min_op_level] .- vPowGenShut[g,w,t]*pow_gen[g, :rep_capacity]*min(pow_gen_var[H_w[w][t], pow_gen[g, :resource]],max(pow_gen[g,:min_op_level],pow_gen[g,:ramp_dn]))<=0
    )
    @constraint(SP, cTherPowGenRampDnFirst[g in G_ther, w in W], 
                vPowGenFirst[g,w] .- vPowGen[g,w,1] .- pow_gen[g, :rep_capacity]*pow_gen[g,:ramp_dn] *(vPowGenOnline[g,w,1].-vPowGenStart[g,w,1])
                .+ pow_gen[g, :rep_capacity]*vPowGenStart[g,w,1]*pow_gen[g,:min_op_level] .- vPowGenShut[g,w,1]*pow_gen[g, :rep_capacity]*min(pow_gen_var[H_w[w][1], pow_gen[g, :resource]],max(pow_gen[g,:min_op_level],pow_gen[g,:ramp_dn]))<=0
    )

    @constraint(SP, cTherPowGenRampUp[g in G_ther, w in W, t in 2:length(T)], 
                vPowGen[g,w,t] .- vPowGen[g,w,t-1] .- pow_gen[g, :rep_capacity]*pow_gen[g,:ramp_up] *(vPowGenOnline[g,w,t].-vPowGenStart[g,w,t])
                .+ pow_gen[g, :rep_capacity]*vPowGenShut[g,w,t]*pow_gen[g,:min_op_level] .- vPowGenStart[g,w,t]*pow_gen[g, :rep_capacity]*min(pow_gen_var[H_w[w][t], pow_gen[g, :resource]],max(pow_gen[g,:min_op_level],pow_gen[g,:ramp_up]))<=0
    )   

    @constraint(SP, cTherPowGenRampUpFirst[g in G_ther, w in W], 
                vPowGen[g,w,1] .- vPowGenFirst[g,w] .- pow_gen[g, :rep_capacity]*pow_gen[g,:ramp_up] *(vPowGenOnline[g,w,1].-vPowGenStart[g,w,1])
                .- pow_gen[g, :rep_capacity]*vPowGenShut[g,w,1]*pow_gen[g,:min_op_level] .+ vPowGenStart[g,w,1]*pow_gen[g, :rep_capacity]*min(pow_gen_var[H_w[w][1], pow_gen[g, :resource]],max(pow_gen[g,:min_op_level],pow_gen[g,:ramp_up]))<=0
    )
    #Minimum up/donw time for thermal generators
    @constraint(SP, cMinUpTimePowGen[g in G_ther,w in W,t in T], sum(vPowGenStart[g,w,tt] for tt in intersect(T, (t - pow_gen[g, :up_time]):t)) .- vPowGenOnline[g,w,t]<= 0)
    @constraint(SP, cMinDnTimePowGen[g in G_ther,w in W,t in T], sum(vPowGenShut[g,w,tt] for tt in intersect(T, (t - pow_gen[g, :down_time]):t)) .- eAvailPowGenUnit[g] .+vPowGenOnline[g,w,t]<= 0)
    # Power Storage #
    @constraint(SP, cPowStoBalance[s in S, w in W, t in 2:length(T)], vPowSOC[s,w,t] == (1-pow_gen[s,:etta_self_dis])*vPowSOC[s,w,t-1] + pow_gen[s,:etta_cha]*vPowStoCha[s,w,t] - (1/pow_gen[s,:etta_dis])*vPowStoDis[s,w,t])
    @constraint(SP, cPowStoBalanceFirst[s in S, w in W], vPowSOC[s,w,1] == (1-pow_gen[s,:etta_self_dis])*vPowSOCFirst[s,w] + pow_gen[s,:etta_cha]*vPowStoCha[s,w,1] .- (1/pow_gen[s,:etta_dis])*vPowStoDis[s,w,1])
    @constraint(SP, cPowStoMaxDis[s in S, w in W, t in 2:length(T)], vPowStoDis[s,w,t] .- pow_gen[s,:etta_dis]*vPowSOC[s,w,t-1]<= 0)
    @constraint(SP, cPowStoMaxDisFirst[s in S, w in W], vPowStoDis[s,w,1] .- pow_gen[s,:etta_dis]*vPowSOCFirst[s,w]<= 0)
    @constraint(SP, cPowStoMaxCha[s in S, w in W, t in T], vPowStoCha[s,w,t] .- eAvailPowStoCap[s]<=0)
    @constraint(SP, cPowSOCCycle[s in S, w in W], vPowSOCFirst[s,w] == vPowSOC[s,w,168])
    @constraint(SP, cPowStoSOCMax[s in S, w in W, t in T], vPowSOC[s,w,t] .- eAvailPowStoCap[s]*pow_gen[s, :max_op_level]<=0)
    @constraint(SP, cPowStoSOCMin[s in S, w in W, t in T], eAvailPowStoCap[s]*pow_gen[s, :min_op_level] .- vPowSOC[s,w,t]<=0)
    # Power Transmission #
    @constraints(SP, begin
                    cMaxPowFlowOut[l in L, w in W, t in T],  vPowFlow[l,w,t] .- eAvailPowTraCap[l]<=0 
                    cMaxPowFlowIn[l in L, w in W, t in T],  -eAvailPowTraCap[l] .- vPowFlow[l,w,t] <=0
    end)

    # H2 Balance constraint
    @constraint(SP, cH2Balance[z in Z, w in W, t in T],
      sum((vH2Gen[h,w,t] .- eH2GenEvap[h,w,t])*(hsc_gen[h, :zone]==z ? 1 : 0) for h in H) .- sum(H2_Network[i,z]*vH2Flow[i,w,t] for i in I) .- 0.5*sum(abs(H2_Network[i,z])*vH2Flow[i,w,t]*hsc_pipelines[i, :pipe_loss_coeff] for i in I) .+
      sum((vH2StoDis[s,w,t] - vH2StoCha[s,w,t])*(hsc_gen[s, :zone]==z ? 1 : 0) for s in Q) + vH2NSD[z,w,t] + vH2Slack[z,w,t] >= H2_D[w,t,z]
    )
    # H2 Generation constraint
    @constraint(SP, cMaxH2GenVar[h in H_dis, w in W, t in T], vH2Gen[h,w,t] .- eAvailH2GenCap[h]*hsc_gen[h, :max_op_level]*hsc_gen_var[H_w[w][t], hsc_gen[h, :resource]] <= 0)
    @constraint(SP, cMinH2GenVar[h in H_dis, w in W, t in T], eAvailH2GenCap[h]*hsc_gen[h,:min_op_level]*hsc_gen_var[H_w[w][t], hsc_gen[h, :resource]]  .- vH2Gen[h,w,t] <= 0)
    # H2 Thermal units constraints
    @constraint(SP, cMaxH2GenTher[h in H_ther, w in W, t in T], vH2Gen[h,w,t] .- hsc_gen[h, :max_op_level] * hsc_gen[h, :rep_capacity] * vH2GenOnline[h,w,t]*hsc_gen_var[H_w[w][t], hsc_gen[h, :resource]]<= 0)
    @constraint(SP, cMinH2GenTher[h in H_ther, w in W, t in T], hsc_gen[h, :min_op_level] * hsc_gen[h, :rep_capacity] * vH2GenOnline[h,w,t]*hsc_gen_var[H_w[w][t], hsc_gen[h, :resource]] .- vH2Gen[h,w,t]<= 0)
    @constraint(SP, cH2OnlineUnits[h in H_ther, w in W, t in T], vH2GenOnline[h,w,t] .- eAvailH2GenUnit[h]<= 0)
    @constraint(SP, cH2StartLimits[h in H_ther, w in W, t in T], vH2GenStart[h,w,t] .- eAvailH2GenUnit[h] .+ vH2GenOnline[h,w,t]<= 0 )
    @constraint(SP, cH2ShutLimits[h in H_ther, w in W, t in T], vH2GenShut[h,w,t] .- vH2GenOnline[h,w,t]<= 0)   
# H2 thermal units cyclic constraints
    @constraint(SP, cH2UnitOnlineCon[h in H_ther, w in W, t in 2:length(T)], vH2GenOnline[h,w,t] .- vH2GenOnline[h,w,t-1] == vH2GenStart[h,w,t].- vH2GenShut[h,w,t])
    @constraint(SP, cH2UnitOnlineConFirst[h in H_ther, w in W], vH2GenOnline[h,w,1] .- vH2GenOnlineFirst[h,w] == vH2GenStart[h,w,1] .- vH2GenShut[h,w,1])
    @constraint(SP, cH2GenOnlineConCycle[h in H_ther, w in W], vH2GenOnlineFirst[h,w] == vH2GenOnline[h,w,168])
    # Min Up and Down time for Thermal H2 generators
    @constraint(SP, cMinUpTimeH2Gen[h in H_ther, w in W, t in T], sum(vH2GenStart[h,w,tt] for tt in intersect(T, (t - hsc_gen[h, :up_time]):t)) .- vH2GenOnline[h,w,t]<= 0)
    @constraint(SP, cMinDnTimeH2Gen[h in H_ther, w in W, t in T], sum(vH2GenShut[h,w,tt] for tt in intersect(T, (t - hsc_gen[h, :down_time]):t)) .- eAvailH2GenUnit[h] .+ vH2GenOnline[h,w,t]<= 0)
    # Ramp constraints for dispatachable units
    @constraint(SP, cH2GenRampUp[g in H_dis, w in W, t in 2:length(T)], vH2Gen[g,w,t]-vH2Gen[g,w,t-1] .- hsc_gen[g, :ramp_up_percentage] * eAvailH2GenCap[g]<= 0) 
    @constraint(SP, cH2GenRampUpFirst[g in H_dis, w in W], vH2Gen[g,w,1]-vH2GenFirst[g,w] .- hsc_gen[g, :ramp_up_percentage]*eAvailH2GenCap[g]<=0)
    @constraint(SP, cH2GenRampDn[g in H_dis, w in W, t in 2:length(T)], vH2Gen[g,w,t-1]-vH2Gen[g,w,t] .- hsc_gen[g, :ramp_down_percentage] * eAvailPowGenCap[g]<= 0)
    @constraint(SP, cH2GenRampDnFirst[g in H_dis, w in W], vH2GenFirst[g,w]-vH2Gen[g,w,1] .- hsc_gen[g, :ramp_down_percentage] * eAvailPowGenCap[g]<= 0)
    @constraint(SP, cH2GenCycle[g in H_dis, w in W], vH2GenFirst[g,w] == vH2Gen[g,w,168])
    #Ramp Constraint for Thermal units
    @constraint(SP, cTherH2GenRampDn[h in H_ther, w in W, t in 2:length(T)], 
            vH2Gen[h,w,t-1] .- vH2Gen[h,w,t] .- hsc_gen[h, :rep_capacity]*hsc_gen[h,:ramp_down_percentage] *(vH2GenOnline[h,w,t].-vH2GenStart[h,w,t])
            .+ hsc_gen[h, :rep_capacity]*vH2GenStart[h,w,t]*hsc_gen[h,:min_op_level] .- vH2GenShut[h,w,t]*hsc_gen[h, :rep_capacity]*min(hsc_gen_var[H_w[w][t], hsc_gen[h, :resource]],max(hsc_gen[h,:min_op_level],hsc_gen[h,:ramp_down_percentage]))<=0
    )
    @constraint(SP, cTherH2GenRampDnFirst[h in H_ther, w in W], 
                vH2GenFirst[h,w] .- vH2Gen[h,w,1] .- hsc_gen[h, :rep_capacity]*hsc_gen[h,:ramp_down_percentage] *(vH2GenOnline[h,w,1].-vH2GenStart[h,w,1])
                .+ hsc_gen[h, :rep_capacity]*vH2GenStart[h,w,1]*hsc_gen[h,:min_op_level] .- vH2GenShut[h,w,1]*hsc_gen[h, :rep_capacity]*min(hsc_gen_var[H_w[w][1], hsc_gen[h, :resource]],max(hsc_gen[h,:min_op_level],hsc_gen[h,:ramp_down_percentage]))<=0
    )

    @constraint(SP, cTherH2GenRampUp[h in H_ther, w in W, t in 2:length(T)], 
                vH2Gen[h,w,t] .- vH2Gen[h,w,t-1] .- hsc_gen[h, :rep_capacity]*hsc_gen[h,:ramp_up_percentage] *(vH2GenOnline[h,w,t].-vH2GenStart[h,w,t])
                .+ hsc_gen[h, :rep_capacity]*vH2GenShut[h,w,t]*hsc_gen[h,:min_op_level] .- vH2GenStart[h,w,t]*hsc_gen[h, :rep_capacity]*min(hsc_gen_var[H_w[w][t], hsc_gen[h, :resource]],max(hsc_gen[h,:min_op_level],hsc_gen[h,:ramp_up_percentage]))<=0
    )

    @constraint(SP, cTherH2GenRampUpFirst[g in H_ther, w in W], 
                vH2Gen[g,w,1] .- vH2GenFirst[g,w] .- hsc_gen[g, :rep_capacity]*hsc_gen[g,:ramp_up_percentage] *(vH2GenOnline[g,w,1].-vH2GenStart[g,w,1])
                .- hsc_gen[g, :rep_capacity]*vH2GenShut[g,w,1]*hsc_gen[g,:min_op_level] .+ vH2GenStart[g,w,1]*hsc_gen[g, :rep_capacity]*min(hsc_gen_var[H_w[w][1], hsc_gen[g, :resource]],max(hsc_gen[g,:min_op_level],hsc_gen[g,:ramp_up_percentage]))<=0
    )
    # H2 Storage constraint
    @constraint(SP, cH2StoBalance[s in Q, w in W, t in 2:length(T)], vH2StoSOC[s,w,t] == (1-hsc_gen[s, :etta_self_dis])*vH2StoSOC[s,w,t-1] + vH2StoCha[s,w,t]*hsc_gen[s,:etta_cha] - (1/hsc_gen[s,:etta_dis])*vH2StoDis[s,w,t])
    @constraint(SP, cH2StoBalanceFirst[s in Q, w in W], vH2StoSOC[s,w,1] == (1-hsc_gen[s, :etta_self_dis])*vH2StoSOCFirst[s,w] + vH2StoCha[s,w,1]*hsc_gen[s,:etta_cha] - (1/hsc_gen[s,:etta_dis])*vH2StoDis[s,w,1])
    @constraint(SP, cMaxH2StoSOC[s in Q, w in W, t in T], vH2StoSOC[s,w,t] .- hsc_gen[s,:h2stor_max_level]*eAvailH2StoCap[s]<= 0)
    @constraint(SP, cMinH2StoSOC[s in Q, w in W, t in T], hsc_gen[s,:h2stor_min_level]*eAvailH2StoCap[s] .- vH2StoSOC[s,w,t]<= 0 )
    #@constraint(SP, cMaxH2StoChar[s in Q,w in W,t in T], vH2StoCha[s,w,t] .- eAvailH2StoCompCap[s] <= 0)
    # H2 Transmission constraints
    @constraints(SP, begin 
                      cMaxH2PipeFlowOut[i in I, w in W, t in T], vH2Flow[i,w,t] .- eAvailH2Pipe[i]*hsc_pipelines[i, :max_op_level] <= 0
                      cMaxH2PipeFlowIn[i in I, w in W, t in T], - vH2Flow[i,w,t] .- eAvailH2Pipe[i]*hsc_pipelines[i, :max_op_level] <= 0
    end)
    @constraints(SP, begin
                      cMaxH2PipeFlowOutComp[i in I,w in W,t in T], vH2Flow[i,w,t] .- eAvailH2PipeCompCap[i] <= 0
                      cMaxH2PipeFlowInComp[i in I, w in W, t in T], -vH2Flow[i,w,t] .- eAvailH2PipeCompCap[i] <= 0
    end)

    # Policy constraints
    @constraint(SP, cPowNSD[z in Z, w in W, t in T], vPowNSD[z,w,t] .- zones[z, :pow_nsd_share]*pow_D[w,t,z] <= 0 )
    @constraint(SP, cH2NSD[z in Z, w in W, t in T], vH2NSD[z,w,t] .- zones[z, :hsc_nsd_share]*H2_D[w, t,z] <= 0)

    #Emission constraint
    @constraint(SP, cZonalEmissionCapByWeek[z in Z, w in W], eZonalEmissionCapByWeek[z,w] .- vExtraEmissionByZone[z,w] .- value(MP[:vMaxEmissionByWeek][z,w])<= 0)    

    optimize!(SP)

#=

    if termination_status(SP) == MOI.OPTIMAL
        UB_candidate = objective_value(SP) + objective_value(MP)
        UB = min(UB, UB_candidate)
        
        pi_cMaxPowGen[g in G_ren] = sum(dual(cMaxPowGen[g,w,t])* pow_gen[g, :rep_capacity]*pow_gen_var[H_w[w][t], pow_gen[g, :resource]] for w in W, t in T)
        pi_cPowOnlineUnits[g in G_ther] = sum((dual(cPowOnlineUnits[g,w,t])+dual(cPowStartLimits[g,w,t]))*pow_gen[g, :rep_capacity] for w in W, t in T)
        pi_cPowGenRamp[g in G_ren] = sum((dual(cPowGenRampUp[g,w,t])+dual(cPowGenRampDn[g,w,t])) * pow_gen[g, :ramp_up]*pow_gen[g, :rep_capacity] for w in W, t in 2:length(T))
        pi_cPowStoSOC[s in S] = sum(dual(cPowStoSOCMax[s,w,t])*pow_gen[s, :rep_capacity]*pow_gen[s, :max_op_level] - dual(cPowStoSOCMin)*pow_gen[s, :rep_capacity]*pow_gen[s, :min_op_level] for w in W, t in 2:length(T))
        pi_cMaxPowFlow[l in L] = sum(dual(cMaxPowFlowOut[l,w,t])-dual(cMaxPowFlowIn[l,w,t]) for w in W, t in T)

        pi_cMaxH2GenVar[h in H_dis] = sum((dual(cMaxH2GenVar[h,w,t])-dual(cMinH2GenVar[h,w,t]))*hsc_gen[h, :max_op_level]*hsc_gen_var[H_w[w][t], hsc_gen[h, :resource]]*hsc_gen[h, :rep_capacity] for w in W, t in T)
        pi_cH2OnlineUnits[h in H_ther] = sum((dual(cH2OnlineUnits[h,w,t])+dual(cH2StartLimits[h,w,t]))*hsc_gen[h, :rep_capacity] for w in W, w in T)
        pi_cH2GenRamp[h in H_dis] = sum((dual(cH2GenRampUp[h,w,t])*hsc_gen[h, :ramp_up_percentage] + dual(cH2GenRampDn[h,w,t]*hsc_gen[h, :ramp_down_percentage]))*hsc_gen[h, :rep_capacity] for w in W, t in 2:length(T))
        pi_cH2StoSOC[s in Q] = sum((dual(cMaxH2StoSOC[s,w,t])*hsc_gen[s,:h2stor_max_level] - dual(cMinH2StoSOC[s,w,t])*hsc_gen[s,:h2stor_min_level])*hsc_gen[s, :rep_capacity] for w in W, t in 2:length(T))
        pi_cMaxH2StoChar[s in Q] = sum(dual(cMaxH2StoChar[s,w,t]) for w in W, t in T)
        pi_cH2PipeFlow[i in I] = sum((dual(cMaxH2PipeFlowOut[i,w,t])-dual(cMaxH2PipeFlowIn[i,w,t]))*hsc_pipelines[i, :max_op_level] for w in W, t in T)
        pi_cH2PipeFlowComp[i in I] = sum(dual(cMaxH2PipeFlowOutComp[i,w,t]) - dual(cMaxH2PipeFlowIn[i,w,t]) for w in W, i in I)
        
        @constraint(MP, BD_Opt_Cut,
            etta >= objective_value(SP)+ 
            sum(pi_cMaxPowGen[g]*(vNewPowGenCap[g] - vRetPowGenCap[g] - value(MP[:vNewPowGenCap][g]) + value(MP[:vRetPowGenCap][g])) for g in G_ren) +
            sum(pi_cPowOnlineUnits[g]*(vNewPowGenCap[g]-vRetPowGenCap[g] - value(MP[:vNewPowGenCap][g]) + value(MP[:vRetPowGenCap][g])) for g in G_ther) +
            sum(pi_cPowGenRamp[g]*(vNewPowGenCap[g]-vRetPowGenCap[g] - value(MP[:vNewPowGenCap][g]) + value(MP[:vRetPowGenCap][g])) for g in G_ren) +
            sum(pi_cPowStoSOC[s]*(vNewPowStoCap[s]-vRetPowStoCap[s] - value(MP[:vNewPowGenCap][s]) + values(MP[:vRetPowStoCap][s])) for s in S) +
            sum(pi_cMaxPowFlow[l]*(vNewPowTraCap[l] - value(MP[:vNewPowTraCap][l]) for l in L)) +
            
            sum(pi_cMaxH2GenVar[h]*(vNewH2GenCap[h]-vRetH2GenCap[h] - value(MP[:vNewH2GenCap][h]) + value(MP[:vRetH2GenCap][h])) for h in H_dis) +
            sum(pi_cH2OnlineUnits[h]*(vNewH2GenCap[h]-vRetH2GenCap[h] - value(MP[:vNewH2GenCap][h]) + value(MP[:vRetH2GenCap][h])) for h in H_ther) +
            sum(pi_cH2GenRamp[h]*(vNewH2GenCap[h]-vRetH2GenCap[h] - value(MP[:vNewH2GenCap][h]) + value(MP[:vRetH2GenCap][h])) for h in H_dis) +
            sum(pi_cH2StoSOC[s]*(vNewH2StoCap[s]-vRetH2StoCap[s] - value(MP[:vNewH2GenCap][s]) + values(MP[:vRetH2StoCap][s])) for s in Q) +
            sum(pi_cMaxH2StoChar[s]*(vNewH2StoCompCap[s]-vRetH2StoCompCap[s] - value(MP[:vNewH2StoCompCap][s]) + value(MP[:vRetH2StoCompCap][s])) for s in Q) +
            sum(pi_cH2PipeFlow[i]*(vNewH2Pipe[i]-vRetH2Pipe[i] - value(MP[:vNewH2Pipe][i]) + value(MP[:vRetH2Pipe][i])) for i in I) +
            sum(pi_cH2PipeFlowComp[i]*(vNewH2PipeCompCap[i]-vRetH2PipeCompCap[i] - value(MP[:vNewH2PipeCompCap][i]) + value(MP[:vRetH2PipeCompCap][i])) for i in I)
        )
        
    else termination_status(SP) == MOI.INFEASIBLE
        
compute_conflict!(SP)

conflict_constraints = []

        for (func_type, set_type) in list_of_constraint_types(SP)
            for con in all_constraints(SP, func_type, set_type)
                if MOI.get(SP, MOI.ConstraintConflictStatus(), con) == MOI.IN_CONFLICT
                    push!(conflict_constraints, con)
                    println("Conflict identified: ", con)
                end
            end
        end
        BD_Feas_Cut = 0
        if !isempty(conflict_constraints)
            if cMaxPowGen in conflict_constraints
                BD_Feas_Cut += sum(vNewPowGenCap[g] - vRetPowGenCap[g] - value(MP[:vNewPowGenCap][g]) + value(MP[:vRetPowGenCap][g]) for g in G_ren, w in W, t in T)
            end
            if cPowOnlineUnits in conflict_constraints
                BD_Feas_Cut += sum(vNewPowGenCap[g]-vRetPowGenCap[g] - value(MP[:vNewPowGenCap][g]) + value(MP[:vRetPowGenCap][g]) for g in G_ther, w in W, t in T)
            end
            if cPowStartLimits in conflict_constraints
                BD_Feas_Cut += sum(vNewPowGenCap[g]-vRetPowGenCap[g] - value(MP[:vNewPowGenCap][g]) + value(MP[:vRetPowGenCap][g]) for g in G_ther, w in W, t in T) 
            end
            if cPowGenRampUp in conflict_constraints
                BD_Feas_Cut += sum(vNewPowGenCap[g]-vRetPowGenCap[g] - value(MP[:vNewPowGenCap][g]) + value(MP[:vRetPowGenCap][g]) for g in G_ren, w in W, t in 2:length(T)) 
            end
            if cPowGenRampDn in conflict_constraints
                BD_Feas_Cut += sum(vNewPowGenCap[g]-vRetPowGenCap[g] - value(MP[:vNewPowGenCap][g]) + value(MP[:vRetPowGenCap][g]) for g in G_ren, w in W, t in 2:length(T))
            end
            if cPowStoSOCMax in conflict_constraints
                BD_Feas_Cut += sum(vNewPowStoCap[s]-vRetPowStoCap[s] - value(MP[:vNewPowGenCap][s]) + values(MP[:vRetPowStoCap][s]) for s in S, w in W, t in 2:length(T)) 
            end
            if cPowStoSOCMin in conflict_constraints
                BD_Feas_Cut += sum(vNewPowStoCap[s]-vRetPowStoCap[s] - value(MP[:vNewPowGenCap][s]) + values(MP[:vRetPowStoCap][s]) for s in S, w in W, t in 2:length(T)) 
            end
            if cMaxPowFlowOut in conflict_constraints
                BD_Feas_Cut +=sum(vNewPowTraCap[l] - value(MP[:vNewPowTraCap][l]) for l in L, w in W, t in T)
            end
            if cMaxPowFlowIn in conflict_constraints
                BD_Feas_Cut +=sum(vNewPowTraCap[l] - value(MP[:vNewPowTraCap][l]) for l in L, w in W, t in T) 
            end
            if cMaxH2GenVar in conflict_constraints
                BD_Feas_Cut +=sum(vNewH2GenCap[h]-vRetH2GenCap[h] - value(MP[:vNewH2GenCap][h]) + value(MP[:vRetH2GenCap][h]) for h in H_dis, w in W, t in T) 
            end
            if cMinH2GenVar in conflict_constraints
                BD_Feas_Cut +=sum(vNewH2GenCap[h]-vRetH2GenCap[h] - value(MP[:vNewH2GenCap][h]) + value(MP[:vRetH2GenCap][h]) for h in H_dis, w in W, t in T)
            end
            if cH2OnlineUnits in conflict_constraints
                BD_Feas_Cut +=sum(vNewH2GenCap[h]-vRetH2GenCap[h] - value(MP[:vNewH2GenCap][h]) + value(MP[:vRetH2GenCap][h]) for h in H_ther, w in W, t in T) 
            end
            if cH2StartLimits in conflict_constraints
                BD_Feas_Cut +=sum(vNewH2GenCap[h]-vRetH2GenCap[h] - value(MP[:vNewH2GenCap][h]) + value(MP[:vRetH2GenCap][h]) for h in H_ther, w in W, t in T) 
            end
            if cH2GenRampUp in conflict_constraints
                BD_Feas_Cut +=sum(vNewH2GenCap[h]-vRetH2GenCap[h] - value(MP[:vNewH2GenCap][h]) + value(MP[:vRetH2GenCap][h]) for h in H_dis, w in W, t in 2:length(T)) 
            end
            if cH2GenRampDn in conflict_constraints
                BD_Feas_Cut +=sum(vNewH2GenCap[h]-vRetH2GenCap[h] - value(MP[:vNewH2GenCap][h]) + value(MP[:vRetH2GenCap][h]) for h in H_dis, w in W, t in 2:length(T)) 
            end
            if cMaxH2StoSOC in conflict_constraints
                BD_Feas_Cut +=sum(vNewH2StoCap[s]-vRetH2StoCap[s] - value(MP[:vNewH2GenCap][s]) + values(MP[:vRetH2StoCap][s]) for s in Q, w in W, t in 2:length(T)) 
            end
            if cMinH2StoSOC in conflict_constraints
                BD_Feas_Cut +=sum(vNewH2StoCap[s]-vRetH2StoCap[s] - value(MP[:vNewH2GenCap][s]) + values(MP[:vRetH2StoCap][s]) for s in Q, w in W, t in 2:length(T)) 
            end
            if cMaxH2StoChar in conflict_constraints
                BD_Feas_Cut +=sum(vNewH2StoCompCap[s]-vRetH2StoCompCap[s] - value(MP[:vNewH2StoCompCap][s]) + value(MP[:vRetH2StoCompCap][s]) for s in Q, w in W, t in T) 
            end
            if cMaxH2PipeFlowOut in conflict_constraints
                BD_Feas_Cut +=sum(vNewH2Pipe[i]-vRetH2Pipe[i] - value(MP[:vNewH2Pipe][i]) + value(MP[:vRetH2Pipe][i]) for i in I, w in W, t in T) 
            end
            if cMaxH2PipeFlowIn in conflict_constraints
                BD_Feas_Cut +=sum(vNewH2Pipe[i]-vRetH2Pipe[i] - value(MP[:vNewH2Pipe][i]) + value(MP[:vRetH2Pipe][i]) for i in I, w in W, t in T) 
            end
            if cMaxH2PipeFlowOutComp in conflict_constraints
                BD_Feas_Cut +=sum(vNewH2PipeCompCap[i]-vRetH2PipeCompCap[i] - value(MP[:vNewH2PipeCompCap][i]) + value(MP[:vRetH2PipeCompCap][i]) for i in I, w in W, t in T) 
            end
            if cMaxH2PipeFlowInComp in conflict_constraints
                BD_Feas_Cut +=sum(vNewH2PipeCompCap[i]-vRetH2PipeCompCap[i] - value(MP[:vNewH2PipeCompCap][i]) + value(MP[:vRetH2PipeCompCap][i]) for i in I, w in W, t in T) 
            end
            @constraint(MP, BD_Feas_Cut>= 1)
        else
        # Fallback: no conflict info or no specific constraint found
        # We can use a no-good cut to prevent this exact combination of first-stage decisions.
            @constraint(MP, sum((vNewPowGenCap[g] - vRetPowGenCap[g]) == (master_cap_new[g] - master_cap_ret[g]) for g in G) <= length(G) - 1)
            println("Added no-good cut to rule out current first-stage solution.")
        end
        
    end

    gap = UB - LB
    if gap <= tolerence * abs(UB)
        println("Converged: Gap = $(gap)")
        break
    end
#end

if iteration == max_iter && (UB - LB) > tolerence * abs(UB)
    @warn "Benders did not converge within $max_iter iterations (final gap=$(UB-LB))."
end

compute_conflict!(SP)
        conflict_constraints = []

        for (func_type, set_type) in list_of_constraint_types(SP)
            for con in all_constraints(SP, func_type, set_type)
                if MOI.get(SP, MOI.ConstraintConflictStatus(), con) == MOI.IN_CONFLICT
                    push!(conflict_constraints, con)
                    println("Conflict identified: ", con)
                end
            end
        end

        BD_Feas_Cut = 0
        if !isempty(conflict_constraints)
            for g in G_ren, w in W, t in T
                if cMaxPowGen in conflict_constraints
                    BD_Feas_Cut += sum(vNewPowGenCap[g] - vRetPowGenCap[g] - value(MP[:vNewPowGenCap][g]) + value(MP[:vRetPowGenCap][g]))
                end
            end
            if cPowOnlineUnits in conflict_constraints
                BD_Feas_Cut += sum(vNewPowGenCap[g]-vRetPowGenCap[g] - value(MP[:vNewPowGenCap][g]) + value(MP[:vRetPowGenCap][g]) for g in G_ther, w in W, t in T)
            end
            if cPowStartLimits in conflict_constraints
                BD_Feas_Cut += sum(vNewPowGenCap[g]-vRetPowGenCap[g] - value(MP[:vNewPowGenCap][g]) + value(MP[:vRetPowGenCap][g]) for g in G_ther, w in W, t in T) 
            end
            if cPowGenRampUp in conflict_constraints
                BD_Feas_Cut += sum(vNewPowGenCap[g]-vRetPowGenCap[g] - value(MP[:vNewPowGenCap][g]) + value(MP[:vRetPowGenCap][g]) for g in G_ren, w in W, t in 2:length(T)) 
            end
            if cPowGenRampDn in conflict_constraints
                BD_Feas_Cut += sum(vNewPowGenCap[g]-vRetPowGenCap[g] - value(MP[:vNewPowGenCap][g]) + value(MP[:vRetPowGenCap][g]) for g in G_ren, w in W, t in 2:length(T))
            end
            if cPowStoSOCMax in conflict_constraints
                BD_Feas_Cut += sum(vNewPowStoCap[s]-vRetPowStoCap[s] - value(MP[:vNewPowGenCap][s]) + values(MP[:vRetPowStoCap][s]) for s in S, w in W, t in 2:length(T)) 
            end
            if cPowStoSOCMin in conflict_constraints
                BD_Feas_Cut += sum(vNewPowStoCap[s]-vRetPowStoCap[s] - value(MP[:vNewPowGenCap][s]) + values(MP[:vRetPowStoCap][s]) for s in S, w in W, t in 2:length(T)) 
            end
            if cMaxPowFlowOut in conflict_constraints
                BD_Feas_Cut +=sum(vNewPowTraCap[l] - value(MP[:vNewPowTraCap][l]) for l in L, w in W, t in T)
            end
            if cMaxPowFlowIn in conflict_constraints
                BD_Feas_Cut +=sum(vNewPowTraCap[l] - value(MP[:vNewPowTraCap][l]) for l in L, w in W, t in T) 
            end
            if cMaxH2GenVar in conflict_constraints
                BD_Feas_Cut +=sum(vNewH2GenCap[h]-vRetH2GenCap[h] - value(MP[:vNewH2GenCap][h]) + value(MP[:vRetH2GenCap][h]) for h in H_dis, w in W, t in T) 
            end
            if cMinH2GenVar in conflict_constraints
                BD_Feas_Cut +=sum(vNewH2GenCap[h]-vRetH2GenCap[h] - value(MP[:vNewH2GenCap][h]) + value(MP[:vRetH2GenCap][h]) for h in H_dis, w in W, t in T)
            end
            if cH2OnlineUnits in conflict_constraints
                BD_Feas_Cut +=sum(vNewH2GenCap[h]-vRetH2GenCap[h] - value(MP[:vNewH2GenCap][h]) + value(MP[:vRetH2GenCap][h]) for h in H_ther, w in W, t in T) 
            end
            if cH2StartLimits in conflict_constraints
                BD_Feas_Cut +=sum(vNewH2GenCap[h]-vRetH2GenCap[h] - value(MP[:vNewH2GenCap][h]) + value(MP[:vRetH2GenCap][h]) for h in H_ther, w in W, t in T) 
            end
            if cH2GenRampUp in conflict_constraints
                BD_Feas_Cut +=sum(vNewH2GenCap[h]-vRetH2GenCap[h] - value(MP[:vNewH2GenCap][h]) + value(MP[:vRetH2GenCap][h]) for h in H_dis, w in W, t in 2:length(T)) 
            end
            if cH2GenRampDn in conflict_constraints
                BD_Feas_Cut +=sum(vNewH2GenCap[h]-vRetH2GenCap[h] - value(MP[:vNewH2GenCap][h]) + value(MP[:vRetH2GenCap][h]) for h in H_dis, w in W, t in 2:length(T)) 
            end
            if cMaxH2StoSOC in conflict_constraints
                BD_Feas_Cut +=sum(vNewH2StoCap[s]-vRetH2StoCap[s] - value(MP[:vNewH2GenCap][s]) + values(MP[:vRetH2StoCap][s]) for s in Q, w in W, t in 2:length(T)) 
            end
            if cMinH2StoSOC in conflict_constraints
                BD_Feas_Cut +=sum(vNewH2StoCap[s]-vRetH2StoCap[s] - value(MP[:vNewH2GenCap][s]) + values(MP[:vRetH2StoCap][s]) for s in Q, w in W, t in 2:length(T)) 
            end
            if cMaxH2StoChar in conflict_constraints
                BD_Feas_Cut +=sum(vNewH2StoCompCap[s]-vRetH2StoCompCap[s] - value(MP[:vNewH2StoCompCap][s]) + value(MP[:vRetH2StoCompCap][s]) for s in Q, w in W, t in T) 
            end
            if cMaxH2PipeFlowOut in conflict_constraints
                BD_Feas_Cut +=sum(vNewH2Pipe[i]-vRetH2Pipe[i] - value(MP[:vNewH2Pipe][i]) + value(MP[:vRetH2Pipe][i]) for i in I, w in W, t in T) 
            end
            if cMaxH2PipeFlowIn in conflict_constraints
                BD_Feas_Cut +=sum(vNewH2Pipe[i]-vRetH2Pipe[i] - value(MP[:vNewH2Pipe][i]) + value(MP[:vRetH2Pipe][i]) for i in I, w in W, t in T) 
            end
            for i in I, w in W, t in T
                if cMaxH2PipeFlowOutComp[i,w,t] in conflict_constraints
                    BD_Feas_Cut += vNewH2PipeCompCap[i]-vRetH2PipeCompCap[i] - value(MP[:vNewH2PipeCompCap][i]) + value(MP[:vRetH2PipeCompCap][i])
                end
                if cMaxH2PipeFlowInComp[i,w,t] in conflict_constraints
                    BD_Feas_Cut +=sum(vNewH2PipeCompCap[i]-vRetH2PipeCompCap[i] - value(MP[:vNewH2PipeCompCap][i]) + value(MP[:vRetH2PipeCompCap][i])) 
                end
            end
        end
=#