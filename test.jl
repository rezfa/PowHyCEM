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

max_iter = 1000
tolerence = 1e-3

for k in 1:max_iter

    println("Iteration:", k)
    
end
