using JuMP, Gurobi, DataFrames, CSV, Plots, Measures, Profile

const SCALE = 1e9

# ── Output writers ─────────────────────────────────────────────────────────────
include("src/Write_Capacity.jl")
include("src/Write_Land_use.jl")
include("src/Write_NSD.jl")
include("src/Write_Curtailment.jl")
include("src/Write_Flows.jl")
include("src/Write_Emissions.jl")
include("src/Write_Storage_Data.jl")
include("src/Write_Generations.jl")
include("src/Write_costs.jl")
include("src/Write_LCOH.jl")
include("src/Write_demand_profiles.jl")

# ── Model ──────────────────────────────────────────────────────────────────────
include("src/Config.jl")          # Algorithm parameters
include("src/Data_loading.jl")    # CSV loading, sets, preprocessing
include("src/Plotting.jl")        # Convergence plot initialisation
include("src/Master_problem.jl")  # Master problem (MP) definition
include("src/Sub_problem.jl")     # Sub-problems (SP) definition
include("src/Benders_loop.jl")    # Benders decomposition iteration
include("src/Output.jl")          # Save results to CSV / figures
