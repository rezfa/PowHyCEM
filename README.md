<img width="1563" height="1563" alt="PowHyCEM (2)" src="https://github.com/user-attachments/assets/b8840311-f433-46c5-b963-784e626345ec" />


# PowHyCEM

PowHyCEM is a spatially aware integrated power–hydrogen capacity expansion model developed for long-term energy system planning. The model jointly optimizes investment, retirement, and hourly operation of power and hydrogen infrastructure while accounting for spatial land constraints, sector coupling, emissions, storage, and transmission.

This repository contains the model and input data used in:

**“When land becomes the bottleneck: Spatial feasibility determines decarbonization pathways.”**

## Repository structure

```text
PowHyCEM/
├── PowHyCEM_BD.jl        # Main execution file
├── Input_Data/            # Model input CSV files
└── src/                   # Model formulation, Benders algorithm, and output routines


The main components under src/ include:

Config.jl – algorithm and solver settings
Data_loading.jl – input-data loading and preprocessing
Master_problem.jl – investment/master problem
Sub_problem.jl – weekly operational subproblems
Benders_loop.jl – regularized temporal multi-cut Benders algorithm
Output.jl and Write_*.jl – result processing and export
Requirements

The model was developed and tested using:

Julia 1.9.2
Gurobi 12.0.1
JuMP
CSV
DataFrames
Plots
Measures

A valid Gurobi license is required.

Required Julia packages can be installed from the Julia package manager:
