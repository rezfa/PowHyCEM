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
