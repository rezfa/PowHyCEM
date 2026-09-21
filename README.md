<img width="1563" height="1563" alt="PowHyCEM (2)" src="https://github.com/user-attachments/assets/b8840311-f433-46c5-b963-784e626345ec" />

# PowHyCEM

**PowHyCEM** is a spatially aware integrated power–hydrogen capacity expansion model for long-term energy system planning.

The model jointly optimizes investment, retirement, and hourly operation of power and hydrogen infrastructure while accounting for:

- power–hydrogen sector coupling;
- renewable and thermal generation;
- electricity transmission and hydrogen pipelines;
- electricity and hydrogen storage;
- emissions and carbon mitigation costs;
- hourly operational constraints; and
- spatial land availability and cumulative infrastructure land use.

This repository contains the model implementation and input data associated with the study:

> **When land becomes the bottleneck: Spatial feasibility determines decarbonization pathways**

---

## Repository Structure

```text
PowHyCEM/
│
├── PowHyCEM_BD.jl          # Main execution file
├── Input_Data/             # Input datasets for the case study
│
└── src/
    ├── Config.jl           # Solver and decomposition settings
    ├── Data_loading.jl     # Data loading and preprocessing
    ├── Master_problem.jl   # Investment/master problem
    ├── Sub_problem.jl      # Weekly operational subproblems
    ├── Benders_loop.jl     # Benders decomposition algorithm
    ├── Output.jl           # Result processing
    └── Write_*.jl          # Output-writing routines
```

---

## Requirements

The model was developed and tested using:

| Software | Version |
|---|---|
| Julia | 1.9.2 |
| Gurobi | 12.0.1 |
| JuMP | Julia optimization package |
| CSV | Julia data package |
| DataFrames | Julia data package |
| Plots | Julia plotting package |
| Measures | Julia plotting utility |

A valid **Gurobi license** is required.

The required Julia packages can be installed using:

```julia
using Pkg

Pkg.add([
    "JuMP",
    "Gurobi",
    "CSV",
    "DataFrames",
    "Plots",
    "Measures"
])
```

---

## Running PowHyCEM

### 1. Clone the repository

```bash
git clone https://github.com/rezfa/PowHyCEM.git
cd PowHyCEM
```

### 2. Check the input data

The datasets used for the case study are provided in:

```text
Input_Data/
```

For reproduction of the submitted case study, the input files can be used without modification.

### 3. Run the model

From the root directory of the repository, execute:

```bash
julia PowHyCEM_BD.jl
```

The main script performs the complete workflow:

1. loads and preprocesses the input datasets;
2. creates the investment master problem;
3. creates the weekly operational subproblems;
4. solves the model using regularized temporal multi-cut Benders decomposition;
5. checks convergence of the upper and lower bounds; and
6. processes and exports the model results.

---

## Model Configuration

Solver settings and parameters of the decomposition algorithm can be inspected or modified in:

```text
src/Config.jl
```

For reproduction of the results reported in the paper, reviewers are advised to retain the default configuration and input datasets provided in this repository.

---

## Input Data

The `Input_Data/` directory contains the main datasets required by PowHyCEM, including information on:

- power generation technologies;
- hydrogen production technologies;
- hourly electricity demand;
- hourly hydrogen demand;
- renewable generation availability;
- electricity transmission infrastructure;
- hydrogen pipeline infrastructure;
- fuel and technology costs; and
- zone-specific spatial characteristics.

The input assumptions correspond to those described in the manuscript and Supplementary Information.

---

## Computational Approach

PowHyCEM represents hourly system operation over a complete year.

To maintain tractability, the operational problem is decomposed into **52 weekly subproblems** using a regularized temporal multi-cut Benders decomposition approach.

The main problem contains long-term planning decisions, while the weekly subproblems represent detailed system operation. The weekly subproblems can be solved in parallel.

The numerical experiments reported in the paper were performed on the **Hábrók high-performance computing cluster at the University of Groningen**.

The model can also be executed on other systems with Julia and Gurobi installed; however, computational time depends strongly on:

- processor performance;
- number of available CPU cores;
- Gurobi configuration; and
- problem size.

Running the complete 14-zone, hourly case study on a standard personal computer can therefore require substantially more computation time than on an HPC system.

---

## Reproducing the Paper Results

To reproduce the submitted case study:

```text
1. Clone this repository
2. Install Julia and the required packages
3. Install and activate Gurobi
4. Keep the provided Input_Data files unchanged
5. Keep the default model configuration
6. Run PowHyCEM_BD.jl
```

The mathematical formulation of the model, spatial methodology, scenario assumptions, and decomposition algorithm are described in detail in the manuscript and its Supplementary Information.

---

## Model Outputs

The output routines contained in `src/` process the main model results, including:

- installed and retired generation capacity;
- power and hydrogen production;
- storage investment and operation;
- electricity transmission;
- hydrogen pipeline flows;
- emissions;
- system costs;
- land use; and
- other system-level and zone-level indicators.

---

## Reviewer Guidance

Reviewers interested primarily in reproducing the numerical experiments can use the following workflow:

```bash
git clone https://github.com/rezfa/PowHyCEM.git
cd PowHyCEM
julia PowHyCEM_BD.jl
```

No modification of the source code or input datasets is required to reproduce the baseline implementation.

For inspection of individual components:

```text
Input_Data/              → Input assumptions and datasets
src/Config.jl            → Solver and algorithm settings
src/Master_problem.jl    → Capacity-expansion decisions
src/Sub_problem.jl       → Hourly operational formulation
src/Benders_loop.jl      → Decomposition algorithm
PowHyCEM_BD.jl           → Complete model workflow
```

---

## Citation

If you use PowHyCEM in academic work, please cite:

> R. Fardi Asrami, A.T.D. Perera, C. Zuidema, and E. Ursavas,  
> **“When land becomes the bottleneck: Spatial feasibility determines decarbonization pathways.”**

Citation information will be updated following publication.

---

## Contact

For questions regarding the model, input data, or reproduction of the results:

**Reza Fardi Asrami**  
Department of Operations  
Faculty of Economics and Business  
University of Groningen  

📧 r.fardi.asrami@rug.nl
