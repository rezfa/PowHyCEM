

function load_data(datadir::AbstractString)

    # ---------- 1.  raw tables ------------------------------------------------
    pow_gen       = CSV.read(joinpath(datadir, "Powe_Gen_Data.csv"),  DataFrame)
    pow_gen_var   = CSV.read(joinpath(datadir, "Pow_Gen_Var.csv"),    DataFrame)
    pow_load      = CSV.read(joinpath(datadir, "Pow_Load.csv"),       DataFrame)
    pow_lines     = CSV.read(joinpath(datadir, "Pow_Network.csv"),    DataFrame)

    hsc_gen       = CSV.read(joinpath(datadir, "HSC_Gen_Data.csv"),   DataFrame)
    hsc_gen_var   = CSV.read(joinpath(datadir, "HSC_Gen_Var.csv"),    DataFrame)
    hsc_load      = CSV.read(joinpath(datadir, "HSC_load.csv"),       DataFrame)
    hsc_pipelines = CSV.read(joinpath(datadir, "HSC_Pipelines.csv"),  DataFrame)

    fuel          = CSV.read(joinpath(datadir, "Fuels_data.csv"),     DataFrame)
    zones         = CSV.read(joinpath(datadir, "Zone_Char.csv"),      DataFrame)

    # ---------- 2.  make headers lower-case ----------------------------------
    for f in (pow_gen, pow_load, pow_lines,
              hsc_gen, hsc_load, hsc_pipelines,
              fuel,    zones)
        rename!(f, lowercase.(names(f)))
    end

    # ---------- 3.  convenience data frames ----------------------------------
    dfGen    = pow_gen[pow_gen.gen_type .>= 1, :]
    dfStor   = pow_gen[pow_gen.stor_type .> 0, :]
    dfH2Gen  = hsc_gen[hsc_gen.h_gen_type .>= 1, :]
    dfH2Stor = hsc_gen[hsc_gen.h_stor    .>= 1, :]

    # ---------- 4.  core sets -------------------------------------------------
    G      = dfGen.r_id
    S      = dfStor.r_id
    Z      = zones.zones
    W      = 1:52
    T      = 1:168
    hours_per_week = length(T)
    H_w    = [((w-1)*hours_per_week+1):(w*hours_per_week) for w in W]

    G_ther = dfGen[dfGen.gen_type .== 1, :r_id]
    G_ren  = dfGen[dfGen.gen_type .>  1, :r_id]
    V      = dfGen[dfGen.gen_type .== 3, :r_id]
    L      = pow_lines.network_lines

    H      = dfH2Gen[dfH2Gen.h_gen_type .>= 1, :r_id]
    Q      = dfH2Stor.r_id
    I      = hsc_pipelines.hsc_pipelines
    H_ther = hsc_gen[hsc_gen.h_gen_type .== 1, :r_id]
    H_dis  = hsc_gen[hsc_gen.h_gen_type .== 2, :r_id]

    # ---------- 5.  weekly demand tensors ------------------------------------
    start_p = findfirst(==("load_mw_z1"), names(pow_load))
    pow_demand = Dict{Int,Matrix}()
    for w in W
        pow_demand[w] = Matrix(pow_load[H_w[w], start_p:start_p+length(Z)-1])
    end

    start_h  = findfirst(==("load_hsc_tonne_z1"), names(hsc_load))
    h2_demand = Dict{Int,Matrix}()
    for w in W
        h2_demand[w] = Matrix(hsc_load[H_w[w], start_h:start_h+length(Z)-1])
    end

    # ---------- 6.  fuel price & emission dictionaries -----------------------
    fuels        = names(fuel)[2:end]
    costs        = Matrix(fuel[2:end, 2:end])
    CO2_content  = fuel[1, 2:end]              # ton CO₂ / MMBtu
    fuel_costs = Dict(fuels[i] => costs[:, i]       for i in eachindex(fuels))
    fuel_CO2   = Dict(fuels[i] => CO2_content[i]    for i in eachindex(fuels))

    # ---------- 7.  network incidence matrices -------------------------------
    col_p       = findfirst(==("z1"), names(pow_lines))
    Pow_Network = pow_lines[1:length(L), col_p:col_p+length(Z)-1]

    col_h       = findfirst(==("z1"), names(hsc_pipelines))
    H2_Network  = hsc_pipelines[1:length(I), col_h:col_h+length(Z)-1]

    # ---------- 8.  bundle & return ------------------------------------------
    return (; pow_gen, pow_gen_var, pow_load, pow_lines,
            hsc_gen, hsc_gen_var, hsc_load, hsc_pipelines,
            fuel, zones,
            dfGen, dfStor, dfH2Gen, dfH2Stor,
            G, S, Z, W, T, H_w,
            G_ther, G_ren, V, L,
            H, Q, I, H_ther, H_dis,
            pow_demand, h2_demand,
            fuels, fuel_costs, fuel_CO2, CO2_content,
            Pow_Network, H2_Network)
end
