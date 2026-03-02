function write_curtailment_by_zone(datadir,
    PowCrt_vals::Dict{Tuple{Int,Int,Int},Float64},
    H2Crt_vals::Dict{Tuple{Int,Int,Int},Float64},
    Z::AbstractVector{<:Integer},
    W::AbstractVector{<:Integer},
    T::AbstractVector{<:Integer}
)
    results_dir = joinpath(datadir, "Results")
    # Sort zones
    zones = sort(collect(Z))
    total_hours = length(collect(W)) * length(T)

    # Build the “Hour” column labels
    labels = Vector{String}(undef, total_hours + 1)
    labels[1] = "Annual Sum"
    for idx in 2:(total_hours + 1)
        labels[idx] = "t$(idx - 1)"
    end

    # Prepare per‐zone columns
    pow_columns = Dict{Symbol, Vector{Float64}}()
    h2_columns  = Dict{Symbol, Vector{Float64}}()

    for z in zones
        # 1) Compute annual sum
        sum_pow = 0.0
        sum_h2  = 0.0
        for w in W, t in T
            sum_pow += PowCrt_vals[(z, w, t)]
            sum_h2  += H2Crt_vals[(z, w, t)]
        end

        # 2) Build vector length = total_hours + 1
        pow_vec = Vector{Float64}(undef, total_hours + 1)
        h2_vec  = Vector{Float64}(undef, total_hours + 1)

        pow_vec[1] = sum_pow
        h2_vec[1]  = sum_h2

        # Fill hourly curtailment
        idx = 2
        for w in W, t in T
            pow_vec[idx] = PowCrt_vals[(z, w, t)]
            h2_vec[idx]  = H2Crt_vals[(z, w, t)]
            idx += 1
        end

        pow_columns[Symbol("Zone_$z")] = pow_vec
        h2_columns[Symbol("Zone_$z")]  = h2_vec
    end

    # Construct DataFrames with first column “Hour” = labels, then per‐zone
    df_pow = DataFrame(merge(Dict(:Hour => labels), pow_columns))
    df_h2  = DataFrame(merge(Dict(:Hour => labels), h2_columns))

    CSV.write(joinpath(results_dir, "08_Power_Curtailment.csv"), df_pow)
    CSV.write(joinpath(results_dir, "09_H2_Curtailment.csv"),  df_h2)
    
end