function write_nsd_by_zone(
    PowNSD_vals::Dict{Tuple{Int,Int,Int},Float64},
    H2NSD_vals::Dict{Tuple{Int,Int,Int},Float64},
    Z::AbstractVector{<:Integer},
    W::AbstractVector{<:Integer},
    T::AbstractVector{<:Integer}
)
    # Ensure zones sorted
    zones = sort(collect(Z))
    hours_per_week = length(T)
    num_hours = length(collect(W)) * hours_per_week

    # Build the “Zone” label column: length = num_hours + 1
    labels = Vector{String}(undef, num_hours + 1)
    labels[1] = "Annual Sum"
    for idx in 2:(num_hours+1)
        labels[idx] = "t$(idx-1)"
    end

    # Prepare per‐zone NSD columns
    pow_columns = Dict{Symbol, Vector{Float64}}()
    h2_columns  = Dict{Symbol, Vector{Float64}}()

    for z in zones
        # (1) Compute annual sum
        sum_pow = 0.0
        sum_h2  = 0.0
        for w in W, t in T
            sum_pow += PowNSD_vals[(z, w, t)]
            sum_h2  += H2NSD_vals[(z, w, t)]
        end

        # (2) Build a vector of length num_hours+1
        pow_vec = Vector{Float64}(undef, num_hours + 1)
        h2_vec  = Vector{Float64}(undef, num_hours + 1)

        pow_vec[1] = sum_pow
        h2_vec[1]  = sum_h2

        # Fill rows 2…(num_hours+1) with hourly NSD
        idx = 2
        for w in W, t in T
            pow_vec[idx] = PowNSD_vals[(z, w, t)]
            h2_vec[idx]  = H2NSD_vals[(z, w, t)]
            idx += 1
        end

        pow_columns[Symbol("Zone_$z")] = pow_vec
        h2_columns[Symbol("Zone_$z")]  = h2_vec
    end

    # Construct DataFrames with first column “Zone” = labels, then per‐zone columns
    df_pow = DataFrame(merge(Dict(:Zone => labels), pow_columns))
    df_h2  = DataFrame(merge(Dict(:Zone => labels), h2_columns))

    CSV.write("06_Power_NSD.csv", df_pow)
    CSV.write("07_H2_NSD.csv",  df_h2)

end