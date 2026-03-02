function write_demand_profiles_wide(datadir::AbstractString,
    pow_demand_static::Dict{Int,Matrix},
    h2_demand_static ::Dict{Int,Matrix},
    Pow_D_vals::Dict{Tuple{Int,Int,Int},Float64},
    H2_D_vals ::Dict{Tuple{Int,Int,Int},Float64},
    Z::AbstractVector,
    W::AbstractVector,
    T::AbstractVector)

    # ------------------------------------------------------------------
    results_dir = joinpath(datadir, "Results")
    mkpath(results_dir)

    # internal helper --------------------------------------------------
    function build_wide(static_dict, dynamic_dict)
        total_h = length(W) * length(T)
        df      = DataFrame(abs_hour = 1:total_h)

        # add two columns per zone
        for z in Z
            df[!, Symbol("Z$(z)_static")]   = zeros(Float64, total_h)
            df[!, Symbol("Z$(z)_realised")] = zeros(Float64, total_h)
        end

        row = 1
        for w in W, t in T
            for z in Z
                df[row, Symbol("Z$(z)_static")]   = static_dict[w][t, z]
                df[row, Symbol("Z$(z)_realised")] = dynamic_dict[(z, w, t)]
            end
            row += 1
        end
        return df
    end

    # ------------------------------------------------------------------
    df_pow = build_wide(pow_demand_static, Pow_D_vals)
    CSV.write(joinpath(results_dir, "26_Power_Demand.csv"), df_pow)

    df_h2  = build_wide(h2_demand_static,  H2_D_vals)
    CSV.write(joinpath(results_dir, "27_H2_Demand.csv"),   df_h2)

    println("✅  wrote 26_Power_Demand.csv and 27_H2_Demand.csv (wide format)")
end