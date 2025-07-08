function write_generation_profiles(datadir,
    pow_gen, hsc_gen,
    G, H,
    PowGen_vals, H2Gen_vals,
    W, T
)
    results_dir = joinpath(datadir, "Results")
    # Total number of hours in the year
    total_hours = length(collect(W)) * length(T)

    # Build the first‐column labels: ["Zone", "Sum", "t1", …, "t<total_hours>"]
    labels = Vector{String}(undef, total_hours + 2)
    labels[1] = "Zone"
    labels[2] = "Sum"
    for idx in 3:(total_hours + 2)
        labels[idx] = "t$(idx - 2)"
    end

    #
    # 1) Power generation by source
    #
    df_p_gen = DataFrame(Resource = labels)
    for g in sort(collect(G))
        name   = pow_gen[g, :resource]
        zone   = pow_gen[g, :zone]

        # Sum of generation over all weeks/hours
        sum_gen = 0.0
        for w in W, t in T
            sum_gen += PowGen_vals[(g, w, t)]
        end

        # Build column vector of length total_hours + 2
        col = Vector{Any}(undef, total_hours + 2)
        col[1] = zone
        col[2] = sum_gen

        idx = 3
        for w in W, t in T
            col[idx] = PowGen_vals[(g, w, t)]
            idx += 1
        end

        # Use the resource name as column header
        df_p_gen[!, Symbol(name)] = col
    end
    CSV.write(joinpath(results_dir,"19_Power_Gen.csv"), df_p_gen)

    #
    # 2) H₂ generation by source
    #
    df_h2_gen = DataFrame(Resource = labels)
    for h in sort(collect(H))
        name   = hsc_gen[h, :resource]
        zone   = hsc_gen[h, :zone]

        # Sum of generation over all weeks/hours
        sum_gen = 0.0
        for w in W, t in T
            sum_gen += H2Gen_vals[(h, w, t)]
        end

        col = Vector{Any}(undef, total_hours + 2)
        col[1] = zone
        col[2] = sum_gen

        idx = 3
        for w in W, t in T
            col[idx] = H2Gen_vals[(h, w, t)]
            idx += 1
        end

        df_h2_gen[!, Symbol(name)] = col
    end
    CSV.write(joinpath(results_dir,"20_H2_Gen.csv"), df_h2_gen)


end