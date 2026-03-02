function write_storage_profiles(datadir,
    pow_gen, hsc_gen,
    S, Q,
    PowStoCha_vals, PowStoDis_vals, PowStoSOC_vals,
    H2StoCha_vals, H2StoDis_vals, H2StoSOC_vals,
    W, T
)
    results_dir = joinpath(datadir, "results")
    # Total number of hours in the year
    total_hours = length(collect(W)) * length(T)

    # Labels for Charge/Discharge files: ["Name", "Zone", "Sum", "t1", …, "t<total_hours>"]
    labels_cd = Vector{String}(undef, total_hours + 3)
    labels_cd[1] = "Name"
    labels_cd[2] = "Zone"
    labels_cd[3] = "Sum"
    for idx in 4:(total_hours + 3)
        labels_cd[idx] = "t$(idx - 3)"
    end

    # Labels for SOC files: ["Name", "Zone", "t1", …, "t<total_hours>"]
    labels_soc = Vector{String}(undef, total_hours + 2)
    labels_soc[1] = "Name"
    labels_soc[2] = "Zone"
    for idx in 3:(total_hours + 2)
        labels_soc[idx] = "t$(idx - 2)"
    end

    #
    # 1) Power Charge
    #
    df_p_ch = DataFrame(Storage = labels_cd)
    for s in sort(collect(S))
        name   = pow_gen[s, :resource]
        zone   = pow_gen[s, :zone]
        # Sum of charge over all weeks/hours
        sum_ch = 0.0
        for w in W, t in T
            sum_ch += PowStoCha_vals[(s, w, t)]
        end
        # Build column vector
        col = Vector{Any}(undef, total_hours + 3)
        col[1] = name
        col[2] = zone
        col[3] = sum_ch
        idx = 4
        for w in W, t in T
            col[idx] = PowStoCha_vals[(s, w, t)]
            idx += 1
        end
        # Use the resource name (not the ID) as column header
        df_p_ch[!, Symbol(name)] = col
    end
    CSV.write(joinpath(results_dir,"13_Power_Charge.csv"), df_p_ch)

    #
    # 2) Power Discharge
    #
    df_p_dc = DataFrame(Storage = labels_cd)
    for s in sort(collect(S))
        name   = pow_gen[s, :resource]
        zone   = pow_gen[s, :zone]
        sum_dc = 0.0
        for w in W, t in T
            sum_dc += PowStoDis_vals[(s, w, t)]
        end
        col = Vector{Any}(undef, total_hours + 3)
        col[1] = name
        col[2] = zone
        col[3] = sum_dc
        idx = 4
        for w in W, t in T
            col[idx] = PowStoDis_vals[(s, w, t)]
            idx += 1
        end
        df_p_dc[!, Symbol(name)] = col
    end
    CSV.write(joinpath(results_dir,"14_Power_Discharge.csv"), df_p_dc)

    #
    # 3) Power SOC
    #
    df_p_soc = DataFrame(Storage = labels_soc)
    for s in sort(collect(S))
        name = pow_gen[s, :resource]
        zone = pow_gen[s, :zone]
        col  = Vector{Any}(undef, total_hours + 2)
        col[1] = name
        col[2] = zone
        idx = 3
        for w in W, t in T
            col[idx] = PowStoSOC_vals[(s, w, t)]
            idx += 1
        end
        df_p_soc[!, Symbol(name)] = col
    end
    CSV.write(joinpath(results_dir,"15_Power_SOC.csv"), df_p_soc)

    #
    # 4) H₂ Charge
    #
    df_h_ch = DataFrame(Storage = labels_cd)
    for q in sort(collect(Q))
        name   = hsc_gen[q, :resource]
        zone   = hsc_gen[q, :zone]
        sum_ch = 0.0
        for w in W, t in T
            sum_ch += H2StoCha_vals[(q, w, t)]
        end
        col = Vector{Any}(undef, total_hours + 3)
        col[1] = name
        col[2] = zone
        col[3] = sum_ch
        idx = 4
        for w in W, t in T
            col[idx] = H2StoCha_vals[(q, w, t)]
            idx += 1
        end
        df_h_ch[!, Symbol(name)] = col
    end
    CSV.write(joinpath(results_dir,"16_H2_Charge.csv"), df_h_ch)

    #
    # 5) H₂ Discharge
    #
    df_h_dc = DataFrame(Storage = labels_cd)
    for q in sort(collect(Q))
        name   = hsc_gen[q, :resource]
        zone   = hsc_gen[q, :zone]
        sum_dc = 0.0
        for w in W, t in T
            sum_dc += H2StoDis_vals[(q, w, t)]
        end
        col = Vector{Any}(undef, total_hours + 3)
        col[1] = name
        col[2] = zone
        col[3] = sum_dc
        idx = 4
        for w in W, t in T
            col[idx] = H2StoDis_vals[(q, w, t)]
            idx += 1
        end
        df_h_dc[!, Symbol(name)] = col
    end
    CSV.write(joinpath(results_dir,"17_H2_Discharge.csv"), df_h_dc)

    #
    # 6) H₂ SOC
    #
    df_h_soc = DataFrame(Storage = labels_soc)
    for q in sort(collect(Q))
        name = hsc_gen[q, :resource]
        zone = hsc_gen[q, :zone]
        col  = Vector{Any}(undef, total_hours + 2)
        col[1] = name
        col[2] = zone
        idx = 3
        for w in W, t in T
            col[idx] = H2StoSOC_vals[(q, w, t)]
            idx += 1
        end
        df_h_soc[!, Symbol(name)] = col
    end
    CSV.write(joinpath(results_dir,"18_H2_SOC.csv"), df_h_soc)

end