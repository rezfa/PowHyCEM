function write_capacity_files()
    # ─────────────────────────────────────────────────────────────────────
    # 1) POWER CAPACITY (generators G and storages S)
    # ─────────────────────────────────────────────────────────────────────
    power_list = Vector{Tuple{Any,String,Float64,Float64,Float64}}()
    for g in G
        header  = pow_gen[g, :resource]
        startv  = pow_gen[g, :existing_cap]
        caprep  = pow_gen[g, :rep_capacity]
        newv    = value(vNewPowGenCap[g]) * caprep
        retv    = value(vRetPowGenCap[g]) * caprep
        totv    = value(eTotPowGenCap[g])
        push!(power_list, (g, header, startv, newv, retv))
    end
    for s in S
        header  = pow_gen[s, :resource]
        startv  = pow_gen[s, :existing_cap_mwh]
        caprep  = pow_gen[s, :rep_capacity]
        newv    = value(vNewPowStoCap[s]) * caprep
        retv    = value(vRetPowStoCap[s]) * caprep
        totv    = value(eTotPowStoCap[s])
        push!(power_list, (s, header, startv, newv, retv))
    end

    Np = length(power_list)
    attrs = ["Zone", "Starting_Cap", "New_Cap", "Retired_Cap", "Total_Cap"]
    colnames_p = ["Attribute"; [item[2] for item in power_list]...]
    names_sym_p = Symbol.(colnames_p)

    mat_p = Matrix{Any}(undef, 5, Np+1)
    mat_p .= missing
    df_p = DataFrame(mat_p, names_sym_p)

    for i in 1:5
        df_p[i, :Attribute] = attrs[i]
    end

    for (j, (r, header, startv, newv, retv)) in enumerate(power_list)
        col = Symbol(header)
        if r in G || r in S
            # Zone for both gens and storages
            z = pow_gen[r, :zone]
        else
            z = -1
        end
        totv = startv + newv - retv

        df_p[1, col] = z
        df_p[2, col] = startv
        df_p[3, col] = newv
        df_p[4, col] = retv
        df_p[5, col] = totv
    end

    CSV.write("01_Capacity_Power.csv", df_p)

    # ─────────────────────────────────────────────────────────────────────
    # 2) H₂ CAPACITY (generators H and storages Q)
    # ─────────────────────────────────────────────────────────────────────
    h2_list = Vector{Tuple{Any,String,Float64,Float64,Float64}}()
    for h in H
        header  = hsc_gen[h, :resource]
        startv  = hsc_gen[h, :existing_cap_tonne_p_hr]
        caprep  = hsc_gen[h, :rep_capacity]
        newv    = value(vNewH2GenCap[h]) * caprep
        retv    = value(vRetH2GenCap[h]) * caprep
        push!(h2_list, (h, header, startv, newv, retv))
    end
    for q in Q
        header  = hsc_gen[q, :resource]
        startv  = hsc_gen[q, :existing_cap_tonne]
        caprep  = hsc_gen[q, :rep_capacity]
        newv    = value(vNewH2StoCap[q]) * caprep
        retv    = value(vRetH2StoCap[q]) * caprep
        push!(h2_list, (q, header, startv, newv, retv))
    end
    for q in Q
        header  = string(hsc_gen[q, :resource], "_Compressor")
        startv  = hsc_gen[q, :existing_cap_comp_tonne_hr]
        newv    = value(vNewH2StoCompCap[q]) 
        retv    = value(vRetH2StoCompCap[q]) 
        push!(h2_list, (q, header, startv, newv, retv))
    end

    Nh = length(h2_list)
    colnames_h = ["Attribute"; [item[2] for item in h2_list]...]
    names_sym_h = Symbol.(colnames_h)

    mat_h = Matrix{Any}(undef, 5, Nh+1)
    mat_h .= missing
    df_h = DataFrame(mat_h, names_sym_h)

    for i in 1:5
        df_h[i, :Attribute] = attrs[i]
    end

    for (j, (r, header, startv, newv, retv)) in enumerate(h2_list)
        col = Symbol(header)
        z   = hsc_gen[r, :zone]    # same for H or Q
        totv = startv + newv - retv

        df_h[1, col] = z
        df_h[2, col] = startv
        df_h[3, col] = newv
        df_h[4, col] = retv
        df_h[5, col] = totv
    end

    CSV.write("02_Capacity_H2.csv", df_h)

    println("✅ Wrote Capacity_Power.csv and Capacity_H2.csv")
end

function write_line_capacity_files()
    # ───────────────────────────
    # Power lines
    # ───────────────────────────
    df_pl = DataFrame(
        Line    = Int[],
        Existing= Float64[],
        New     = Float64[],
        Total   = Float64[]
    )

    for l in L
        existing = pow_lines[l, :existing_transmission_cap_mw]
        newv     = value(vNewPowTraCap[l])
        # If no retirement variable exists for lines, set retired = 0.0
        total   = existing + newv 
        push!(df_pl, (l, existing, newv, total))
    end

    CSV.write("03_Line_Capacity_Power.csv", df_pl)


    # ───────────────────────────
    # H₂ pipelines
    # ───────────────────────────
    df_hp = DataFrame(
        Pipe    = Int[],
        Existing= Float64[],
        New     = Float64[],
        Retired = Float64[],
        Total   = Float64[],
        Max_Pipe_Capacity = Float64[],
        Total_Compressor_Cap = Float64[]
    )

    for i in I
        existing = hsc_pipelines[i, :existing_num_pipes]
        newv     = value(vNewH2Pipe[i])
        retired  = value(vRetH2Pipe[i])
        total    = existing + newv - retired
        total_cap = total * hsc_pipelines[i, :max_pipe_cap_tonne]
        pipe_comp = hsc_pipelines[i, :existing_comp_cap_tonne_hr] + value(vNewH2PipeCompCap[i]) - value(vRetH2PipeCompCap[i])
        push!(df_hp, (i, existing, newv, retired, total, total_cap, pipe_comp))
    end

    CSV.write("04_Pipe_Capacity_H2.csv", df_hp)

end