function write_land_use_data(datadir::AbstractString)
    results_dir = joinpath(datadir, "Results")
    mkpath(results_dir)

    df = DataFrame(
        zone               = Int[],
        # TOTAL land (spatial feasibility, wake/setbacks)
        total_used_km2     = Float64[],
        total_freed_km2    = Float64[],
        net_total_km2      = Float64[],
        maximum_km2        = Float64[],
        available_km2      = Float64[],
        unused_vs_max_km2  = Float64[],
        pct_total_of_max   = Float64[],
    )

    for z in Z
        # ---------- POWER GENERATION (TOTAL only) ----------
        new_powgen_total = 0.0; ret_powgen_total = 0.0
        for g in G
            if pow_gen[g,:zone] == z
                cap = pow_gen[g,:rep_capacity]
                lut = pow_gen[g,:total_land_use_km2_p_cap]
                new_powgen_total += value(vNewPowGenCap[g]) * cap * lut
                ret_powgen_total += value(vRetPowGenCap[g]) * cap * lut
            end
        end

        # ---------- POWER STORAGE (TOTAL only) ----------
        new_powsto_total = 0.0; ret_powsto_total = 0.0
        for s in S
            if pow_gen[s,:zone] == z
                cap = pow_gen[s,:rep_capacity]
                lut = pow_gen[s,:total_land_use_km2_p_cap]
                new_powsto_total += value(vNewPowStoCap[s]) * cap * lut
                ret_powsto_total += value(vRetPowStoCap[s]) * cap * lut
            end
        end

        # ---------- H2 GENERATION (TOTAL only) ----------
        new_h2gen_total = 0.0; ret_h2gen_total = 0.0
        for h in H
            if hsc_gen[h,:zone] == z
                cap = hsc_gen[h,:rep_capacity]
                lut = hsc_gen[h,:total_land_use_km2_p_cap]
                new_h2gen_total += value(vNewH2GenCap[h]) * cap * lut
                ret_h2gen_total += value(vRetH2GenCap[h]) * cap * lut
            end
        end

        # ---------- H2 STORAGE (TOTAL only) ----------
        new_h2sto_total = 0.0; ret_h2sto_total = 0.0
        for q in Q
            if hsc_gen[q,:zone] == z
                # land_use columns are per storage capacity unit
                lut = hsc_gen[q,:total_land_use_km2_p_cap]
                new_h2sto_total += value(vNewH2StoCap[q]) * lut
                ret_h2sto_total += value(vRetH2StoCap[q]) * lut
            end
        end

        # ---------- H2 PIPELINES (TOTAL only; split 50/50 to incident zones) ----------
        new_h2pipe_total = 0.0; ret_h2pipe_total = 0.0
        for i in I
            if abs(H2_Network[i,z]) == 1
                dist = hsc_pipelines[i,:distance]
                lutL = hsc_pipelines[i,:total_land_use_km2_p_length]
                # multiply land-use-per-km by distance ONCE; split 50/50
                new_h2pipe_total += 0.5 * value(vNewH2Pipe[i]) * dist * lutL
                ret_h2pipe_total += 0.5 * value(vRetH2Pipe[i]) * dist * lutL
            end
        end

        # ---------- AGGREGATION ----------
        total_used  = new_powgen_total + new_powsto_total + new_h2gen_total + new_h2sto_total + new_h2pipe_total
        total_freed = ret_powgen_total + ret_powsto_total + ret_h2gen_total + ret_h2sto_total + ret_h2pipe_total
        net_total   = total_used - total_freed

        maximum   = zones[z,:maximum_available_land]
        available = zones[z,:available_land]   # kept for info in the CSV
        unused_vs_max = max(0.0, maximum - net_total)
        pct_total_max = maximum > 0 ? 100 * net_total / maximum : 0.0

        push!(df, (z, total_used, total_freed, net_total, maximum, available, unused_vs_max, pct_total_max))
    end

    CSV.write(joinpath(results_dir, "05_Land_use.csv"), df)
    return df
end