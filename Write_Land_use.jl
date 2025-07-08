function write_land_use_data(datadir)

    results_dir = joinpath(datadir, "results")
    df_sum = DataFrame(
        zone            = Int[],
        total_used      = Float64[],  # sum of new land
        total_freed     = Float64[],  # sum of retired land
        unused_land     = Float64[],  # available − (used − freed)
        pct_used        = Float64[]   # 100 * (used − freed) / available
    )

    for z in Z
        # Power generation
        new_powgen  = 0.0
        ret_powgen  = 0.0
        for g in G
            if pow_gen[g, :zone] == z
                caprep = pow_gen[g, :rep_capacity]
                lu     = pow_gen[g, :land_use_km2_p_cap]
                new_powgen += value(vNewPowGenCap[g]) * caprep * lu
                ret_powgen += value(vRetPowGenCap[g]) * caprep * lu
            end
        end

        # Power storage
        new_powsto  = 0.0
        ret_powsto  = 0.0
        for s in S
            if pow_gen[s, :zone] == z
                caprep = pow_gen[s, :rep_capacity]
                lu     = pow_gen[s, :land_use_km2_p_cap]
                new_powsto += value(vNewPowStoCap[s]) * caprep * lu
                ret_powsto += value(vRetPowStoCap[s]) * caprep * lu
            end
        end

        # H₂ generation
        new_h2gen  = 0.0
        ret_h2gen  = 0.0
        for h in H
            if hsc_gen[h, :zone] == z
                caprep = hsc_gen[h, :rep_capacity]
                lu     = hsc_gen[h, :land_use_km2_p_cap]
                new_h2gen += value(vNewH2GenCap[h]) * caprep * lu
                ret_h2gen += value(vRetH2GenCap[h]) * caprep * lu
            end
        end

        # H₂ storage
        new_h2sto  = 0.0
        ret_h2sto  = 0.0
        for q in Q
            if hsc_gen[q, :zone] == z
                caprep = hsc_gen[q, :rep_capacity]
                lu     = hsc_gen[q, :land_use_km2_p_cap]
                new_h2sto += value(vNewH2StoCap[q]) * caprep * lu
                ret_h2sto += value(vRetH2StoCap[q]) * caprep * lu
            end
        end

        # H₂ pipelines
        new_h2pipe = 0.0
        ret_h2pipe = 0.0
        for i in I
            if abs(H2_Network[i, z]) == 1
                dist  = hsc_pipelines[i, :distance]
                lup   = hsc_pipelines[i, :land_use_km2_p_length]
                new_h2pipe += 0.5 * value(vNewH2Pipe[i]) * dist * lup
                ret_h2pipe += 0.5 * value(vRetH2Pipe[i]) * dist * lup
            end
        end

        total_used   = new_powgen + new_powsto + new_h2gen + new_h2sto + new_h2pipe
        total_freed  = ret_powgen + ret_powsto + ret_h2gen + ret_h2sto + ret_h2pipe
        net_land     = total_used - total_freed

        avail        = zones[z, :available_land]
        unused       = avail - net_land
        pct          = avail > 0 ? 100 * net_land / avail : 0.0

        push!(df_sum, (z, total_used, total_freed, unused, pct))
    end

    # Write the summary to CSV
    CSV.write(joinpath(results_dir,"05_Land_use.csv"), df_sum)

    return df_sum
end