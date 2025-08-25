function write_land_use_data(datadir::AbstractString;
    vLandConversion=nothing)

    results_dir = joinpath(datadir, "Results")
    mkpath(results_dir)

    df = DataFrame(
        zone                   = Int[],
        # ---- TOTAL land (spatial feasibility, wake/setbacks) ----
        total_used_km2         = Float64[],
        total_freed_km2        = Float64[],
        net_total_km2          = Float64[],
        maximum_km2            = Float64[],
        unused_vs_max_km2      = Float64[],
        pct_total_of_max       = Float64[],
        # ---- DIRECT land (costed land conversion) ----
        direct_used_km2        = Float64[],
        direct_freed_km2       = Float64[],
        net_direct_km2         = Float64[],
        available_km2          = Float64[],
        converted_km2          = Float64[],
        land_conv_cost_eur     = Float64[],
    )

    # helper to read axis safely if a JuMP container is passed
    _ax1(x) = try; axes(x, 1); catch; nothing; end
    ax_conv = vLandConversion === nothing ? nothing : _ax1(vLandConversion)

    for z in Z
        # ---------- POWER GENERATION ----------
        new_powgen_total = 0.0; ret_powgen_total = 0.0
        new_powgen_direct = 0.0; ret_powgen_direct = 0.0
        for g in G
            if pow_gen[g,:zone] == z
                cap = pow_gen[g,:rep_capacity]
                lut = pow_gen[g,:total_land_use_km2_p_cap]
                lud = pow_gen[g,:direct_land_use_km2_p_cap]
                new_powgen_total  += value(vNewPowGenCap[g]) * cap * lut
                ret_powgen_total  += value(vRetPowGenCap[g]) * cap * lut
                new_powgen_direct += value(vNewPowGenCap[g]) * cap * lud
                ret_powgen_direct += value(vRetPowGenCap[g]) * cap * lud
            end
        end

        # ---------- POWER STORAGE ----------
        new_powsto_total = 0.0; ret_powsto_total = 0.0
        new_powsto_direct = 0.0; ret_powsto_direct = 0.0
        for s in S
            if pow_gen[s,:zone] == z
                cap = pow_gen[s,:rep_capacity]
                lut = pow_gen[s,:total_land_use_km2_p_cap]
                lud = pow_gen[s,:direct_land_use_km2_p_cap]
                new_powsto_total  += value(vNewPowStoCap[s]) * cap * lut
                ret_powsto_total  += value(vRetPowStoCap[s]) * cap * lut
                new_powsto_direct += value(vNewPowStoCap[s]) * cap * lud
                ret_powsto_direct += value(vRetPowStoCap[s]) * cap * lud
            end
        end

        # ---------- H2 GENERATION ----------
        new_h2gen_total = 0.0; ret_h2gen_total = 0.0
        new_h2gen_direct = 0.0; ret_h2gen_direct = 0.0
        for h in H
            if hsc_gen[h,:zone] == z
                cap = hsc_gen[h,:rep_capacity]
                lut = hsc_gen[h,:total_land_use_km2_p_cap]
                lud = hsc_gen[h,:direct_land_use_km2_p_cap]
                new_h2gen_total  += value(vNewH2GenCap[h]) * cap * lut
                ret_h2gen_total  += value(vRetH2GenCap[h]) * cap * lut
                new_h2gen_direct += value(vNewH2GenCap[h]) * cap * lud
                ret_h2gen_direct += value(vRetH2GenCap[h]) * cap * lud
            end
        end

        # ---------- H2 STORAGE ----------
        new_h2sto_total = 0.0; ret_h2sto_total = 0.0
        new_h2sto_direct = 0.0; ret_h2sto_direct = 0.0
        for q in Q
            if hsc_gen[q,:zone] == z
                # (H2 storage 'cap' is in tonnes; land_use columns are per capacity unit)
                lut = hsc_gen[q,:total_land_use_km2_p_cap]
                lud = hsc_gen[q,:direct_land_use_km2_p_cap]
                new_h2sto_total  += value(vNewH2StoCap[q]) * lut
                ret_h2sto_total  += value(vRetH2StoCap[q]) * lut
                new_h2sto_direct += value(vNewH2StoCap[q]) * lud
                ret_h2sto_direct += value(vRetH2StoCap[q]) * lud
            end
        end

        # ---------- H2 PIPELINES (split 50/50 to incident zones) ----------
        new_h2pipe_total = 0.0; ret_h2pipe_total = 0.0
        new_h2pipe_direct = 0.0; ret_h2pipe_direct = 0.0
        for i in I
            if abs(H2_Network[i,z]) == 1
                dist = hsc_pipelines[i,:distance]
                lutL = hsc_pipelines[i,:total_land_use_km2_p_length]
                ludL = hsc_pipelines[i,:total_land_use_km2_p_length]
                # multiply land-use-per-km by distance ONCE; split 50/50 between the two zones
                new_h2pipe_total  += 0.5 * value(vNewH2Pipe[i]) * dist * lutL
                ret_h2pipe_total  += 0.5 * value(vRetH2Pipe[i]) * dist * lutL
                new_h2pipe_direct += 0.5 * value(vNewH2Pipe[i]) * dist * ludL
                ret_h2pipe_direct += 0.5 * value(vRetH2Pipe[i]) * dist * ludL
            end
        end

        # ---------- AGGREGATION ----------
        total_used  = new_powgen_total + new_powsto_total + new_h2gen_total + new_h2sto_total + new_h2pipe_total
        total_freed = ret_powgen_total + ret_powsto_total + ret_h2gen_total + ret_h2sto_total + ret_h2pipe_total
        net_total   = total_used - total_freed

        direct_used  = new_powgen_direct + new_powsto_direct + new_h2gen_direct + new_h2sto_direct + new_h2pipe_direct
        direct_freed = ret_powgen_direct + ret_powsto_direct + ret_h2gen_direct + ret_h2sto_direct + ret_h2pipe_direct
        net_direct   = direct_used - direct_freed

        maximum   = zones[z,:maximum_available_land]
        available = zones[z,:available_land]

        unused_vs_max = max(0.0, maximum - net_total)
        pct_total_max = maximum > 0 ? 100 * net_total / maximum : 0.0

        # conversion actually used (if JuMP var was passed); else compute none
        converted = (ax_conv !== nothing && z ∈ ax_conv) ? value(vLandConversion[z]) : max(0.0, net_direct - available)
        conv_cost = converted * zones[z,:land_conversion_cost_p_km2]

        push!(df, (z,
            total_used, total_freed, net_total, maximum, unused_vs_max, pct_total_max,
            direct_used, direct_freed, net_direct, available, converted, conv_cost))
    end

    CSV.write(joinpath(results_dir, "05_Land_use.csv"), df)
    return df
end