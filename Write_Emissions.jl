function write_emissions_detail(
    SP_models::AbstractVector{<:Model},
    Z::AbstractVector{<:Integer},
    W::AbstractVector{<:Integer}
    )

    nW = length(W)
    week_labels = Vector{String}(undef, nW + 1)
    week_labels[1] = "Sum"
    for (i, w) in enumerate(W)
        week_labels[i + 1] = "W$(w)"
    end
    zone_cols = Dict{Symbol, Vector{Float64}}()
    for z in Z
        col = Vector{Float64}(undef, nW + 1)
        # Sum over weeks for zone z
        sum_z = 0.0
        for w in W
            val = value(SP_models[w][:eEmissionByWeekZone][z])
            sum_z += val
        end
        col[1] = sum_z
        # Fill per-week entries
        for (i, w) in enumerate(W)
            col[i + 1] = value(SP_models[w][:eEmissionByWeekZone][z])
        end
        zone_cols[Symbol("Zone_$(z)")] = col
    end

    # 3) Build the "Weekly_Emission", "Extra_Emission", "Emission_Cost" columns
    weekly_em_col  = Vector{Float64}(undef, nW + 1)
    extra_em_col   = Vector{Float64}(undef, nW + 1)
    cost_em_col    = Vector{Float64}(undef, nW + 1)

    # Compute sums first
    sum_weekly = 0.0
    sum_extra  = 0.0
    sum_cost   = 0.0
    for w in W
        sum_weekly += value(SP_models[w][:eEmissionByWeek])
        sum_extra  += value(SP_models[w][:vExtraEmission])
        sum_cost   += value(SP_models[w][:eEmissionCost])
    end
    weekly_em_col[1] = sum_weekly
    extra_em_col[1]  = sum_extra
    cost_em_col[1]   = sum_cost

    # Fill per-week data
    for (i, w) in enumerate(W)
        weekly_em_col[i + 1] = value(SP_models[w][:eEmissionByWeek])
        extra_em_col[i + 1]  = value(SP_models[w][:vExtraEmission])
        cost_em_col[i + 1]   = value(SP_models[w][:eEmissionCost])
    end

    # 4) Assemble into a DataFrame
    # Columns in order: "Week", then each "Zone_<z>" for z∈Z, then 
    # "Weekly_Emission", "Extra_Emission", "Emission_Cost"
    df = DataFrame(Week = week_labels)
    for z in Z
        df[!, Symbol("Zone_$(z)")] = zone_cols[Symbol("Zone_$(z)")]
    end
    df[!, :Weekly_Emission] = weekly_em_col
    df[!, :Extra_Emission]  = extra_em_col
    df[!, :Emission_Cost]   = cost_em_col

    # 5) Write to CSV
    CSV.write("12_Emissions.csv", df)

end