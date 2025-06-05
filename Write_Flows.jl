function write_flows(
    PowFlow_vals::Dict{Tuple{Int,Int,Int},Float64},
    H2FlowPos_vals::Dict{Tuple{Int,Int,Int},Float64},
    H2FlowNeg_vals::Dict{Tuple{Int,Int,Int},Float64},
    pow_lines,
    hsc_pipelines,
    L::AbstractVector{<:Integer},
    I::AbstractVector{<:Integer},
    W::AbstractVector{<:Integer},
    T::AbstractVector{<:Integer}
)
    # Total number of hours in the year
    total_hours = length(collect(W)) * length(T)

    # Build the “row labels” for the first column:
    #  Row 1 = "Path Name", Row 2 = "Sum", Rows 3… = "t1", "t2", …, "t<total_hours>"
    labels = Vector{String}(undef, total_hours + 2)
    labels[1] = "Path Name"
    labels[2] = "Sum"
    for idx in 3:(total_hours + 2)
        labels[idx] = "t$(idx - 2)"
    end

    # ─────────────────────────────────────────────────────────────────────
    # 1) POWER FLOWS
    # ─────────────────────────────────────────────────────────────────────
    # Start by creating a DataFrame with a single column :Line = labels
    df_pow = DataFrame(Line = labels)

    # For each line l, construct a column of length total_hours+2:
    #   [ path_name, sum_flow, hourly flows... ]
    for l in sort(collect(L))
        # (a) Path name string
        path_name = pow_lines[l, :transmission_path_name]

        # (b) Annual sum of flows on line l
        sum_flow = 0.0
        for w in W, t in T
            sum_flow += PowFlow_vals[(l, w, t)]
        end

        # (c) Build the column vector
        col_data = Vector{Any}(undef, total_hours + 2)
        col_data[1] = path_name
        col_data[2] = sum_flow

        idx = 3
        for w in W, t in T
            col_data[idx] = PowFlow_vals[(l, w, t)]
            idx += 1
        end

        # Insert this column into df_pow, with header = string(l)
        df_pow[!, Symbol(string(l))] = col_data
    end

    CSV.write("10_Power_Flows.csv", df_pow)


    # ─────────────────────────────────────────────────────────────────────
    # 2) H₂ FLOWS
    # ─────────────────────────────────────────────────────────────────────
    df_h2 = DataFrame(Pipe = labels)

    for i in sort(collect(I))
        # (a) Path name string
        path_name = hsc_pipelines[i, :pipe_path]

        # (b) Annual net flow = sum(pos − neg)
        sum_flow = 0.0
        for w in W, t in T
            sum_flow += (H2FlowPos_vals[(i, w, t)] - H2FlowNeg_vals[(i, w, t)])
        end

        # (c) Build the column vector
        col_data = Vector{Any}(undef, total_hours + 2)
        col_data[1] = path_name
        col_data[2] = sum_flow

        idx = 3
        for w in W, t in T
            col_data[idx] = H2FlowPos_vals[(i, w, t)] - H2FlowNeg_vals[(i, w, t)]
            idx += 1
        end

        df_h2[!, Symbol(string(i))] = col_data
    end

    CSV.write("11_H2_Flows.csv", df_h2)

end