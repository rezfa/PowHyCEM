# ─── Benders Decomposition Loop ──────────────────────────────────────────────

# ---- Helper: toggle investment variable integrality ----------------------
function toggle_integrality!(make_integer::Bool)
    containers = (
        vNewPowGenCap,    vRetPowGenCap,
        vNewPowStoCap,    vRetPowStoCap,
        vNewPowTraCap,
        vNewH2GenCap,     vRetH2GenCap,
        vNewH2StoCap,     vRetH2StoCap,
        vNewH2Pipe,       vRetH2Pipe,
        vNewH2PipeCompCap, vRetH2PipeCompCap,
        vNewH2StoCompCap,  vRetH2StoCompCap,
    )
    for c in containers, v in c
        if make_integer
            !JuMP.is_integer(v) && JuMP.set_integer(v)
        else
            JuMP.is_integer(v)  && JuMP.unset_integer(v)
        end
    end
    INTEGER_MODE[] = make_integer
    return
end

toggle_integrality!(false)

# ---- Solve initial Master Problem to get a starting investment plan ------
optimize!(MP)
@assert termination_status(MP) == MOI.OPTIMAL
LB = objective_value(MP)
println("Initial MP objective (LB) = ", round(LB, digits=2))

# ---- Result storage dictionaries -----------------------------------------
coupling = Dict{Int, Dict{Symbol, Any}}()
for w in W
    coupling[w] = Dict(
        :gencap     => ConstraintRef[],
        :genunit    => ConstraintRef[],
        :stocap     => ConstraintRef[],
        :tracap     => ConstraintRef[],
        :h2gen      => ConstraintRef[],
        :h2genunit  => ConstraintRef[],
        :h2stocap   => ConstraintRef[],
        :h2stocomp  => ConstraintRef[],
        :h2pipe     => ConstraintRef[],
        :h2pipecomp => ConstraintRef[],
        :emission   => ConstraintRef[],
        :h2SOCfirst => ConstraintRef[],
        :h2SOClast  => ConstraintRef[],
    )
end

PowNSD_vals    = Dict{Tuple{Int,Int,Int}, Float64}()
H2NSD_vals     = Dict{Tuple{Int,Int,Int}, Float64}()
PowCrt_vals    = Dict{Tuple{Int,Int,Int}, Float64}()
H2Crt_vals     = Dict{Tuple{Int,Int,Int}, Float64}()
PowGen_vals    = Dict{Tuple{Int,Int,Int}, Float64}()
H2Gen_vals     = Dict{Tuple{Int,Int,Int}, Float64}()
PowStoCha_vals = Dict{Tuple{Int,Int,Int}, Float64}()
H2StoCha_vals  = Dict{Tuple{Int,Int,Int}, Float64}()
PowStoDis_vals = Dict{Tuple{Int,Int,Int}, Float64}()
H2StoDis_vals  = Dict{Tuple{Int,Int,Int}, Float64}()
PowStoSOC_vals = Dict{Tuple{Int,Int,Int}, Float64}()
H2StoSOC_vals  = Dict{Tuple{Int,Int,Int}, Float64}()
PowFlow_vals   = Dict{Tuple{Int,Int,Int}, Float64}()
H2FlowPos_vals = Dict{Tuple{Int,Int,Int}, Float64}()
H2FlowNeg_vals = Dict{Tuple{Int,Int,Int}, Float64}()
Pow_D_vals     = Dict{Tuple{Int,Int,Int}, Float64}()
H2_D_vals      = Dict{Tuple{Int,Int,Int}, Float64}()

# ---- Initialise coupling constraints from first MP solution --------------
for w in W
    cc = coupling[w]

    cAvailPowGenCap  = @constraint(SP_models[w], [g in G],     SP_models[w][:eAvailPowGenCap][g]  == value(eTotPowGenCap[g]))
    cAvailPowGenUnit = @constraint(SP_models[w], [g in G_ther], SP_models[w][:eAvailPowGenUnit][g] == pow_gen[g, :num_units] + (value(vNewPowGenCap[g]) - value(vRetPowGenCap[g])))
    cAvailPowStoCap  = @constraint(SP_models[w], [s in S],     SP_models[w][:eAvailPowStoCap][s]  == value(eTotPowStoCap[s]))
    cAvailPowTraCap  = @constraint(SP_models[w], [l in L],     SP_models[w][:eAvailPowTraCap][l]  == value(vNewPowTraCap[l]) + pow_lines[l, :existing_transmission_cap_mw])
    cAvailH2GenCap   = @constraint(SP_models[w], [h in H_dis], SP_models[w][:eAvailH2GenCap][h]   == value(eTotH2GenCap[h]))
    cAvailH2GenUnit  = @constraint(SP_models[w], [h in H_ther], SP_models[w][:eAvailH2GenUnit][h] == hsc_gen[h, :num_units] + (value(vNewH2GenCap[h]) - value(vRetH2GenCap[h])))
    cAvailH2StoCap   = @constraint(SP_models[w], [s in Q],     SP_models[w][:eAvailH2StoCap][s]   == value(eTotH2StoCap[s]))
    cAvailH2Pipe     = @constraint(SP_models[w], [i in I],     SP_models[w][:eAvailH2Pipe][i]     == value(eTotH2Pipe[i]))
    cAvailH2PipeComp = @constraint(SP_models[w], [i in I],     SP_models[w][:eAvailH2PipeCompCap][i] == value(eTotH2PipeCompCap[i]))
    cAvailH2StoComp  = @constraint(SP_models[w], [s in Q],     SP_models[w][:eAvailH2StoCompCap][s]  == value(eTotH2StoCompCap[s]))
    cEmission        = @constraint(SP_models[w],               SP_models[w][:eMaxEmissionByWeek]   == value(vMaxEmissionByWeek[w]))
    cH2SOCFirst      = @constraint(SP_models[w], [s in Q],     SP_models[w][:eH2SOCFirst][s]       == value(vH2SOCFirst[s, w]))
    cH2SOCLast       = @constraint(SP_models[w], [s in Q],     SP_models[w][:eH2SOCLast][s]        == value(vH2SOCLast[s, w]))

    cc[:gencap]     = Dict(g => cAvailPowGenCap[g]   for g in G)
    cc[:genunit]    = Dict(g => cAvailPowGenUnit[g]  for g in G_ther)
    cc[:stocap]     = Dict(s => cAvailPowStoCap[s]   for s in S)
    cc[:tracap]     = Dict(l => cAvailPowTraCap[l]   for l in L)
    cc[:h2gen]      = Dict(h => cAvailH2GenCap[h]    for h in H_dis)
    cc[:h2genunit]  = Dict(h => cAvailH2GenUnit[h]   for h in H_ther)
    cc[:h2stocap]   = Dict(s => cAvailH2StoCap[s]    for s in Q)
    cc[:h2pipe]     = Dict(i => cAvailH2Pipe[i]      for i in I)
    cc[:h2pipecomp] = Dict(i => cAvailH2PipeComp[i]  for i in I)
    cc[:h2stocomp]  = Dict(s => cAvailH2StoComp[s]   for s in Q)
    cc[:h2SOCfirst] = Dict(s => cH2SOCFirst[s]       for s in Q)
    cc[:h2SOClast]  = Dict(s => cH2SOCLast[s]        for s in Q)
    cc[:emission]   = cEmission
end

# ---- Main Benders iteration loop -----------------------------------------
for k in 1:max_iter

    println("────────────────────────────────────────")
    println(" BENDERS ITERATION $k")
    println("────────────────────────────────────────")

    # Snapshot current MP solution
    vNewPowGenCap_val     = value.(vNewPowGenCap)
    vRetPowGenCap_val     = value.(vRetPowGenCap)
    vNewPowStoCap_val     = value.(vNewPowStoCap)
    vRetPowStoCap_val     = value.(vRetPowStoCap)
    vNewPowTraCap_val     = value.(vNewPowTraCap)
    vNewH2GenCap_val      = value.(vNewH2GenCap)
    vRetH2GenCap_val      = value.(vRetH2GenCap)
    vNewH2StoCap_val      = value.(vNewH2StoCap)
    vRetH2StoCap_val      = value.(vRetH2StoCap)
    vNewH2Pipe_val        = value.(vNewH2Pipe)
    vRetH2Pipe_val        = value.(vRetH2Pipe)
    vNewH2PipeCompCap_val = value.(vNewH2PipeCompCap)
    vRetH2PipeCompCap_val = value.(vRetH2PipeCompCap)
    vNewH2StoCompCap_val  = value.(vNewH2StoCompCap)
    vRetH2StoCompCap_val  = value.(vRetH2StoCompCap)
    vMaxEmissionByWeek_val = value.(vMaxEmissionByWeek)
    vH2SOCFirst_val       = value.(vH2SOCFirst)
    vH2SOCLast_val        = value.(vH2SOCLast)

    # Solve all sub-problems in parallel
    Threads.@threads for w in W
        optimize!(SP_models[w])
    end

    for w in W
        @assert termination_status(SP_models[w]) == MOI.OPTIMAL
    end

    total_sp_cost  = sum(objective_value(SP_models[w]) for w in W)
    invest_cost    = value(MP_obj)
    UB_candidate   = invest_cost + total_sp_cost
    global UB      = min(UB, UB_candidate)

    println(" → Total SP cost  = ", round(total_sp_cost, digits=2))
    println(" → Candidate UB   = ", round(UB_candidate,  digits=2),
            " (Best UB = ",         round(UB,            digits=2), ")")

    # Add Benders optimality cuts to MP
    for w in W
        @constraint(MP,
            theta[w] >= objective_value(SP_models[w]) +
            sum(dual(coupling[w][:gencap][g])    * pow_gen[g, :rep_capacity] *
                (vNewPowGenCap[g] - vRetPowGenCap[g] - vNewPowGenCap_val[g] + vRetPowGenCap_val[g])
                for g in G_ren) +
            sum(dual(coupling[w][:genunit][g])   *
                (vNewPowGenCap[g] - vRetPowGenCap[g] - vNewPowGenCap_val[g] + vRetPowGenCap_val[g])
                for g in G_ther) +
            sum(dual(coupling[w][:stocap][s])    * pow_gen[s, :rep_capacity] *
                (vNewPowStoCap[s] - vRetPowStoCap[s] - vNewPowStoCap_val[s] + vRetPowStoCap_val[s])
                for s in S) +
            sum(dual(coupling[w][:tracap][l])    *
                (vNewPowTraCap[l] - vNewPowTraCap_val[l])
                for l in L) +
            sum(dual(coupling[w][:h2gen][h])     * hsc_gen[h, :rep_capacity] *
                (vNewH2GenCap[h] - vRetH2GenCap[h] - vNewH2GenCap_val[h] + vRetH2GenCap_val[h])
                for h in H_dis) +
            sum(dual(coupling[w][:h2genunit][h]) *
                (vNewH2GenCap[h] - vRetH2GenCap[h] - vNewH2GenCap_val[h] + vRetH2GenCap_val[h])
                for h in H_ther) +
            sum(dual(coupling[w][:h2stocap][s])  *
                (vNewH2StoCap[s] - vRetH2StoCap[s] - vNewH2StoCap_val[s] + vRetH2StoCap_val[s])
                for s in Q) +
            sum(dual(coupling[w][:h2pipe][i])    *
                (vNewH2Pipe[i] - vRetH2Pipe[i] - vNewH2Pipe_val[i] + vRetH2Pipe_val[i])
                for i in I) +
            sum(dual(coupling[w][:h2pipecomp][i]) *
                (vNewH2PipeCompCap[i] - vRetH2PipeCompCap[i] - vNewH2PipeCompCap_val[i] + vRetH2PipeCompCap_val[i])
                for i in I) +
            sum(dual(coupling[w][:h2stocomp][s]) *
                (vNewH2StoCompCap[s] - vRetH2StoCompCap[s] - vNewH2StoCompCap_val[s] + vRetH2StoCompCap_val[s])
                for s in Q) +
            sum(dual(coupling[w][:h2SOCfirst][s]) *
                (vH2SOCFirst[s, w] - vH2SOCFirst_val[s, w])
                for s in Q) +
            sum(dual(coupling[w][:h2SOClast][s])  *
                (vH2SOCLast[s, w] - vH2SOCLast_val[s, w])
                for s in Q) +
            dual(coupling[w][:emission]) * (vMaxEmissionByWeek[w] - vMaxEmissionByWeek_val[w])
        )
    end

    set_objective_function(MP, MP_obj + sum(theta[w] for w in W))
    optimize!(MP)
    println("MP status after adding cuts: ", termination_status(MP))
    @assert termination_status(MP) == MOI.OPTIMAL
    global LB = objective_value(MP)
    println(" → Updated MP objective (LB) = ", round(LB, digits=2))
    println(" → Gap = ", round((UB - LB) * 100 / abs(LB), digits=2), "%")

    # === Continuous (LP relaxation) stage ====================================
    if !INTEGER_MODE[]
        rhs = LB + ALPHA * (UB - LB)
        if LEVELSET_CON[] === nothing
            LEVELSET_CON[] = @constraint(MP, MP_obj + sum(theta[w] for w in W) <= rhs)
        else
            JuMP.set_normalized_rhs(LEVELSET_CON[], rhs)
        end

        JuMP.set_objective_function(MP, φ)
        optimize!(MP)

        # Switch to MILP stage when gap < 1 %
        if (UB - LB) / abs(LB) < 0.01
            println("→  switching to integer mode (gap < 1 %)")
            delete(MP, LEVELSET_CON[])
            LEVELSET_CON[] = nothing
            JuMP.set_objective_function(MP, MP_obj + sum(theta[w] for w in W))
            toggle_integrality!(true)
            INTEGER_MODE[] = true
            set_optimizer_attribute(MP, "Method",   2)
            set_optimizer_attribute(MP, "Crossover", 0)
            set_optimizer_attribute(MP, "MIPGap",   1e-2)
            set_optimizer_attribute(MP, "BarConvTol",1e-2)
            optimize!(MP)
        end

    # === Integer (MILP) stage ================================================
    else
        optimize!(MP)

        if (UB - LB) / max(1e-6, abs(LB)) <= tolerence
            println("Converged (gap = ", UB - LB, "). Optimal investment plan found.")

            # Fix sub-problems to optimal investment plan and re-solve
            for w in W
                for g in G
                    set_normalized_rhs(coupling[w][:gencap][g], value(eTotPowGenCap[g]))
                end
                for g in G_ther
                    set_normalized_rhs(coupling[w][:genunit][g],
                        pow_gen[g, :num_units] + (value(vNewPowGenCap[g]) - value(vRetPowGenCap[g])))
                end
                for s in S
                    set_normalized_rhs(coupling[w][:stocap][s], value(eTotPowStoCap[s]))
                end
                for l in L
                    set_normalized_rhs(coupling[w][:tracap][l],
                        value(vNewPowTraCap[l]) + pow_lines[l, :existing_transmission_cap_mw])
                end
                for h in H_dis
                    set_normalized_rhs(coupling[w][:h2gen][h], value(eTotH2GenCap[h]))
                end
                for h in H_ther
                    set_normalized_rhs(coupling[w][:h2genunit][h],
                        hsc_gen[h, :num_units] + (value(vNewH2GenCap[h]) - value(vRetH2GenCap[h])))
                end
                for s in Q
                    set_normalized_rhs(coupling[w][:h2stocap][s], value(eTotH2StoCap[s]))
                end
                for i in I
                    set_normalized_rhs(coupling[w][:h2pipe][i], value(eTotH2Pipe[i]))
                end
                for i in I
                    set_normalized_rhs(coupling[w][:h2pipecomp][i], value(eTotH2PipeCompCap[i]))
                end
                for s in Q
                    set_normalized_rhs(coupling[w][:h2stocomp][s], value(eTotH2StoCompCap[s]))
                end
                for s in Q
                    set_normalized_rhs(coupling[w][:h2SOCfirst][s], value(vH2SOCFirst[s, w]))
                end
                for s in Q
                    set_normalized_rhs(coupling[w][:h2SOClast][s], value(vH2SOCLast[s, w]))
                end
                set_normalized_rhs(coupling[w][:emission], value(vMaxEmissionByWeek[w]))
            end

            Threads.@threads for w in W
                optimize!(SP_models[w])
            end

            # Collect final solution values
            for w in W
                sp = SP_models[w]
                for z in Z, t in T
                    PowNSD_vals[(z, w, t)] = value(sp[:vPowNSD][z, t])
                    H2NSD_vals[(z, w, t)]  = value(sp[:vH2NSD][z, t])
                    PowCrt_vals[(z, w, t)] = value(sp[:vPowCrt][z, t])
                    H2Crt_vals[(z, w, t)]  = value(sp[:vH2Crt][z, t])
                    Pow_D_vals[(z, w, t)]  = value(sp[:pow_D][t, z])
                    H2_D_vals[(z, w, t)]   = value(sp[:H2_D][t, z])
                end
                for g in G, t in T
                    PowGen_vals[(g, w, t)] = value(sp[:vPowGen][g, t])
                end
                for h in H, t in T
                    H2Gen_vals[(h, w, t)] = value(sp[:vH2Gen][h, t])
                end
                for s in S, t in T
                    PowStoCha_vals[(s, w, t)] = value(sp[:vPowStoCha][s, t])
                    PowStoDis_vals[(s, w, t)] = value(sp[:vPowStoDis][s, t])
                    PowStoSOC_vals[(s, w, t)] = value(sp[:vPowSOC][s, t])
                end
                for s in Q, t in T
                    H2StoCha_vals[(s, w, t)] = value(sp[:vH2StoCha][s, t])
                    H2StoDis_vals[(s, w, t)] = value(sp[:vH2StoDis][s, t])
                    H2StoSOC_vals[(s, w, t)] = value(sp[:vH2StoSOC][s, t])
                end
                for l in L, t in T
                    PowFlow_vals[(l, w, t)] = value(sp[:vPowFlow][l, t])
                end
                for i in I, t in T
                    H2FlowPos_vals[(i, w, t)] = value(sp[:vH2FlowPos][i, t])
                    H2FlowNeg_vals[(i, w, t)] = value(sp[:vH2FlowNeg][i, t])
                end
            end
            break
        end
    end

    push!(iters, k)
    push!(LB_hist, LB / SCALE)
    push!(UB_hist, UB / SCALE)

    if k % 10 == 0 || (UB - LB) / abs(LB + eps()) <= tolerence
        plot!(plt, iters, LB_hist, color=:blue, label=false)
        plot!(plt, iters, UB_hist, color=:red,  label=false)
        display(plt)
        println("Iter $k ▶  LB = $(round(LB/1e6, digits=3))e6, UB = $(round(UB/1e6, digits=3))e6, ",
                "rel gap = $(round((UB - LB) / abs(LB + eps()) * 100, digits=3))%")
    end
    k += 1

    # Update coupling constraint RHS for next iteration
    for w in W
        for g in G
            set_normalized_rhs(coupling[w][:gencap][g], value(eTotPowGenCap[g]))
        end
        for g in G_ther
            set_normalized_rhs(coupling[w][:genunit][g],
                pow_gen[g, :num_units] + (value(vNewPowGenCap[g]) - value(vRetPowGenCap[g])))
        end
        for s in S
            set_normalized_rhs(coupling[w][:stocap][s], value(eTotPowStoCap[s]))
        end
        for l in L
            set_normalized_rhs(coupling[w][:tracap][l],
                value(vNewPowTraCap[l]) + pow_lines[l, :existing_transmission_cap_mw])
        end
        for h in H_dis
            set_normalized_rhs(coupling[w][:h2gen][h], value(eTotH2GenCap[h]))
        end
        for h in H_ther
            set_normalized_rhs(coupling[w][:h2genunit][h],
                hsc_gen[h, :num_units] + (value(vNewH2GenCap[h]) - value(vRetH2GenCap[h])))
        end
        for s in Q
            set_normalized_rhs(coupling[w][:h2stocap][s], value(eTotH2StoCap[s]))
        end
        for i in I
            set_normalized_rhs(coupling[w][:h2pipe][i], value(eTotH2Pipe[i]))
        end
        for i in I
            set_normalized_rhs(coupling[w][:h2pipecomp][i], value(eTotH2PipeCompCap[i]))
        end
        for s in Q
            set_normalized_rhs(coupling[w][:h2stocomp][s], value(eTotH2StoCompCap[s]))
        end
        for s in Q
            set_normalized_rhs(coupling[w][:h2SOCfirst][s], value(vH2SOCFirst[s, w]))
        end
        for s in Q
            set_normalized_rhs(coupling[w][:h2SOClast][s], value(vH2SOCLast[s, w]))
        end
        set_normalized_rhs(coupling[w][:emission], value(vMaxEmissionByWeek[w]))
    end

end
