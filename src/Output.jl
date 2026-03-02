# ─── Post-solve: Save results ────────────────────────────────────────────────

results_dir = joinpath(datadir, "Results")
savefig(plt, joinpath(results_dir, "Convergence_LB_UB.png"))

println("Done: Iterations = $k, final gap = ", round(UB - LB, digits=4))
if UB - LB <= tolerence
    println("Optimal objective = ", round(UB, digits=2))
else
    println("Stopped without full convergence. Current best solution cost = ", round(UB, digits=2))
end

write_capacity_files(datadir)
write_line_capacity_files(datadir)
write_land_use_data(datadir)
write_nsd_by_zone(datadir, PowNSD_vals, H2NSD_vals, Z, W, T)
write_curtailment_by_zone(datadir, PowCrt_vals, H2Crt_vals, Z, W, T)
write_flows(datadir, PowFlow_vals, H2FlowPos_vals, H2FlowNeg_vals, pow_lines, hsc_pipelines, L, I, W, T)
write_emissions_detail(datadir, SP_models, Z, W)
write_storage_profiles(datadir, pow_gen, hsc_gen, S, Q,
    PowStoCha_vals, PowStoDis_vals, PowStoSOC_vals,
    H2StoCha_vals,  H2StoDis_vals,  H2StoSOC_vals, W, T)
write_generation_profiles(datadir, pow_gen, hsc_gen, G, H, PowGen_vals, H2Gen_vals, W, T)
write_costs_files(datadir, SP_models, pow_gen, hsc_gen, G, S, H, Q,
    PowGen_vals, PowStoCha_vals, H2Gen_vals, H2StoCha_vals, W, T)
write_line_costs_power(datadir, pow_lines, L)
write_h2_pipe_costs(datadir, hsc_pipelines, I, H2FlowPos_vals, H2FlowNeg_vals, W, T)
write_LCOH(datadir, SP_models, hsc_gen, H, Q, H2Gen_vals, H2StoCha_vals,
    H2FlowPos_vals, H2FlowNeg_vals, W, T, fuel_costs)
write_demand_profiles_wide(datadir, pow_demand, h2_demand, Pow_D_vals, H2_D_vals, Z, W, T)
