# ─── Convergence Plot Initialisation ─────────────────────────────────────────

iters    = Int[]
LB_hist  = Float64[]
UB_hist  = Float64[]
gap_hist = Float64[]

iter_max  = isempty(iters) ? 0 : last(iters)
tick_step = 1

plt = plot(
    size       = (800, 480),
    xlabel     = "Iteration",
    ylabel     = "Cost (billion €)",
    xticks     = (1:tick_step:iter_max),
    legend     = :topright,
    ylim       = (0, 100),
    grid       = :dot,
    margin     = 8mm,
    guidefont  = font("Times", 11),
    tickfont   = font("Times", 9),
    legendfont = font("Times", 9),
)

# Initialise empty series (legend labels registered here)
plot!(plt, iters, LB_hist, label="LB", color=:blue)
plot!(plt, iters, UB_hist, label="UB", color=:red)
