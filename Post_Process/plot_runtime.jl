using CSV, DataFrames, Colors, Measures
using Plots


df = CSV.read("/Users/rez/Documents/Engineering/Coding/Julia/Model_running_time.csv",
              DataFrame; missingstring = ["", "NA"])

zones = df.Zones 

log10mins(col) = map(x -> ismissing(x) ? missing : log10(x), col)

t_mono         = log10mins(df.Monolithic)          # seconds → minutes
t_tempbd       = log10mins(df.Temp_BD)
t_tempbd_lds   = log10mins(df.Temp_BD_LDS)   

iters_tempbd = df.iter_temp
iters_tempbd_lds = df.iter_LDS

avg_tempbd_sec     = df.Temp_BD          ./ iters_tempbd
avg_tempbd_lds_sec = df.Temp_BD_LDS      ./ iters_tempbd_lds

###################################################
runtime_zone = plot(size      = (800, 450),
           margins = 8mm,
           xlabel    = "Number of zones",
           ylabel    = "Log of Convergence time",
           legend    = :bottomright,
           grid      = :false,
           ylim      = (1.5, 6),
           xticks    = zones,
           guidefont  = font("Times", 11),
           tickfont   = font("Times", 9),
           legendfont = font("Times", 9))

hline!(runtime_zone, [log10(360_000)];
           linestyle = :dash,
           linewidth = 2,
           color     = :black,
           label     = "100 h limit")

plot!(runtime_zone, zones, t_mono;
      color      = :red,
      linewidth  = 2,
      marker     = :diamond,
      marker_size = 2,
      label      = "Monolithic")

plot!(runtime_zone, zones, t_tempbd;
      color      = :"blue",
      linewidth  = 2,
      marker     = :square,
      marker_size = 2,
      label      = "TempBD")

plot!(runtime_zone, zones, t_tempbd_lds;
      color      = :green,
      linewidth  = 2,
      marker     = :circle,
      marker_size = 2,
      label      = "TempBD + LDS")

######################################
iter_zone = plot(zones, iters_tempbd;
           margin = 8mm,
          size       = (800, 450),
          xlabel     = "Number of zones",
          ylabel     = "Iterations to convergence",
          legend     = :topleft,
          grid       = :false,
          xticks     = zones,
          ylim      = (0, 1200),
          guidefont  = font("Times",11),
          tickfont   = font("Times",9),
          legendfont = font("Times",9),
          color      = :blue,
          linewidth  = 2,
          marker_size = 2,
          marker     = :square,
          label      = "Temporal BD without LDS")
plot!(iter_zone, zones, iters_tempbd_lds;
          color      = :green,
          linewidth  = 2,
          marker     = :circle,
          marker_size = 2,
          label      = "Temporal BD with LDS")
#########################################
avg_iter = plot(zones, avg_tempbd_sec;
           margin = 8mm,
          size       = (800, 450),
          xlabel     = "Number of zones",
          ylabel     = "Average time per iteration (s)",
          legend     = :topleft,
          grid       = :false,
          xticks     = zones,
          ylim      = (0, 250),
          guidefont  = font("Times",11),
          tickfont   = font("Times",9),
          legendfont = font("Times",9),
          color      = :blue,
          linewidth  = 2,
          marker     = :square,
          marker_size = 2,
          label      = "TempBD without LDS")
plot!(avg_iter, zones, avg_tempbd_lds_sec;
      color     = :green,
      linewidth = 2,
      marker    = :circle,
      marker_size = 2,
      label     = "TempBD with LDS")
          
savefig(runtime_zone, "fig_runtime_vs_zones.png")
savefig(iter_zone, "fig_iterations_vs_zones.png")
savefig(avg_iter, "fig_avg_iter_time_sec_vs_zones.png")