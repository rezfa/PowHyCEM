using CSV, DataFrames, CairoMakie          # plotting
using GeoMakie, GeoJSON, HTTP              # basemap + geodata
using CategoricalArrays

using CSV, DataFrames, CairoMakie          # plotting
using GeoMakie, GeoJSON, HTTP              # basemap + geodata
using CategoricalArrays

const BASE_DIR = "/Users/rez/Documents/Engineering/Coding/Julia/Revised Input/untitled folder"
const CAP_POWER = joinpath(BASE_DIR, "01_Capacity_Power.csv")
const CAP_H2    = joinpath(BASE_DIR, "02_Capacity_H2.csv")
const CENTROIDS = joinpath(BASE_DIR, "Centroids_zones.csv")
zones    = 1:14

function tidy_total(path::AbstractString)
    df = CSV.read(path, DataFrame)

    # 1) Find the row index for "Total_Cap"
    rc = findfirst(df.Attribute .== "Total_Cap")
    @assert rc !== nothing "No Total_Cap row found in $path"

    # 2) Identify the data columns whose header ends in _<digits>
    hdrs = names(df)
    # skip the first column (Attribute), look at hdrs[2:end]
    data_hdrs = hdrs[2:end]
    idxs = findall(h -> occursin(r"_\d+$", h), data_hdrs)

    # 3) Extract the values from that row, adjusting for the +1 offset
    #    since df[:, 1] is Attribute, df[:, 2:end] correspond to data_hdrs
    vals = Vector{Float64}(undef, length(idxs))
    for (i, j) in enumerate(idxs)
        vals[i] = df[rc, j+1]   # +1 because data_hdrs begin at column 2
    end

    # 4) Parse tech names and zone numbers from the headers
    parts = rsplit.(data_hdrs[idxs], "_"; limit=2)
    tech  = first.(parts)
    zone  = parse.(Int, last.(parts))

    return DataFrame(zone = zone, tech = tech, cap = vals)
end

long_pwr  = tidy_total(CAP_POWER)
long_h2   = tidy_total(CAP_H2)

wide_pwr  = unstack(long_pwr, :zone, :tech, :cap; fill = 0.0)
wide_h2   = unstack(long_h2,  :zone, :tech, :cap; fill = 0.0)

pow_cols = [:zone; Symbol.(pow_order)]
h2_cols  = [:zone; Symbol.(h2_order)]

cent = CSV.read(coords, DataFrame)
select!(cent, [:zone, :NUTS_ID])


power_pie = wide_pwr[:, pow_cols]     # 14 × 13  (zone + 12 techs)
h2_pie    = wide_h2[:,  h2_cols]

power_pie = leftjoin(power_pie, cent; on = :zone)
h2_pie = leftjoin(h2_pie, cent; on = :zone)


CSV.write("power_pie_attributes.csv", power_pie)
CSV.write("h2_pie_attributes.csv",    h2_pie)
