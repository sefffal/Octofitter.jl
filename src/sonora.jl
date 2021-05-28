
# function sonora()

# end
using DelimitedFiles, NamedTupleTools, ScatteredInterpolation
function load_table(fname=joinpath(@__DIR__, "sonora_flux_table.txt"))

    headers = open(fname, lock=false, read=true) do f
        readline(f)
        readline(f)
        readline(f)
        readline(f)
        readline(f)
        h1 = readline(f)
        h2 = readline(f)
        headers_1 = strip.(split(h1, '|'))
        headers_2 = strip.(split(h2, r"  +"))

        # Starting after 6
        headers_2_1_indices = [
            2, 2, 2, 2, 2, 2, 2,
            3, 3, 3, 
            4, 4, 4,
            5, 5, 5, 5,
            6,6,6,6,
            7, 7, 7, 7,
        ]

        combined = vcat(headers_2[1:6], headers_1[headers_2_1_indices] .* '_' .*  headers_2[7:end])
        # Normalize headers
        combined = replace.(combined, '\''=>'′')
        combined = replace.(combined, ' '=>"")
        combined = replace.(combined, "/"=>'_')
        combined = replace.(combined, "2MASS"=>"TwoMASS")
        combined = replace.(combined, r"[^\w′_]"=>"")
        return combined
    end

    data = readdlm(fname, String, skipstart=7, header=false)
    data = [
        try
            typeof(d) <:AbstractString ? parse(Float64, replace(d, '*'=>"")) : d
        catch err
            NaN
        end
        for d in data
    ]
    # Return simple table
    return namedtuple(headers, eachcol(data));
end

# Prepare a set of 2D interpolations from Teff, mass -> flux in all the different bands
# in the sonora models table.
# The results are returned in mJy at 10px (not logged)


local sonora_table
function __init__()
    @info "Loading Sonora model table"
    global sonora_table = load_table()
end

function make_itp(itp)
    function (Teff,mass)
        if 0.53 ≤ mass ≤ 98 && 200 ≤ Teff ≤ 2400
            in = @SArray[log10(Teff),mass]
            return 10^only(evaluate(itp, in))
        else
            return NaN
        end
    end
end

export sonora_interpolator
function sonora_interpolator(key)

    points = vcat(
        log10.((sonora_table.Teff))',
        sonora_table.mass',
    )

    # filtered_keys = filter(keys(sonora)) do key
    #     key ∉ (:Teff, :logg, :mass, :R_Rsun, :Y, :logKzz)
    # end

    # interpolators = map(filtered_keys) do key
    #     itp = ScatteredInterpolation.interpolate(ScatteredInterpolation.ThinPlate(), points, getproperty(sonora, key))
    #     return make_itp(itp)
    # end

    itp = ScatteredInterpolation.interpolate(ScatteredInterpolation.ThinPlate(), points, getproperty(sonora_table, key))
    return make_itp(itp)

    # global sonora_flux_interp = namedtuple(filtered_keys, interpolators)
end



# ##

# ##
# using Plots
# theme(:dao)
# scatter(sonora.mass, sonora.Keck_L′, marker_z=sonora.Teff, colorbar_title=raw"$\mathrm{T_{eff}}$", label="")
# xlabel!("mass")
# ylabel!("L′")
# ##
# scatter(sonora.logg, sonora.Keck_L′, marker_z=sonora.Teff, colorbar_title=raw"$\mathrm{T_{eff}}$", label="")
# xlabel!("log g")
# ylabel!("L′")

# ##

# ##
# itp = LinearInterpolation((sonora.logg, sonora.Teff), sonora.Keck_L′)

# plot(itp.(3.5))

# ##

# # ##
# # points = vcat(
# #     sonora.logg',
# #     log10.((sonora.Teff))'
# # )
# # samples = sonora.Keck_Ks
# # itp = ScatteredInterpolation.interpolate(ScatteredInterpolation.ThinPlate(), points, samples);
# # interpolated = evaluate(itp, [4.2; 1200])

# ##
# gs = 3:0.05:5.5
# ts = 200:5:2500
# heatmap(
#     gs, ts, sonora_flux_interp.(Ref(itp_L′), gs, ts')', colorbar_title="L′", label="",
#     background=:black,
#     foreground=:white,
#     fontfamily="",
#     color=:turbo
# )
# xlabel!("log g")
# ylabel!("T_eff")
# scatter!(sonora.logg, sonora.Teff, marker_z=sonora.Keck_L′, color=:turbo, colorbar_title="Ks", label="", markerstrokewidth=1, markerstrokecolor=:white)


# ##

# #=
# So, we have flux in log Jy.

# We ultimately have images in contrast. Wait? Or do we?
# We have the planet brightness already in each band. Should we put th images in terms of Jy??


# =#