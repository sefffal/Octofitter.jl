
# function sonora()

# end
using DelimitedFiles, NamedTupleTools, Interpolations
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

##
sonora = load_table()

##
using Plots
theme(:dao)
scatter(sonora.mass, sonora.Keck_L′, marker_z=sonora.Teff, colorbar_title=raw"$\mathrm{T_{eff}}$", label="")
xlabel!("mass")
ylabel!("L′")
##
scatter(sonora.logg, sonora.Keck_L′, marker_z=sonora.Teff, colorbar_title=raw"$\mathrm{T_{eff}}$", label="")
xlabel!("log g")
ylabel!("L′")

##

##
itp = LinearInterpolation((sonora.logg, sonora.Teff), sonora.Keck_L′)

plot(itp.(3.5))

##
using ScatteredInterpolation
samples = [0.0; 0.5; 0.5; 0.5; 1.0];
points = [0.0 0.0; 0.0 1.0; 1.0 0.0; 0.5 0.5; 1.0 1.0]';
itp = ScatteredInterpolation.interpolate(Multiquadratic(), points, samples);
interpolated = evaluate(itp, [0.6; 0.6])

##
points = vcat(
    sonora.logg',
    log10.((sonora.Teff))'
)
itp_Ks = ScatteredInterpolation.interpolate(ScatteredInterpolation.ThinPlate(), points, sonora.Keck_Ks);
itp_L′ = ScatteredInterpolation.interpolate(ScatteredInterpolation.ThinPlate(), points, sonora.Keck_L′);
function sonora_flux_interp(itp, logg, Teff)
    return only(evaluate(itp, [logg; log10(Teff)]))
end

# ##
# points = vcat(
#     sonora.logg',
#     log10.((sonora.Teff))'
# )
# samples = sonora.Keck_Ks
# itp = ScatteredInterpolation.interpolate(ScatteredInterpolation.ThinPlate(), points, samples);
# interpolated = evaluate(itp, [4.2; 1200])

##
gs = 3:0.05:5.5
ts = 200:5:2500
heatmap(
    gs, ts, sonora_flux_interp.(Ref(itp_L′), gs, ts')', colorbar_title="L′", label="",
    background=:black,
    foreground=:white,
    fontfamily="",
    color=:turbo
)
xlabel!("log g")
ylabel!("T_eff")
scatter!(sonora.logg, sonora.Teff, marker_z=sonora.Keck_L′, color=:turbo, colorbar_title="Ks", label="", markerstrokewidth=1, markerstrokecolor=:white)


##

#=
So, we have flux in log Jy.

We ultimately have images in contrast. Wait? Or do we?
We have the planet brightness already in each band. Should we put th images in terms of Jy??


=#