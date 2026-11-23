# Cut the small fixtures the test suite runs on out of the full products.
# The fixtures are committed; this script is how they were made.
#
#   julia --project=<env with OctofitterGaia and HDF5> make_fixtures.jl <mist grid .h5> <sed nodes .h5>
#
# mist_fixture.h5   [Fe/H] -0.25, 0, +0.25; log age 9.55–9.80; EEPs up to 808
#                   (the end of core helium burning)
# sed_fixture.h5    Teff 3000–6500 K, log g 4.0–5.5, [M/H] -0.5, 0, +0.5
#
# Both keep the layout of the full files, so the same loaders read them.

using HDF5

function cut_mist(src, dst)
    h5open(src) do f
        feh, age = read(f["feh"]), read(f["log_age"])
        js = findall(x -> -0.25 <= x <= 0.25, feh)
        ks = findall(x -> 9.549 <= x <= 9.801, age)
        ne = 809
        h5open(dst, "w") do g
            for k in keys(attributes(f))
                attributes(g)[k] = read(attributes(f)[k])
            end
            attributes(g)["fixture"] = "subset of the full grid: [Fe/H] $(feh[js]), log age $(age[ks]), EEP <= $(ne - 1)"
            g["feh"] = feh[js]
            g["log_age"] = age[ks]
            g["eep_lo"] = read(f["eep_lo"])[ks, js]
            g["eep_hi"] = min.(read(f["eep_hi"])[ks, js], ne)
            for name in ("log_mass", "log_teff", "logg", "log_L", "feh_surf", "MG", "MBP", "MRP", "phase", "filled")
                g[name, chunk=(ne, 1, 1), deflate=9] = read(f[name])[1:ne, ks, js]
            end
        end
    end
end

function cut_sed(src, dst)
    h5open(src) do f
        teff, logg, mh = read(f["teff"]), read(f["logg"]), read(f["mh"])
        is = findall(x -> 3000 <= x <= 6500, teff)
        js = findall(x -> 4.0 <= x <= 5.5, logg)
        ks = findall(x -> -0.5 <= x <= 0.5, mh)
        h5open(dst, "w") do g
            for k in keys(attributes(f))
                attributes(g)[k] = read(attributes(f)[k])
            end
            attributes(g)["fixture"] = "subset of the full node product"
            g["teff"] = teff[is]
            g["logg"] = logg[js]
            g["mh"] = mh[ks]
            g["edges_nm"] = read(f["edges_nm"])
            lf = read(f["log_flux"])[:, is, js, ks]
            g["log_flux", chunk=size(lf), deflate=9] = lf
            g["present"] = read(f["present"])[is, js, ks]
        end
    end
end

cut_mist(ARGS[1], joinpath(@__DIR__, "mist_fixture.h5"))
cut_sed(ARGS[2], joinpath(@__DIR__, "sed_fixture.h5"))
