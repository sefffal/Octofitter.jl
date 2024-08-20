#=
Functions for saving chains to disk and restoring them.

We use a FITS table so that we can store a table along with 
metadata in the header.

We convert unicode column names into ASCII latex expressions
and then convert back when reading them.
=#

using Latexify, MathTeXEngine

using FITSIO
function savechain(fname, chain::MCMCChains.Chains)
    # Need to massage data into the right format on the way out
    # data = permutedims(stack(row for row in Table(chain_octo)))
    tbl = Table(chain)
    FITS(fname, "w") do fits
        
        info = chain.info
        ks = collect(string.(keys(info)))
        vals = convert(Vector{Any}, collect(values(info)))
        h = nothing

        # Expand any vectors into multiple headers
        for k in ks
            i = findfirst(==(k),ks)
            if typeof(vals[i]) <: AbstractArray{<:Number}
                for ki in eachindex(vals[i])
                    push!(ks, "$(k)_$ki")
                    push!(vals, vals[i][ki])
                end
                vals[i] = "ARRAY"
            end
        end

        normalize(x::Number) = x
        normalize(x::AbstractString) = string(x)
        normalize(x::Nothing) = x
        normalize(x::Any) = filter(isascii, string(x)[1:min(255,end)])
        if length(ks) > 0
            h = FITSHeader(
                ks,
                normalize.(vals),
                ["chain-info" for _ in vals]
            )
        end

        display(h)

        write(fits,zeros(0), header=h)
        col_titles = ["iteration"; "chain"; string.(keys(chain))]

        write(
            fits,
            String.(map(title->all(isascii,title) ? title : latexify(title), col_titles)),
            [collect(getproperty(tbl,k)) for k in propertynames(tbl)];
            hdutype=FITSIO.TableHDU,
        )
    end
end


function coltitle_restorer(title::AbstractString)
    replacements = Pair{String,Char}[]
    for m in eachmatch(r"(\\[^\W_]+)[}_]", title)
        if !haskey(MathTeXEngine.command_definitions, m.captures[1])
            @error "Could not restore column" title m.captures[1] 
            error()
        end
        push!(
            replacements,
            m.captures[1] =>
                Char(MathTeXEngine.command_definitions[m.captures[1]][1].args[1])
        )
    end
    replace(title, replacements..., '$'=>"", '{'=>"", '}'=>"")
end

function loadchain(fname)
    return FITS(fname) do fits
        tbl = Table(fits[2])
        coltitles_serialized = propertynames(tbl)
        coltitles_deserialized = map(coltitle_restorer, string.(coltitles_serialized))
        out = map(coltitles_serialized, coltitles_deserialized) do colin, colout
            Symbol(colout) => getproperty(tbl, colin)
        end
        input =  Table(;out...)
        # Now we have the original table back with proper column names

        # Convert back into MCMCChains 
        chain_nums = input.chain
        unique_chain_nums = sort(unique(chain_nums))
        cols = setdiff(propertynames(input), (:iteration, :chain))
        T = typeof(sum(collect(first(input))))
        data = zeros(
            T,
            size(input, 1) ÷ length(unique_chain_nums),
            length(cols),
            length(unique_chain_nums)
        )
        for (j, key) in enumerate(cols)
            col = getproperty(input, key)
            for i in eachindex(input)
                data[
                    input.iteration[i],
                    j,
                    input.chain[i]
                ] = col[i]
            end
        end
        chn =  MCMCChains.Chains(data, cols)
        
        # And restore metadata
        h = read_header(fits[1])
        meta = Dict{Symbol,Any}()
        for k in keys(h)
            if get_comment(h,k) == "chain-info"
                if h[k] != "ARRAY"
                    meta[Symbol(lowercase(k))] = h[k]
                else
                    arrdat = []
                    i = 0
                    while true
                        i+=1
                        ki = "$(k)_$i"
                        if haskey(h,ki)
                            push!(arrdat, h[ki])
                        else
                            break
                        end
                    end
                    meta[Symbol(lowercase(k))] = identity.(arrdat)
                end
            end
        end
        chn = MCMCChains.setinfo(chn, NamedTuple(meta))
        return chn
    end
end
