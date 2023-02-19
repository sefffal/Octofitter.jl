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
    # @show eltype(data)
    tbl = Table(chain)
    FITS(fname, "w") do fits
        
        info = chain.info
        ks = collect(string.(keys(info)))
        vals = collect(values(info))
        h = nothing
        if length(ks) > 0
            h = FITSHeader(
                ks,
                vals,
                ["chain-info" for _ in values(info)]
            )
        end

        write(fits,zeros(0,0), header=h)
        col_titles = ["iteration"; "chain"; string.(keys(chain))]

        write(
            fits,
            String.(map(title->all(isascii,title) ? title : latexify(title), col_titles)),
            [collect(getproperty(tbl,k)) for k in propertynames(tbl)];
            hdutype=FITSIO.TableHDU,
            # header=nothing
        )
    end
end


function coltitle_restorer(title::AbstractString)
    replacements = Pair{String,Char}[]
    for m in eachmatch(r"(\\\w+)}", title)
        push!(
            replacements,
            m.captures[1] =>
                Char(MathTeXEngine.command_to_canonical[m.captures[1]])
        )
    end
    title, replacements
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
                meta[Symbol(lowercase(k))] = h[k]
            end
        end
        chn = MCMCChains.setinfo(chn, NamedTuple(meta))
        return chn
    end
end
