# Run one manual page's `@example` blocks locally, the way Documenter would.
#
#     julia --project=docs docs/run-page.jl cross-validation.md
#     julia --project=docs docs/run-page.jl rv.md extract-phot-astrom.md
#
# A full `docs/make.jl` takes hours because it runs real MCMC on every tutorial.
# When you have broken (or are fixing) *one* page, this runs just that page:
# same block semantics as Documenter — `@example <name>` blocks with the same
# name share a module, `# hide` lines still execute, ```julia fences do not run —
# so a block that passes here passes in the build.
#
# Exits non-zero if any block throws, printing the block and the exception in
# the same shape the Documenter log uses.

const SRC = joinpath(@__DIR__, "src")

isempty(ARGS) && error("usage: julia --project=docs docs/run-page.jl <page.md> [...]")

"""
Fenced blocks of a page, as (language-info, code) pairs, in order.
"""
function fenced_blocks(md::AbstractString)
    blocks = Tuple{String,String}[]
    lines = split(md, '\n')
    i = 1
    while i <= length(lines)
        if startswith(lines[i], "```")
            info = strip(lines[i][4:end])
            j = i + 1
            while j <= length(lines) && !startswith(lines[j], "```")
                j += 1
            end
            push!(blocks, (String(info), join(lines[i+1:min(j - 1, length(lines))], '\n')))
            i = j + 1
        else
            i += 1
        end
    end
    return blocks
end

failures = Tuple{String,String,Any}[]

for page in ARGS
    path = isabspath(page) ? page : joinpath(SRC, page)
    isfile(path) || error("no such page: $path")
    name = basename(path)
    println("\n", "="^78, "\n== ", name, "\n", "="^78)

    mods = Dict{String,Module}()
    nrun = 0
    for (info, code) in fenced_blocks(read(path, String))
        # Documenter executes `@example` and `@setup`; `@repl` is line-by-line
        # but the same module rules apply. Everything else (```julia, ```) is
        # displayed only.
        kind = split(info)[1:min(1, length(split(info)))]
        isempty(kind) && continue
        k = first(kind)
        k in ("@example", "@setup", "@repl") || continue
        blockname = length(split(info)) > 1 ? split(info)[2] : "__anon__$(nrun)"
        mod = get!(mods, blockname) do
            m = Module(Symbol("__page__", replace(name, r"\W" => "_"), "__", blockname))
            Core.eval(m, :(eval(x) = Core.eval($m, x)))
            m
        end
        nrun += 1
        println("\n-- block $nrun ($info) --")
        try
            Core.eval(mod, Meta.parseall(code))
        catch err
            println(stderr, "\n!! FAILED: ", name, " block $nrun ($info)")
            println(stderr, code)
            println(stderr, "  exception = ")
            showerror(stderr, err, catch_backtrace())
            println(stderr)
            push!(failures, (name, code, err))
            break   # Documenter keeps going, but later blocks depend on this one
        end
    end
    println("\n== $name: ran $nrun block(s)")
end

println("\n", "="^78)
if isempty(failures)
    println("OK — every executed block ran.")
else
    println("FAILED — $(length(failures)) block(s):")
    for (page, code, err) in failures
        println("  $page: ", first(split(strip(code), '\n')), " …  → ", sprint(showerror, err))
    end
    exit(1)
end
