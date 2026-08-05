# Fast pre-flight for the manual, to be run *before* `docs/make.jl`.
#
# A full Documenter build runs real MCMC on ~20 tutorial pages and takes
# hours. The overwhelming majority of what a large refactor breaks in the docs
# is not the science — it is a renamed symbol inside a code fence, a `@ref`
# pointing at a name that no longer exists, or a page added to `docs/src/`
# and never listed in `make.jl`. All of that is decidable statically, in a few
# seconds. This script decides it.
#
#     julia --project=docs docs/check-examples.jl
#     julia --project=docs docs/check-examples.jl --names   # also resolve @refs
#
# What it checks:
#
#   1. Every fenced code block on every page **parses** (`Meta.parse` over the
#      whole block, so an unbalanced paren or a stray `end` is caught). Blocks
#      are not executed — this says nothing about whether they *run*.
#   2. `make.jl`'s page list and `docs/src/*.md` agree: every listed page
#      exists (a missing one is a hard Documenter error), and every unlisted
#      page is reported (with `pagesonly=true` it is silently not built, which
#      is fine deliberately and a bug accidentally).
#   3. Every `@id` anchor is unique, and every `[…](@ref anchor)` names one
#      that exists.
#   4. Every `[`name`](@ref)` and every `@docs` entry names something that
#      actually resolves in the loaded packages (only with `--names`, which
#      needs the packages importable; without it the name checks are skipped
#      and everything else still runs).
#   5. Links to other pages (`](something.md)`) point at files that exist and
#      are in the nav.
#
# Exit code is non-zero if anything in 1–3 or 5 fails, so it works in CI.

const DOCS = @__DIR__
const SRC = joinpath(DOCS, "src")
const CHECK_NAMES = "--names" in ARGS

failures = String[]
warnings = String[]
fail!(msg) = push!(failures, msg)
warn!(msg) = push!(warnings, msg)

pages() = sort([replace(relpath(joinpath(root, f), SRC), '\\' => '/')
                for (root, _, files) in walkdir(SRC) for f in files if endswith(f, ".md")])

# ---------------------------------------------------------------------------
# make.jl's page list
# ---------------------------------------------------------------------------

"""Every `"…" => "page.md"` target in make.jl, in order."""
function listed_pages()
    txt = read(joinpath(DOCS, "make.jl"), String)
    # Only the `pages = [...]` argument, so a stray string elsewhere in the
    # file cannot be mistaken for a page.
    i = findfirst("pages = [", txt)
    isnothing(i) && (i = findfirst("pages=[", txt))
    isnothing(i) && error("could not find the `pages = [` list in make.jl")
    depth = 0
    j = last(i)
    for k in last(i):lastindex(txt)
        c = txt[k]
        c == '[' && (depth += 1)
        if c == ']'
            depth -= 1
            depth == 0 && (j = k; break)
        end
    end
    block = txt[first(i):j]
    return [m.captures[1] for m in eachmatch(r"=>\s*\"([^\"]+\.md)\"", block)]
end

# ---------------------------------------------------------------------------
# Fenced code blocks
# ---------------------------------------------------------------------------

struct Block
    page::String
    line::Int
    lang::String
    code::String
end

function codeblocks(page)
    out = Block[]
    lines = readlines(joinpath(SRC, page))
    i = 1
    while i <= length(lines)
        m = match(r"^\s*```+\s*(\S.*)?$", lines[i])
        if !isnothing(m) && !isnothing(m.captures[1])
            lang = strip(m.captures[1])
            start = i
            i += 1
            buf = String[]
            while i <= length(lines) && isnothing(match(r"^\s*```+\s*$", lines[i]))
                push!(buf, lines[i])
                i += 1
            end
            push!(out, Block(page, start, lang, join(buf, "\n")))
        end
        i += 1
    end
    return out
end

"""
Should this fence be parsed as Julia? `@example`, `@repl`, `@setup`,
`@eval` and plain `julia` are Julia; `@meta`, `@docs`, `@autodocs`,
`@raw`, `@index`, `@contents` and non-Julia languages are not.
"""
function isjulia(lang)
    l = lowercase(lang)
    startswith(l, "@example") && return true
    startswith(l, "@repl") && return true
    startswith(l, "@setup") && return true
    startswith(l, "@eval") && return true
    (l == "julia" || startswith(l, "julia ")) && return true
    return false
end

"""
`@repl` blocks are a sequence of independent top-level inputs, and `# hide`
suffixes are Documenter's, not Julia's — strip them before parsing.
"""
function parsecheck(b::Block)
    code = join([replace(l, r"#\s*hide\s*$" => "") for l in split(b.code, "\n")], "\n")
    isempty(strip(code)) && return nothing
    if startswith(lowercase(b.lang), "@repl")
        # Parse statement by statement: a REPL block may legitimately contain
        # several complete expressions with no separator.
        pos = 1
        while pos <= lastindex(code)
            ex, pos = try
                Meta.parse(code, pos; greedy=true, raise=true)
            catch err
                return sprint(showerror, err)
            end
            isnothing(ex) && break
        end
        return nothing
    end
    ex = try
        Meta.parseall(code)
    catch err
        return _oneline(sprint(showerror, err))
    end
    # `parseall` returns an expression containing `:error`/`:incomplete` heads
    # rather than throwing.
    for a in ex.args
        if a isa Expr && (a.head === :error || a.head === :incomplete)
            return _oneline(a.args[1] isa Exception ? sprint(showerror, a.args[1]) : string(a))
        end
    end
    return nothing
end

"""The first line of a ParseError that actually says what is wrong."""
function _oneline(s)
    for l in split(s, '\n')
        t = strip(l)
        (isempty(t) || t == "ParseError:" || startswith(t, "# Error @")) && continue
        return t
    end
    return first(split(s, '\n'))
end

"""
Is this block deliberately illustrative rather than runnable? Plain ```julia
fences in the manual routinely show REPL/Pkg prompts and elide bodies with
`...`; those never parsed and are not meant to. Only `@example`-family blocks,
which Documenter actually executes, are held to parsing.
"""
function is_illustrative(code)
    for l in split(code, '\n')
        t = strip(l)
        (startswith(t, "pkg>") || startswith(t, "julia>") || startswith(t, "shell>") ||
         startswith(t, "(@v") || occursin(r"(^|[\s(,;=])\.\.\.([\s)\]},;]|$)", t)) && return true
    end
    return false
end

# ---------------------------------------------------------------------------
# Anchors and cross-references
# ---------------------------------------------------------------------------

anchors_in(page) = [m.captures[1] for m in
    eachmatch(r"\(@id\s+([^)\s]+)\s*\)", read(joinpath(SRC, page), String))]

refs_in(page) = [(m.captures[1], m.captures[2]) for m in
    eachmatch(r"\[([^\]]*)\]\(@ref\s*([^)\s]*)\s*\)", read(joinpath(SRC, page), String))]

filelinks_in(page) = [m.captures[1] for m in
    eachmatch(r"\]\(([^)\s#]+\.md)(?:#[^)]*)?\)", read(joinpath(SRC, page), String))]

"""Names inside `@docs` / `@autodocs` fences, one per line."""
function docs_entries(page)
    out = Tuple{Int,String}[]
    for b in codeblocks(page)
        startswith(lowercase(b.lang), "@docs") || continue
        for (k, l) in enumerate(split(b.code, "\n"))
            s = strip(l)
            (isempty(s) || startswith(s, "#")) && continue
            push!(out, (b.line + k, s))
        end
    end
    return out
end

# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------

all_pages = pages()
listed = listed_pages()

println("Pages in docs/src: $(length(all_pages));  listed in make.jl: $(length(listed))")

for p in listed
    isfile(joinpath(SRC, p)) ||
        fail!("make.jl lists `$p`, which does not exist in docs/src/ (hard Documenter error).")
end
for p in all_pages
    p in listed || warn!("`$p` exists but is not in make.jl's page list; with `pagesonly=true` it will not be built.")
end
if length(listed) != length(unique(listed))
    for p in unique(listed)
        count(==(p), listed) > 1 && fail!("make.jl lists `$p` more than once.")
    end
end

# 1. parse every Julia code fence
nblocks = Ref(0)
for p in all_pages
    p in listed || continue
    for b in codeblocks(p)
        isjulia(b.lang) || continue
        nblocks[] += 1
        err = parsecheck(b)
        isnothing(err) && continue
        msg = "$(b.page):$(b.line)  ```$(b.lang) does not parse: $err"
        # Documenter executes the `@example` family, so a parse error there is
        # a build failure. A plain ```julia fence is display-only; hold it to
        # parsing unless it is obviously illustrative (`...`, a Pkg prompt).
        if lowercase(b.lang) == "julia" || startswith(lowercase(b.lang), "julia ")
            is_illustrative(b.code) || warn!(msg)
        else
            fail!(msg)
        end
    end
end
println("Parsed $(nblocks[]) Julia code blocks across the listed pages.")

# 3. anchors
anchor_owner = Dict{String,String}()
for p in all_pages
    for a in anchors_in(p)
        if haskey(anchor_owner, a)
            p in listed && anchor_owner[a] in listed &&
                fail!("duplicate `@id $a` in `$p` and `$(anchor_owner[a])`.")
        else
            anchor_owner[a] = p
        end
    end
end
for p in listed
    for (text, target) in refs_in(p)
        isempty(target) && continue           # name ref, handled below
        startswith(target, "@") && continue
        haskey(anchor_owner, target) ||
            fail!("$p: `[$text](@ref $target)` — no page defines `(@id $target)`.")
        get(anchor_owner, target, "") in listed ||
            fail!("$p: `[$text](@ref $target)` resolves to `$(anchor_owner[target])`, which is not built.")
    end
end
println("Checked $(length(anchor_owner)) `@id` anchors.")

# 5. inter-page links
for p in listed
    for l in filelinks_in(p)
        tgt = normpath(joinpath(dirname(p), l))
        tgt = replace(tgt, '\\' => '/')
        if !isfile(joinpath(SRC, tgt))
            fail!("$p: link to `$l` — no such file in docs/src/.")
        elseif !(tgt in listed)
            warn!("$p: links to `$l`, which make.jl does not build.")
        end
    end
end

# 4. names — only when the packages can be loaded
if CHECK_NAMES
    @eval using Octofitter, OctofitterRadialVelocity, PlanetOrbits
    mods = Module[Octofitter, OctofitterRadialVelocity, PlanetOrbits]
    for m in (:OctofitterImages, :OctofitterInterferometry)
        try
            @eval using $m
            push!(mods, getfield(Main, m))
        catch
            warn!("could not load $m; its `@docs` entries are unchecked.")
        end
    end

    """Does `entry` (possibly `Mod.name`, `@macro`, or `f(::T)`) resolve in any
    loaded module? Macros are checked by binding, never evaluated — `Core.eval`
    on a `:macrocall` would expand it."""
    function resolves(entry)
        e = strip(split(entry, '(')[1])
        isempty(e) && return false
        if startswith(e, "@")
            sym = Symbol(e)
            return any(m -> isdefined(m, sym), mods)
        end
        ex = try Meta.parse(e) catch; return false end
        for m in mods
            try
                Core.eval(m, ex)
                return true
            catch
            end
        end
        return false
    end

    for p in listed
        for (ln, e) in docs_entries(p)
            resolves(e) ||
                fail!("$p:$ln  `@docs` entry `$e` does not resolve — Documenter fails the build on this.")
        end
        for (text, target) in refs_in(p)
            isempty(target) || continue
            nm = strip(text, ['`', ' ', '[', ']'])
            (isempty(nm) || occursin(' ', nm)) && continue   # prose ref, resolved by anchor text
            resolves(nm) || warn!("$p: `[$text](@ref)` — `$nm` does not resolve as a name.")
        end
    end
    println("Checked `@docs` entries and name `@ref`s against the loaded packages.")
else
    println("Name resolution skipped (pass --names to enable).")
end

# ---------------------------------------------------------------------------

println()
if !isempty(warnings)
    println("WARNINGS ($(length(warnings))):")
    for w in warnings
        println("  - ", w)
    end
    println()
end
if isempty(failures)
    println("OK — no build-breaking problems found.")
else
    println("FAILURES ($(length(failures))):")
    for f in failures
        println("  - ", f)
    end
    exit(1)
end
