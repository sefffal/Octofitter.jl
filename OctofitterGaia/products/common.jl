# Shared by the node-product builders.

# Fill the missing nodes of `a[:, i, j, k]`, in place: two-point linear
# extrapolation along axis 3 (log g) where possible, then along axis 2
# (Teff), repeated until nothing is missing.
function fill_missing!(a, have)
    have = copy(have)
    nT, ng, nm = size(have)
    while !all(have)
        progress = false
        for axis in (3, 2)
            snapshot = copy(have)
            for k in 1:nm, j in 1:ng, i in 1:nT
                snapshot[i, j, k] && continue
                step(d) = axis == 3 ? (i, j + d, k) : (i + d, j, k)
                inb(I) = 1 <= I[1] <= nT && 1 <= I[2] <= ng
                for dir in (1, -1)
                    p1, p2 = step(dir), step(2dir)
                    (inb(p1) && inb(p2) && snapshot[p1...] && snapshot[p2...]) || continue
                    @views a[:, i, j, k] .= 2 .* a[:, p1...] .- a[:, p2...]
                    have[i, j, k] = true
                    progress = true
                    break
                end
            end
            progress && break
        end
        progress || error("could not fill every missing node")
    end
    return a
end
