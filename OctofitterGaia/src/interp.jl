# ---------------------------------------------------
# Interpolation on tabulated axes
# ---------------------------------------------------

# Catmull–Rom stencil on a sorted, possibly non-uniform axis: four node
# indices and their weights, interpolating and C¹, and the fraction `t` of
# the way from the second node to the third. Tangents are the three-point
# parabola derivatives, one-sided at the ends of the axis, where the unused
# index is clamped and carries zero weight. `nothing` outside the axis.
@inline function _stencil(nodes::Vector{Float64}, x)
    n = length(nodes)
    (x >= nodes[1] && x <= nodes[n]) || return nothing
    i = min(searchsortedlast(nodes, x), n - 1)
    @inbounds x1, x2 = nodes[i], nodes[i+1]
    d1 = x2 - x1
    t = (x - x1) / d1
    h00 = (2t - 3) * t^2 + 1
    h10 = ((t - 2) * t + 1) * t
    h01 = (3 - 2t) * t^2
    h11 = (t - 1) * t^2
    w0, w1, w2, w3 = zero(t), h00, h01, zero(t)
    if i > 1
        d0 = x1 - @inbounds nodes[i-1]
        w0 -= h10 * d1 * d1 / (d0 * (d0 + d1))
        w1 += h10 * (d1 - d0) / d0
        w2 += h10 * d0 / (d0 + d1)
    else
        w1 -= h10
        w2 += h10
    end
    if i < n - 1
        d2 = @inbounds(nodes[i+2]) - x2
        w1 -= h11 * d2 / (d1 + d2)
        w2 += h11 * (d2 - d1) / d2
        w3 += h11 * d1 * d1 / (d2 * (d1 + d2))
    else
        w1 -= h11
        w2 += h11
    end
    return (max(i - 1, 1), i, i + 1, min(i + 2, n)), (w0, w1, w2, w3), t
end
