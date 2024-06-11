#=

Correlations


Antoine:
I use only bootstrapping and group data by baseline and date: there are correlations arising from calibration, and they affect all data from the same baseline and date the same way (across spectral channels).

So we should use a correlation matrix between all baseline and date for all spectral channels.


M × B × Λ, where 
m = 1 . . . M is the number of individual measurements,
b = 1 . . . B is the number of baselines and
λ = 1 . . . Λ is the number of spectral channels.

From the coherent flux, we compute the squared visibility amplitudes
T3_mtλ = K · ∠VIS_mbλ

where
t = 1...T is the number of triangles,
K is a stack of M matrices k which encode how the four unique triangles can be formed from the six unique baselines of the VLTI, that is



FOR VIS2:
The correlation matrix is a block matrix consisting of B × B blocks, where each individual block is a Λ × Λ matrix

FOR T3:

C_T3 = [
    Y1 Y2 Y2 ... Y2
    Y2 Y1 Y2 ... Y2
    Y2 Y2 Y1 ... Y2
    Y2 Y2 Y2 ... Y2
    ...
    Y2 Y2 Y2 ... Y1
    ]

Y1 = [
    1 y ... y
    y 1 ... y
    ...    
    y y ... 1
]

Y2 = [
    ±1/3  ±y/3  ±y/3
    ±y/3  ±1/3  ±y/3 
    ...
    ±y/3  ±y/3  ±1/3 
]

For my implementation, I should construct a matrix with this structure analytically. 
ie a struct with a single parameter.
=#


using BlockArrays
function CT3(dat, y)
    # https://www.aanda.org/articles/aa/pdf/2020/12/aa38563-20.pdf

    # The correlation matrix is a block matrix consisting of T × T
    # where each individual block is a Λ×Λ matrix
    # λ = 1 . . . Λ is the number of spectral channels.
    # t = 1...T is the number of triangles,
    
    Λ = length(dat.eff_wave)
    T = size(dat.cps_data,1)
    L = T * Λ

    C = BlockedArray{typeof(y)}(zeros(typeof(y), L,L), fill(Λ,T), fill(Λ, T))

    # populate the block banded matrix
    for t_i in 1:T, t_j in 1:T
        if t_i == t_j
            # Blocks along the diagonal: Y1
            # block = fill(y, Λ, Λ)
            # block[diagind(block)] .= 1
            # C[Block(t_i,t_i)] = block
            block = view(C, Block(t_i,t_i))
            block .= y
            for λ in 1:Λ
                block[λ,λ] = 1
            end
        else
            # Blocks not along the diagonal: Y2
            thesign = -1
            if dat.index_cps3[t_i] == dat.index_cps3[t_j] ||
                dat.index_cps1[t_i] == dat.index_cps1[t_j] ||
                dat.index_cps1[t_i] == dat.index_cps2[t_j] ||
                dat.index_cps2[t_i] == dat.index_cps1[t_j] ||
                dat.index_cps2[t_i] == dat.index_cps2[t_j]
                thesign = +1
            end
            # block = fill(thesign*y/3, Λ, Λ)
            # block[diagind(block)] .= thesign*1/3
            # C[Block(t_i,t_j)] = block
            block = view(C, Block(t_i,t_j))
            block .= thesign*y/3
            for λ in 1:Λ
                block[λ,λ] = thesign/3
            end
        end
    end


    # The  CPs are computed with index_cps1 + index_cps2 - index_cps3
    # The sign is positive if the two triangles share a baseline in parallel
    # direction and negative if they share a baseline in anti-parallel direction.

    # I need to see if, for this block comparing one T3 against another T3 (computed from three baselines each),
    # do they share a baseline +/-?

    # Inidices into u and v give the baseline numbers
    # cp_inds give the baselines to subtract

    return C
end

# using AstroImages
# CT3(vis_like.table[1], 7.4e-2);
# @time C = CT3(vis_like.table[1], 7.4e-2);
# # @descend CT3(vis_like.table[1], 7.4e-2)

# # fig,ax,pl=heatmap(collect(reverse(C,dims=1)), colorrange=(-0.5,0.5), colormap=cgrad(:diverging_bwr_20_95_c54_n256,rev=true), axis=(;aspect=1,autolimitaspect=1))#cgrad(:seaborn_icefire_gradient,rev=true))
# # Colorbar(fig[1,2],pl)
# # fig