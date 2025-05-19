

function getrv(dr)
    if !hasproperty(dr, :radial_velocity)
        0.
    elseif isnothing(dr.radial_velocity)
        0.
    else
        dr.radial_velocity
    end
end




