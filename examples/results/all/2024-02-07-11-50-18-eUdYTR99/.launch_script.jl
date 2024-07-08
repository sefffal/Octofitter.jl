using Serialization
using Pigeons
include(raw"/Users/thompsonw/.julia/dev/Octofitter/examples/distributed-model.jl")
Pigeons.mpi_active_ref[] = true

pt_arguments = 
    try
        Pigeons.deserialize_immutables!(raw"/Users/thompsonw/.julia/dev/Octofitter/examples/results/all/2024-02-07-11-50-18-eUdYTR99/immutables.jls")
        deserialize(raw"/Users/thompsonw/.julia/dev/Octofitter/examples/results/all/2024-02-07-11-50-18-eUdYTR99/.pt_argument.jls")
    catch e
        println("Hint: probably missing dependencies, use the dependencies argument in MPI() or ChildProcess()")
        rethrow(e)
    end

pt = PT(pt_arguments, exec_folder = raw"/Users/thompsonw/.julia/dev/Octofitter/examples/results/all/2024-02-07-11-50-18-eUdYTR99")
pigeons(pt)
