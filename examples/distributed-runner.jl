include("distributed-model.jl")
##
using Pigeons
pt = pigeons(target=model ,n_rounds=8, record=[traces; round_trip; record_default()])
results = Chains(model, pt)
##
octocorner(model, results,small=true)
##
# to run
pt = pigeons(
    target = Pigeons.LazyTarget(MyTargetFlag()),
    record = [traces; round_trip; record_default()],
    on = ChildProcess(
        n_local_mpi_processes = 4,
        dependencies = ["distributed-model.jl"]
    )
)