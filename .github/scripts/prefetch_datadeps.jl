# Pre-fetch the network-backed DataDeps used by the test suite, with retries.
#
# Why this exists (see issue #128): several tests reach out over the network at
# test time -- most painfully `de440.bsp`, a 128 MiB SPK kernel from
# naif.jpl.nasa.gov. DataDeps uses a 30 s connect timeout, and when that timed
# out the test job *errored* and the whole matrix went red even though every
# assertion had passed ("0 failed, 2 errored").
#
# The depot cache (`julia-actions/cache`, which caches `scratchspaces` where
# DataDeps stores its downloads) means the network is normally not touched at
# all. This script covers the remaining cold-cache runs: it fetches each
# dependency up front and retries with backoff, so a single blip doesn't fail
# the job.
#
# Deliberately non-fatal. If a dependency still can't be fetched we emit a
# workflow warning and exit 0, leaving the test run to try again on demand --
# the goal is to remove a spurious failure mode, not to add a new one.

using Octofitter
using Octofitter.DataDeps: resolve

# Only the ones the test suite actually pulls over the network. The G23H
# dependencies are deliberately excluded: they live on a host that is not always
# reachable, and are not needed by the default test run.
const DEPS = ("DE440_Ephemeris", "HGCA_eDR3", "Hipparcos_IAD")

const MAX_ATTEMPTS = 3

failed = String[]
for name in DEPS
    for attempt in 1:MAX_ATTEMPTS
        try
            path = resolve(name, "", @__FILE__)
            println("ok: $name -> $path")
            break
        catch err
            println("failed to resolve $name (attempt $attempt/$MAX_ATTEMPTS): ",
                    sprint(showerror, err))
            if attempt == MAX_ATTEMPTS
                push!(failed, name)
            else
                backoff = 30 * attempt
                println("retrying in $(backoff)s...")
                sleep(backoff)
            end
        end
    end
end

if isempty(failed)
    println("All DataDeps are present.")
else
    # `::warning::` surfaces this in the workflow summary without failing the job.
    println("::warning::Could not pre-fetch DataDeps: ", join(failed, ", "),
            ". The test run will retry these on demand.")
end
