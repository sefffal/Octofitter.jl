#=
Script to generate minimal test data subsets for documentation builds.

Run this script once (with the full data catalogs available) to create small
subset files that can be committed to the repository. This allows docs to build
without downloading large (multi-GB) data files.

Usage (from repository root):
    julia --project=@octo-tests test/data/generate_doc_subsets.jl

This will create:
- test/data/HGCA-test-subset.fits  (HGCA data for gaia_ids used in docs)
- test/data/HIRES-HD22049.txt      (HIRES RV data for HD22049)
=#

using Octofitter
using Octofitter: FITSIO, FITS, Table, @datadep_str
using OctofitterRadialVelocity

# Gaia IDs used in documentation examples
DOC_GAIA_IDS = [
    756291174721509376,    # pma.md, data-simulation.md, astrom-pma-rv.md (HD 91312)
    6166183842771027328,   # limits.md
    5164707970261890560,   # rv.md (epsilon Eridani)
]

function create_hgca_subset()
    println("Creating HGCA test subset...")

    # Load full HGCA catalog
    hgca_path = (datadep"HGCA_eDR3") * "/HGCA_vEDR3.fits"
    hgca_all = FITS(hgca_path, "r") do fits
        Table(fits[2])
    end

    # Find rows for our target gaia_ids
    indices = Int[]
    for gaia_id in DOC_GAIA_IDS
        idx = findfirst(==(gaia_id), hgca_all.gaia_source_id)
        if !isnothing(idx)
            push!(indices, idx)
            println("  Found gaia_id $gaia_id at row $idx")
        else
            @warn "gaia_id $gaia_id not found in HGCA catalog"
        end
    end

    # Extract subset
    subset = hgca_all[indices]

    # Write to FITS file
    output_path = joinpath(@__DIR__, "HGCA-test-subset.fits")
    FITS(output_path, "w") do f
        # Write data as binary table (FITSIO will create primary HDU automatically)
        colnames = String.(collect(propertynames(subset)))
        coldata = [getproperty(subset, col) for col in propertynames(subset)]
        write(f, colnames, coldata)
    end
    println("  Wrote $(length(indices)) rows to $output_path")
end

function create_hires_subset()
    println("Creating HIRES test subset...")

    # Get the HD22049 file from the full catalog
    target = "HD22049"
    catalog = datadep"HIRES_rvs"

    fname = OctofitterRadialVelocity.HIRES_search(target, catalog)

    # Copy to test data directory
    output_path = joinpath(@__DIR__, "HIRES-HD22049.txt")
    cp(fname, output_path, force=true)
    println("  Copied to $output_path")
end

function main()
    println("Generating documentation test data subsets...")
    println("=" ^ 50)

    create_hgca_subset()
    create_hires_subset()

    println("=" ^ 50)
    println("Done! Test data files created in $(dirname(@__FILE__))")
    println("\nThese files should be committed to the repository.")
end

main()
