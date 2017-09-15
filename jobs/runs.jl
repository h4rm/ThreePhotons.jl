using ThreePhotons

function check_git_status()
    status = readstring(`git ls-files --modified`)
    @assert length(status) == 0 "Please commit all changes before running server scripts."

    return readstring(`git rev-parse HEAD`)
end

function jobname(directory, number::Int64=0)
    main = replace(directory, "/", "_")
    suffix = number > 0 ? "_$(number)" : ""
    return "$main$suffix"
end

if contains(readstring(`hostname`), "owl")
    include("environment_owl.jl")
elseif contains(readstring(`hostname`), "hydra")
    include("environment_hydra.jl")
elseif contains(readstring(`hostname`), "gwdu103")
    include("environment_gwdg.jl")
else
    include("environment_local.jl")
end

include("run_generate_histograms.jl")
include("run_optimal.jl")
include("run_parallel_determination.jl")







# run_calculate_correlation_from_images("coliphage_symmetric", "$(ENV["DETERMINATION_DATA"])/exp_data/Coliphage_PR772/amo86615_194_PR772_single.h5", 24, 38, 26, 16, 42)

# run_calculate_correlation_from_images("coliphage_symmetric_N32", "$(ENV["DETERMINATION_DATA"])/exp_data/Coliphage_PR772/amo86615_194_PR772_single.h5", 24, 38, 26, 32, 42)
