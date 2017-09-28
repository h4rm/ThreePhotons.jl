using ThreePhotons
using Compat

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
