###############################################################
#                     SGE-Scheduling
###############################################################

if contains(readstring(`hostname`), "owl")
    include("environment_owl.jl")
elseif contains(readstring(`hostname`), "hydra")
    include("environment_hydra.jl")
elseif contains(readstring(`hostname`), "gwdu103")
    include("environment_gwdg.jl")
else
    include("environment_local.jl")
end

include("../src/utilities.jl")

function jobname(directory, number::Int64=0)
    main = replace(directory, "/", "_")
    suffix = number > 0 ? "_$(number)" : ""
    return "$main$suffix"
end

function check_git_status()
    status = readstring(`git ls-files --modified`)
    @assert length(status) == 0 "Please commit all changes before running server scripts."

    return readstring(`git rev-parse HEAD`)
end

function run_determination(dir::String; histograms::String="", initial_stepsize::Float64=Float64(pi), K::Integer=8, L::Integer=8, optimizer::String="rotate_hierarchical", initial_temperature_factor::Float64=1.0, temperature_decay::Float64=0.99, N::Integer=32, range=1000:1019, fresh::Bool=false, gpu::Bool=true, Ncores::Integer=8, successive_jobs::Integer=1, measure="Bayes", postprocess::Bool=true, stepsizefactor::Float64=1.02, KMAX::Int64=35, rmax::Float64=35.0, force_repostprocess::Bool=false, run_denoise::Bool=false, architecture::String="ivy-bridge|sandy-bridge|haswell|broadwell|skylake", hours::Int64=48, sigma::Float64=0.0, reference_pdb_path::String="", lambda::Float64=0.0)

    julia_script = """
    using ThreePhotons

    #Let's continue
    if isfile("state.dat") && isfile("params.dat")

        state = deserializeFromFile("state.dat")
        params = deserializeFromFile("params.dat")

        if state["state"] == "running"
            params,state = rotation_search(params, state)
            postprocess_run(params, state, "$(reference_pdb_path)", true, 35, $sigma)
        elseif state["state"] == "finished_structure" && $postprocess
            postprocess_run(params, state, "$(reference_pdb_path)", true, 35, $sigma)
        elseif $force_repostprocess
            #TODO: Temporary fix for wrong results in noise runs
            params["rmax"] = $rmax
            state["intensity"].rmax = qmax(params["KMAX"], params["rmax"])
            postprocess_run(params, state, "$(reference_pdb_path)", true, 35, $sigma)
        end

        #Or start a completely new run
    else
        params,state = rotation_search(Dict( "reference_pdb_path"=>"$(reference_pdb_path)", "stepsizefactor"=>$stepsizefactor, "initial_stepsize" => $initial_stepsize, "L"=>$L, "K" =>$K, "N"=>$N, "histograms"=>"$(histograms)", "optimizer"=>$optimizer, "initial_temperature_factor"=>$initial_temperature_factor, "measure"=>"$measure", "temperature_decay"=>$temperature_decay, "LMAX"=>25, "KMAX"=>$KMAX, "rmax"=>$rmax, "lambda"=>$lambda))
        if $postprocess
            postprocess_run(params, state, "$(reference_pdb_path)", true, 35, $sigma)
        end
    end
    """

    for n in range
        launch_job("$dir/$n", Ncores, gpu, julia_script, successive_jobs; architecture=architecture, hours=hours, fresh=fresh)
    end
end

function histogram_name(prefix::String, ppi::Int64, N::Int64, KMAX::Int64, rmax::Float64, max_pictures::Int64, suffix::String, gamma::Float64=0.0, sigma::Float64=1.0)
    return "$(prefix)$(ppi)p_N$(N)_K$(KMAX)_R$(rmax)_P$(max_pictures)$(gamma > 0.0 ? "_G$(gamma)_S$(sigma)" : "")$(suffix)"
end

"""Starts a cluster job for synthetic correlation generation"""
function generate_histograms(; max_triplets::Integer=Integer(0), max_pictures::Integer=Integer(0), N::Integer=32, photons_per_image::Integer=500, incident_photon_variance::Integer = 0, gamma::Float64=0.0, sigma::Float64=1.0, noise_photons::Int64=0, Ncores::Integer=8, batchsize::Integer = Integer(1e4), successive_jobs::Integer=1, prefix::String="correlations_", suffix::String="", use_cube::Bool=true, qcut_ratio::Float64=1.0, K::Integer=35, rmax::Float64=35.0, histogram_method="histogramCorrelationsInPicture_alltoall", structure_pdb_path::String="", number_particles::Int64=1, lambda::Float64=1.0)
    name = ""
    if max_triplets > 0
        name = "$(prefix)_N$(N)_K$(K)_R$(rmax)_T$(max_triplets)$(gamma > 0.0 ? "_G$(gamma)_S$(sigma)" : "")$(suffix)"
    elseif max_pictures > 0
        # name = "$(prefix)_N$(N)_K$(K)_R$(rmax)_P$(max_pictures)$(gamma > 0.0 ? "_G$(gamma)_S$(sigma)" : "")$(suffix)"
        name = histogram_name(prefix, photons_per_image, N, K, rmax, max_pictures, suffix, gamma, sigma)
    end
    number_incident_photons = calculate_incident_photons(photons_per_image)
    julia_script = """
    using ThreePhotons

    volume = 0.0

    if $(use_cube)
        # _,_,volume = createCubicStructure("$(structure_pdb_path)", 4*$K+1, 2.0*$rmax)
        volume = loadCube("$(ENV["DETERMINATION_DATA"])/histograms/synthetic_crambin/intensityCube_high.mrc")
    else
        _,_,intensity = createSphericalHarmonicsStructure("$(structure_pdb_path)", 35, $K, $rmax)
        volume = getSurfaceVolume(intensity)
        volume.radial_interp = false
    end

    generateHistogram(volume; qcut=$(qcut_ratio)*volume.rmax, K=$K, N=$N, max_triplets=$max_triplets, max_pictures=$max_pictures, number_incident_photons=$number_incident_photons, incident_photon_variance=$incident_photon_variance, numprocesses=$(Ncores), file="histo.dat", noise=GaussianNoise($gamma, $sigma, $(noise_photons)), batchsize = $batchsize, histogramMethod=$histogram_method, number_particles=$(number_particles), lambda=$(lambda))
    """

    launch_job("data_generation/$name", Ncores, false, julia_script, successive_jobs, architecture="haswell|broadwell|skylake")
end

"""Calculates the resolution of corresponding (dq,K,L) combination"""
function run_optimal(K2::Int64, K3::Int64, L3::Int64)
    name = "optimal/optimal_K2_$(K2)_K3_$(K3)_L3_$(L3)"

    julia_script = """
    using ThreePhotons
    deltar = 1.0 #2 Angstrom resolution maximum
    qmax2 = pi/deltar
    K2 = $K2
    rmax = K2*deltar
    deltaq = qmax2 / K2

    K3 = $K3
    L3 = $L3
    _,_,intensity = createSphericalHarmonicsStructure("$(ENV["THREEPHOTONS_PATH"])/structures/crambin.pdb", 25, K2, rmax)
    densityCube,fourierCube,intensityCube = createCubicStructure("$(ENV["THREEPHOTONS_PATH"])/structures/crambin.pdb", 2*K2+1, rmax)
    calculate_optimal(intensity,densityCube,fourierCube,intensityCube, K3, L3, ".", 16)

    #Now check the same with retrieving the structure
    checkRotation_basis = complexBasis(L3,32,25)
    c2hist_full = twoPhotons(intensity, checkRotation_basis, intensity.KMAX, true, false)
    start = retrieveSolution(c2hist_full, checkRotation_basis.L, intensity.LMAX, intensity.KMAX, intensity.rmax)
    energyfunc = (volume) -> 0.0
    reference_intensity = deleteTerms(intensity, K3, checkRotation_basis.L)
    rightCoeff,finalFSC = checkRotationSearch(start, reference_intensity, K3, checkRotation_basis.L, checkRotation_basis, save_structures=false, energy = energyfunc, iterations=4.0e4, reduce_stepsize=3000)
    saveCube(rightCoeff, "./retrievedIntensity.mrc")

    #Phase retrieved structure
    res = calculateSC(deleteTerms(rightCoeff, K3,checkRotation_basis.L), densityCube, fourierCube, intensityCube, 16)
    saveCube(res[1], "./retrievedDensity.mrc")
    serializeToFile("./retrieved.dat", res)
    """
    launch_job(name, 8, false, julia_script, 1)
end

"""Distributes the calculation of correlations among many jobs"""
function run_calculate_correlation_from_images(particle_name::String, images_path::String, number_images::Int64, K2::Int64, K3::Int64, N::Int64, number_runs::Int64; Ncores::Int64=8)
    for n in 1:number_runs
        julia_script = """
        using HDF5
        using ThreePhotons
        using Images

        K2 = $K2
        K3 = $K3
        N = $N

        file = h5open("$(images_path)", "r")
        photonConverter = read(file["photonConverter"])
        resized_image_list = [ Images.imresize(convert(Images.Image,convert(Array{Float64},photonConverter["pnccdBack"]["photonCount"][:,:,i])), (2*K2, 2*K2)).data for i=$((n-1)*number_images+1):$(n*number_images)]
        calculate_correlations_in_image(resized_image_list, K2, K3, N)
        """
        launch_job("exp_data/parts/$(particle_name)_$(n)", Ncores, false, julia_script, 1)#, memory="$(Ncores*1.5)G")
    end
end

# run_calculate_correlation_from_images("coliphage_symmetric", "$(ENV["DETERMINATION_DATA"])/exp_data/Coliphage_PR772/amo86615_194_PR772_single.h5", 24, 38, 26, 16, 42)
