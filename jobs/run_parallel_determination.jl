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

function run_set(image_list::Array{Int64}, KMAX::Int64=38, N::Int64=32, L::Int64=18, K::Int64=26, temperature_decay::Float64=0.99998, ppi::Int64=10, include_infinite::Bool=true)
    histograms_finite = Dict( "P$(img)" => histogram_name("parallel/data_generation/SH_", ppi, N, KMAX, K, float(KMAX), img, "") for img in image_list)
    histograms_infinite = include_infinite ? Dict("L20_inf" => "expdata/correlations_N$(N)_K$(KMAX)_L20_inf.dat") : Dict()
    histogram_list = merge(histograms_finite, histograms_infinite)

    for img in keys(histogram_list)
        run_determination("paper_res_vs_pictures_$(ppi)p_KMAX$(KMAX)_N$(N)_K$(K)_L$(L)_$(temperature_decay)/$img", histograms=histogram_list[img], initial_stepsize=pi/4.0, K=K, L=L, KMAX=KMAX, rmax=float(KMAX), optimizer="rotate_all_at_once", initial_temperature_factor=0.1, temperature_decay=temperature_decay, N=N, successive_jobs=3, measure="Bayes", range=1000:1019, postprocess=true, gpu=true, Ncores=20, stepsizefactor=1.01)
    end
end

# run_set(calculate_images_ppi(10), 38, 32, 18, 26, 0.99998, 10)
# run_set(calculate_images_ppi(25), 38, 32, 18, 26, 0.99998, 25, false)
# run_set(calculate_images_ppi(50), 38, 32, 18, 26, 0.99998, 50, false)
# run_set(calculate_images_ppi(100), 38, 32, 18, 26, 0.99998, 100, false)

function run_set_vs_L(image_list::Array{Int64}, KMAX::Int64=38, N::Int64=32, LMAX::Int64=18, K::Int64=26, temperature_decay::Float64=0.99998, ppi::Int64=10)
    histogram_list = Dict(
    "P$(img)" => histogram_name("parallel/data_generation/SH_", ppi, N, KMAX, K, float(KMAX), img, "") for img in images)
    for img in keys(histogram_list)
        for L=2:2:LMAX
            run_determination("paper_res_vs_L_$(ppi)p_KMAX$(KMAX)_N$(N)_K$(K)_$(temperature_decay)/$(img)_L$(L)", histograms=histogram_list[img], initial_stepsize=pi/4.0, K=K, L=L, KMAX=KMAX, rmax=float(KMAX), optimizer="rotate_all_at_once", initial_temperature_factor=0.1, temperature_decay=temperature_decay, N=N, successive_jobs=3, measure="Bayes", range=1000:1019, postprocess=true, gpu=true, Ncores=20, stepsizefactor=1.01, reference_pdb_path="$(ENV["THREEPHOTONS_PATH"])/data/structures/crambin.pdb")
        end
    end
end

# run_set_vs_L([5120000, 20480000, 81920000, 3276800000], 38, 32, 18, 26, 0.99998, 10)

function run_noise_set(sigmas::Array{Float64}, gammas::Array{Float64}, KMAX::Int64=38, N::Int64=32, L::Int64=18, K::Int64=26, temperature_decay::Float64=0.99998, ppi::Int64=10)
    for sigma in sigmas
        for gamma in gammas
            histo_name = histogram_name("parallel/data_generation/SH_", ppi, N, KMAX, K, float(KMAX), 3276800000, "", gamma, sigma)
            run_determination("paper_noise_$(ppi)p_KMAX$(KMAX)_N$(N)_K$(K)_L$(L)_$(temperature_decay)/G$(gamma)_S$(sigma)", histograms=histo_name, initial_stepsize=pi/4.0, K=K, L=L, KMAX=KMAX, rmax=float(KMAX), optimizer="rotate_all_at_once", initial_temperature_factor=0.1, temperature_decay=temperature_decay, N=N, successive_jobs=1, measure="Bayes", range=1000:1019, postprocess=true, gpu=true, Ncores=8, stepsizefactor=1.01, run_denoise=true, sigma=sigma, force_repostprocess=true, reference_pdb_path="$(ENV["THREEPHOTONS_PATH"])/data/structures/crambin.pdb")
        end
    end
end

  # run_noise_set([0.5, 0.75, 1.125, 2.5], collect(0.1:0.1:0.5), 38, 32, 18, 26, 0.99998, 10)


#For coliphage
# run_determination("exp_data/coliphage_determination", histograms="$(ENV["DETERMINATION_DATA"])/output_owl/exp_data/coliphage_symmetric_N32/histo.dat", lambda=0.0, initial_stepsize=pi/4.0, K=26, L=18, KMAX=38, rmax=float(38), optimizer="rotate_all_at_once", initial_temperature_factor=0.1, temperature_decay=0.99998, N=32, successive_jobs=3, measure="Bayes", range=1000:1019, postprocess=false, gpu=true, Ncores=20, stepsizefactor=1.01)

#For multi particle
# run_determination("multi_particle", histograms="$(ENV["DETERMINATION_DATA"])/output_owl/data_generation/multi_2_SH_10p_N32_K38_R38.0_P3276800000/histo.dat", lambda=0.0, initial_stepsize=pi/4.0, K=26, L=18, KMAX=38, rmax=float(38), optimizer="rotate_all_at_once", initial_temperature_factor=0.1, temperature_decay=0.99998, N=32, successive_jobs=3, measure="Bayes", range=1000:1019, postprocess=true, gpu=true, Ncores=20, stepsizefactor=1.01, reference_pdb_path="$(ENV["DETERMINATION_DATA"])/structures/crambin.pdb", force_repostprocess=true, run_denoise=true)

#Repeat original calculations for performance comparison
# run_determination("performance_CUDA", histograms="$(ENV["DETERMINATION_DATA"])/output_owl/data_generation/SH_10p_N32_K38_R38.0_P3276800000/histo.dat", lambda=0.0, initial_stepsize=pi/4.0, K=26, L=18, KMAX=38, rmax=float(38), optimizer="rotate_all_at_once", initial_temperature_factor=0.1, temperature_decay=0.99998, N=32, successive_jobs=1, measure="Bayes", range=1000:1000, postprocess=true, gpu=true, Ncores=20, stepsizefactor=1.01, reference_pdb_path="$(ENV["DETERMINATION_DATA"])/structures/crambin.pdb", force_repostprocess=true, run_denoise=true)

#For Ewald SH structure determination
# run_determination("Ewald_lambda_2.0_SH_10p_N32_K2_38_K3_26_R38.0_P3276800000", histograms="$(ENV["DETERMINATION_DATA"])/output_owl/data_generation/Ewald_lambda_2.0_SH_10p_N32_K2_38_K3_26_R38.0_P3276800000/histo.dat", lambda=2.0, initial_stepsize=pi/4.0, K=26, L=18, KMAX=38, rmax=float(38), optimizer="rotate_all_at_once", initial_temperature_factor=0.1, temperature_decay=0.99998, N=32, successive_jobs=3, measure="Bayes", range=1000:1019, postprocess=true, gpu=true, Ncores=20, stepsizefactor=1.01, reference_pdb_path="$(ENV["DETERMINATION_DATA"])/structures/crambin.pdb")

# run_determination("Ewald_lambda_4.0_alt_infinity", histograms="$(ENV["DETERMINATION_DATA"])/output/data_generation/correlations_N32_K238_K326_L18_inf.dat", lambda=4.0, initial_stepsize=pi/4.0, K=26, L=18, KMAX=38, rmax=float(38), optimizer="rotate_all_at_once", initial_temperature_factor=0.1, temperature_decay=0.99998, N=32, successive_jobs=3, measure="Bayes", range=1000:1019, postprocess=true, gpu=true, Ncores=20, stepsizefactor=1.01, reference_pdb_path="$(ENV["DETERMINATION_DATA"])/structures/crambin.pdb")
