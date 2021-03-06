function run_determination(dir::String; histograms::String="", initial_stepsize::Float64=Float64(pi), K3_range::UnitRange{Int64}=1:8, L::Integer=8, optimizer::String="rotate_hierarchical", initial_temperature_factor::Float64=1.0, temperature_decay::Float64=0.99, N::Integer=32, range=1000:1019, fresh::Bool=false, gpu::Bool=true, Ncores::Integer=8, successive_jobs::Integer=1, measure="Bayes", postprocess::Bool=true, stepsizefactor::Float64=1.02, K2_range::UnitRange{Int64}=1:35, qmax::Float64=1.0, force_repostprocess::Bool=false, run_denoise::Bool=false, architecture::String="", hours::Int64=48, sigma::Float64=0.0, reference_pdb_path::String="", lambda::Float64=0.0, include_negativity::Bool=false, gpu_mem::String="4000M")

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
            params["qmax"] = $qmax
            state["intensity"].rmax = $qmax
            postprocess_run(params, state, "$(reference_pdb_path)", true, 35, $sigma)
        end

        #Or start a completely new run
    else
        params,state = rotation_search(Dict( "reference_pdb_path"=>"$(reference_pdb_path)", "stepsizefactor"=>$stepsizefactor, "initial_stepsize" => $initial_stepsize, "L"=>$L, "K3_range"=>$(K3_range), "N"=>$N, "histograms"=>"$(histograms)", "optimizer"=>$optimizer, "initial_temperature_factor"=>$initial_temperature_factor, "measure"=>"$measure", "temperature_decay"=>$temperature_decay, "LMAX"=>25, "K2_range"=>$(K2_range), "qmax"=>$qmax, "lambda"=>$lambda, "include_negativity"=>$include_negativity))
        if $postprocess
            postprocess_run(params, state, "$(reference_pdb_path)", true, 35, $sigma)
        end
    end
    """

    for n in range
        launch_job("$dir/$n", Ncores, gpu, julia_script, successive_jobs; architecture=architecture, hours=hours, fresh=fresh, gpu_mem=gpu_mem)
    end
end

"""Fits all intensites with respect to the first for consecutive averaging"""
function run_postprocess_coliphage_results( dir::String, range::UnitRange{Int64}, L::Int64)
    for n in 1000:1019
        julia_script = """
        using ThreePhotons
        L=$L
        K3=maximum($range)
        intensity_reference = deserializeFromFile("$(environment_path("$dir"))/1000/state.dat")["intensity"]
        intensity = deserializeFromFile("$(environment_path("$dir"))/$n/state.dat")["intensity"]

        bestfit, bestsc, _, _, _ =  fitStructures_random(deleteTerms(intensity, K3, L), deleteTerms(intensity_reference, K3, L), $range, L, 0.995, nworkers()*8)
        serializeToFile("intensity.dat", bestfit)
        saveCube(bestfit, "fitted_intensity.mrc")
        """
        launch_job("$(dir)_fitted/$n", 8, false, julia_script, 1)
    end
end

function run_average_core_completion_phasing(name::String, center_range::UnitRange{Int64}, range::UnitRange{Int64}, coliphage_correlation_path::String)
    for beta_end in [0.8, 0.85, 0.9, 0.95, 0.99]
        julia_script = """
        using ThreePhotons

        _,_,_,_,c1 = loadHistograms(1,1,"$(environment_path("exp_data/$(coliphage_correlation_path)/histo.dat"))")

        extended_corrected_intensity = complete_core("$(environment_path("exp_data/$name"))", c1, $(center_range), $(range))
        averaged_density = phase_completed_intensity(extended_corrected_intensity, 8, $(beta_end))
        """
        launch_job("exp_data/$(name)_phased_center_$(minimum(center_range))_$(maximum(center_range))/b$(beta_end)", 8, false, julia_script, 1)
    end
end

function run_set(image_list::Array{Int64}, K2_range::UnitRange{Int64}=1:38, N::Int64=32, L::Int64=18, K3_range::UnitRange{Int64}=1:26, temperature_decay::Float64=0.99998, ppi::Int64=10, include_infinite::Bool=true)
    histograms_finite = Dict( "P$(img)" => histogram_name("parallel/data_generation/SH_", ppi, N, maximum(K2_range), maximum(K3_range), float(maximum(K2_range)), img, "") for img in image_list)
    histograms_infinite = include_infinite ? Dict("L20_inf" => "expdata/correlations_N$(N)_K$(maximum(K2_range))_L20_inf.dat") : Dict()
    histogram_list = merge(histograms_finite, histograms_infinite)

    for img in keys(histogram_list)
        run_determination("paper_res_vs_pictures_$(ppi)p_KMAX$(maximum(K2_range))_N$(N)_K$(maximum(K3_range))_L$(L)_$(temperature_decay)/$img", histograms=histogram_list[img], initial_stepsize=pi/4.0, K3_range=K3_range, L=L, K2_range=K2_range, qmax=qmax(maximum(K2_range), float(maximum(K2_range))), optimizer="rotate_all_at_once", initial_temperature_factor=0.1, temperature_decay=temperature_decay, N=N, successive_jobs=3, measure="Bayes", range=1000:1019, postprocess=true, gpu=true, Ncores=20, stepsizefactor=1.01, reference_pdb_path="$(ENV["THREEPHOTONS_PATH"])/data/structures/crambin.pdb", lambda=0.0, include_negativity=false)
    end
end

# run_set(calculate_images_ppi(10), 38, 32, 18, 26, 0.99998, 10)
# run_set(calculate_images_ppi(25), 38, 32, 18, 26, 0.99998, 25, false)
# run_set(calculate_images_ppi(50), 38, 32, 18, 26, 0.99998, 50, false)
# run_set(calculate_images_ppi(100), 38, 32, 18, 26, 0.99998, 100, false)

function run_set_vs_L(image_list::Array{Int64}, K2_range::UnitRange{Int64}=1:38, N::Int64=32, LMAX::Int64=18, K3_range::UnitRange{Int64}=1:26, temperature_decay::Float64=0.99998, ppi::Int64=10)
    histogram_list = Dict(
    "P$(img)" => histogram_name("parallel/data_generation/SH_", ppi, N, maximum(K2_range), maximum(K3_range), float(maximum(K2_range)), img, "") for img in images)
    for img in keys(histogram_list)
        for L=2:2:LMAX
            run_determination("paper_res_vs_L_$(ppi)p_KMAX$(maximum(K2_range))_N$(N)_K$(maximum(K3_range))_$(temperature_decay)/$(img)_L$(L)", histograms=histogram_list[img], initial_stepsize=pi/4.0, K3_range=K3_range, L=L, K2_range=K2_range, qmax=qmax(maximum(K2_range), float(maximum(K2_range))), optimizer="rotate_all_at_once", initial_temperature_factor=0.1, temperature_decay=temperature_decay, N=N, successive_jobs=3, measure="Bayes", range=1000:1019, postprocess=true, gpu=true, Ncores=20, stepsizefactor=1.01, reference_pdb_path="$(ENV["THREEPHOTONS_PATH"])/data/structures/crambin.pdb", lambda=0.0)
        end
    end
end

# run_set_vs_L([5120000, 20480000, 81920000, 3276800000], 38, 32, 18, 26, 0.99998, 10)

function run_noise_set(sigmas::Array{Float64}, gammas::Array{Float64}, K2_range::UnitRange{Int64}=1:38, N::Int64=32, L::Int64=18, K3_range::UnitRange{Int64}=1:26, temperature_decay::Float64=0.99998, ppi::Int64=10)
    for sigma in sigmas
        for gamma in gammas
            histo_name = histogram_name("parallel/data_generation/SH_", ppi, N, maximum(K2_range), maximum(K3_range), float(maximum(K2_range)), 3276800000, "", gamma, sigma)
            run_determination("paper_noise_$(ppi)p_KMAX$(maximum(K2_range))_N$(N)_K$(maximum(K3_range))_L$(L)_$(temperature_decay)/G$(gamma)_S$(sigma)", histograms=histo_name, initial_stepsize=pi/4.0, K3_range=K3_range, L=L, K2_range=K2_range, qmax=qmax(maximum(K2_range), float(maximum(K2_range))), optimizer="rotate_all_at_once", initial_temperature_factor=0.1, temperature_decay=temperature_decay, N=N, successive_jobs=1, measure="Bayes", range=1000:1019, postprocess=true, gpu=true, Ncores=8, stepsizefactor=1.01, run_denoise=true, sigma=sigma, force_repostprocess=true, reference_pdb_path="$(ENV["THREEPHOTONS_PATH"])/data/structures/crambin.pdb")
        end
    end
end

  # run_noise_set([0.5, 0.75, 1.125, 2.5], collect(0.1:0.1:0.5), 38, 32, 18, 26, 0.99998, 10)

#For multi particle
# run_determination("multi_particle", histograms="$(ENV["DETERMINATION_DATA"])/output_owl/data_generation/multi_2_SH_10p_N32_K38_R38.0_P3276800000/histo.dat", lambda=0.0, initial_stepsize=pi/4.0, K3_range=1:26, L=18, K2_range=1:38, qmax=qmax(38, float(38)), optimizer="rotate_all_at_once", initial_temperature_factor=0.1, temperature_decay=0.99998, N=32, successive_jobs=3, measure="Bayes", range=1000:1019, postprocess=true, gpu=true, Ncores=20, stepsizefactor=1.01, reference_pdb_path="$(ENV["DETERMINATION_DATA"])/structures/crambin.pdb", force_repostprocess=true, run_denoise=true)

#Repeat original calculations for performance comparison
# run_determination("performance_CUDA", histograms="$(ENV["DETERMINATION_DATA"])/output_owl/data_generation/SH_10p_N32_K38_R38.0_P3276800000/histo.dat", lambda=0.0, initial_stepsize=pi/4.0, K3_range=1:26, L=18, K2_range=1:38, qmax=qmax(38, float(38)), optimizer="rotate_all_at_once", initial_temperature_factor=0.1, temperature_decay=0.99998, N=32, successive_jobs=1, measure="Bayes", range=1000:1000, postprocess=true, gpu=true, Ncores=20, stepsizefactor=1.01, reference_pdb_path="$(ENV["DETERMINATION_DATA"])/structures/crambin.pdb", force_repostprocess=true, run_denoise=true)

#For Ewald SH structure determination
# run_determination("Ewald_lambda_2.0_SH_10p_N32_K2_38_K3_26_R38.0_P3276800000", histograms="$(ENV["DETERMINATION_DATA"])/output_owl/data_generation/Ewald_lambda_2.0_SH_10p_N32_K2_38_K3_26_R38.0_P3276800000/histo.dat", lambda=2.0, initial_stepsize=pi/4.0, K3_range=1:26, L=18, K2_range=1:38, qmax=qmax(38, float(38)), optimizer="rotate_all_at_once", initial_temperature_factor=0.1, temperature_decay=0.99998, N=32, successive_jobs=3, measure="Bayes", range=1000:1019, postprocess=true, gpu=true, Ncores=20, stepsizefactor=1.01, reference_pdb_path="$(ENV["DETERMINATION_DATA"])/structures/crambin.pdb")

# run_determination("Ewald_lambda_4.0_alt_infinity", histograms="$(ENV["DETERMINATION_DATA"])/output/data_generation/correlations_N32_K238_K326_L18_inf.dat", lambda=4.0, initial_stepsize=pi/4.0, K3_range=1:26, L=18, K2_range=1:38, qmax=qmax(38, float(38)), optimizer="rotate_all_at_once", initial_temperature_factor=0.1, temperature_decay=0.99998, N=32, successive_jobs=3, measure="Bayes", range=1000:1019, postprocess=true, gpu=true, Ncores=20, stepsizefactor=1.01, reference_pdb_path="$(ENV["DETERMINATION_DATA"])/structures/crambin.pdb")

#Coliphage with correct qmax and lmabda, lambda=7.75 A, max_resolution=116 A, qmax= pi / 116
# run_determination("exp_data/coliphage_determination_newhisto", histograms="$(ENV["DETERMINATION_DATA"])/output_owl/exp_data/coliphage/histo.dat", lambda=0.0, initial_stepsize=pi/4.0, K3_range=5:26, L=12, K2_range=5:38, qmax=pi/116.0, optimizer="rotate_all_at_once", initial_temperature_factor=0.1, temperature_decay=0.99998, N=32, successive_jobs=3, measure="Bayes", range=1000:1019, postprocess=false, gpu=true, Ncores=20, stepsizefactor=1.01, include_negativity=false)

#Coliphage with symmetric and L=12
# run_determination("exp_data/coliphage_determination_newhisto_symmetric", histograms="$(ENV["DETERMINATION_DATA"])/output_owl/exp_data/coliphage_symmetric/histo.dat", lambda=0.0, initial_stepsize=pi/3.0, K3_range=6:26, L=12, K2_range=6:38, qmax=pi/116.0, optimizer="rotate_all_at_once", initial_temperature_factor=0.1, temperature_decay=0.99998, N=32, successive_jobs=3, measure="Bayes", range=1000:1019, postprocess=false, gpu=true, Ncores=20, stepsizefactor=1.01, include_negativity=false)

#Coliphage with symmetric and L=12, K3=32
# run_determination("exp_data/coliphage_determination_newhisto_symmetric_K32", histograms="$(ENV["DETERMINATION_DATA"])/output_owl/exp_data/coliphage_symmetric/histo.dat", lambda=0.0, initial_stepsize=pi/3.0, K3_range=6:32, L=12, K2_range=6:38, qmax=pi/116.0, optimizer="rotate_all_at_once", initial_temperature_factor=0.1, temperature_decay=0.99998, N=32, successive_jobs=3, measure="Bayes", range=1000:1019, postprocess=false, gpu=true, Ncores=20, stepsizefactor=1.01, include_negativity=false)

#Coliphage with symmetric and L=12 and low range
# run_determination("exp_data/coliphage_determination_paper_lowrange", histograms="$(ENV["DETERMINATION_DATA"])/output_owl/exp_data/coliphage_K2_38_K3_30_N32/histo.dat", lambda=0.0, initial_stepsize=pi/4.0, K3_range=8:26, L=12, K2_range=8:38, qmax=pi/90.0, optimizer="rotate_all_at_once", initial_temperature_factor=0.1, temperature_decay=0.99998, N=32, successive_jobs=3, measure="Bayes", range=1000:1019, postprocess=false, gpu=true, Ncores=20, stepsizefactor=1.01, include_negativity=false)

#Coliphage with symmetric and L=12 and low range
# run_determination("exp_data/coliphage_determination_paper_lowrange_2", histograms="$(ENV["DETERMINATION_DATA"])/output_owl/exp_data/coliphage_K2_38_K3_30_N32/histo.dat", lambda=0.0, initial_stepsize=pi/4.0, K3_range=7:26, L=12, K2_range=3:38, qmax=pi/90.0, optimizer="rotate_all_at_once", initial_temperature_factor=0.25, temperature_decay=0.99998, N=32, successive_jobs=3, measure="Bayes", range=1000:1019, postprocess=false, gpu=true, Ncores=20, stepsizefactor=1.01, include_negativity=false)

#Coliphage with symmetric and L=12 and high range
# run_determination("exp_data/coliphage_determination_paper_K2_5_38_K3_5_30", histograms="$(ENV["DETERMINATION_DATA"])/output_owl/exp_data/coliphage_K2_38_K3_30_N32/histo.dat", lambda=0.0, initial_stepsize=pi/4.0, K3_range=5:30, L=12, K2_range=5:38, qmax=pi/90.0, optimizer="rotate_all_at_once", initial_temperature_factor=0.2, temperature_decay=0.99998, N=32, successive_jobs=3, measure="Bayes", range=1000:1019, postprocess=false, gpu=true, Ncores=20, stepsizefactor=1.01, include_negativity=false)

#Coliphage with symmetric and L=12 and high range
# run_determination("exp_data/coliphage_determination_paper_K2_7_35_K3_7_30", histograms="$(ENV["DETERMINATION_DATA"])/output_owl/exp_data/coliphage_K2_38_K3_30_N32/histo.dat", lambda=0.0, initial_stepsize=pi/4.0, K3_range=7:30, L=12, K2_range=7:35, qmax=pi/90.0, optimizer="rotate_all_at_once", initial_temperature_factor=0.2, temperature_decay=0.99998, N=32, successive_jobs=3, measure="Bayes", range=1000:1019, postprocess=false, gpu=true, Ncores=20, stepsizefactor=1.01, include_negativity=false)

#Coliphage with symmetric and L=12 and high range
# run_determination("exp_data/coliphage_determination_paper_K2_10_35_K3_10_30", histograms="$(ENV["DETERMINATION_DATA"])/output_owl/exp_data/coliphage_K2_38_K3_30_N32/histo.dat", lambda=0.0, initial_stepsize=pi/4.0, K3_range=10:30, L=12, K2_range=10:35, qmax=pi/90.0, optimizer="rotate_all_at_once", initial_temperature_factor=0.1, temperature_decay=0.99998, N=32, successive_jobs=3, measure="Bayes", range=1000:1019, postprocess=false, gpu=true, Ncores=20, stepsizefactor=1.01, include_negativity=false)

# run_determination("exp_data/coliphage_determination_paper_K2_10_36_K3_10_26", histograms="$(ENV["DETERMINATION_DATA"])/output_owl/exp_data/coliphage_K2_38_K3_30_N32/histo.dat", lambda=0.0, initial_stepsize=pi/4.0, K3_range=10:26, L=12, K2_range=10:36, qmax=pi/90.0, optimizer="rotate_all_at_once", initial_temperature_factor=0.1, temperature_decay=0.99998, N=32, successive_jobs=3, measure="Bayes", range=1000:1019, postprocess=false, gpu=true, Ncores=20, stepsizefactor=1.01, include_negativity=false)

#Coliphage with L=12, K=12
# run_determination("exp_data/coliphage_determination_paper_K2_12_38_K3_12_26_L12", histograms="$(ENV["DETERMINATION_DATA"])/output_owl/exp_data/coliphage_K2_38_K3_30_N32/histo.dat", lambda=0.0, initial_stepsize=pi/4.0, K3_range=12:26, L=12, K2_range=12:38, qmax=pi/90.0, optimizer="rotate_all_at_once", initial_temperature_factor=0.1, temperature_decay=0.99998, N=32, successive_jobs=3, measure="Bayes", range=1000:1019, postprocess=false, gpu=true, Ncores=20, stepsizefactor=1.01, include_negativity=false)

#Coliphage smaller dq for L=14 and L=12
# run_determination("exp_data/coliphage_determination_paper2_K2_12_40_K3_12_30_L14", histograms="$(ENV["DETERMINATION_DATA"])/output_owl/exp_data/coliphage_K2_50_K3_30_N32/histo.dat", lambda=0.0, initial_stepsize=pi/4.0, K3_range=12:30, L=14, K2_range=12:40, qmax=pi/90.0, optimizer="rotate_all_at_once", initial_temperature_factor=0.1, temperature_decay=0.99998, N=32, successive_jobs=3, measure="Bayes", range=1000:1019, postprocess=false, gpu=true, Ncores=8, stepsizefactor=1.01, include_negativity=false, gpu_mem="4000M")

#Coliphage best bet
# run_determination("exp_data/coliphage_determination_rerun_K2_5_38_K3_5_28_L12", histograms="$(ENV["DETERMINATION_DATA"])/output_owl/exp_data/coliphage_K2_40_K3_30_N32_processed/histo.dat", lambda=0.0, initial_stepsize=pi/4.0, K3_range=5:28, L=12, K2_range=5:38, qmax=pi/90.0, optimizer="rotate_all_at_once", initial_temperature_factor=0.1, temperature_decay=0.99998, N=32, successive_jobs=3, measure="Bayes", range=1000:1019, postprocess=false, gpu=true, Ncores=8, stepsizefactor=1.01, include_negativity=false, gpu_mem="4000M")

#Rerun Ewald stuff
# run_determination("Ewald_P327680000_lambda_2.5", histograms="$(ENV["DETERMINATION_DATA"])/output_owl/data_generation/Ewald_large_SH_10p_N32_K2_38_K3_26_R38.0_P327680000_lambda_2.5/histo.dat", lambda=2.5, initial_stepsize=pi/4.0, K3_range=1:26, L=18, K2_range=1:38, qmax=qmax(38, float(38)), optimizer="rotate_all_at_once", initial_temperature_factor=0.1, temperature_decay=0.99998, N=32, successive_jobs=3, measure="Bayes", range=1000:1019, postprocess=true, gpu=true, Ncores=8, stepsizefactor=1.01, reference_pdb_path="$(ENV["DETERMINATION_DATA"])/structures/crambin.pdb", include_negativity=false, gpu_mem="4000M")
