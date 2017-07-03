include("runs.jl")

if ENV_name == "owl" || ENV_name == "gwdg"

  function run_set(image_list::Array{Int64}, KMAX::Int64=38, N::Int64=32, L::Int64=18, K::Int64=26, temperature_decay::Float64=0.99998, ppi::Int64=10, include_infinite::Bool=true)
    histograms_finite = Dict( "P$(img)" => histogram_name("parallel/data_generation/SH_", ppi, N, KMAX, float(KMAX), img, "") for img in image_list)
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
    "P$(img)" => histogram_name("parallel/data_generation/SH_", ppi, N, KMAX, float(KMAX), img, "") for img in images)
    for img in keys(histogram_list)
      for L=2:2:LMAX
        run_determination("paper_res_vs_L_$(ppi)p_KMAX$(KMAX)_N$(N)_K$(K)_$(temperature_decay)/$(img)_L$(L)", histograms=histogram_list[img], initial_stepsize=pi/4.0, K=K, L=L, KMAX=KMAX, rmax=float(KMAX), optimizer="rotate_all_at_once", initial_temperature_factor=0.1, temperature_decay=temperature_decay, N=N, successive_jobs=3, measure="Bayes", range=1000:1019, postprocess=true, gpu=true, Ncores=20, stepsizefactor=1.01, reference_pdb_path="$(ENV["DETERMINATION_PATH"])/data/structures/crambin.pdb")
      end
    end
  end

  # run_set_vs_L([5120000, 20480000, 81920000, 3276800000], 38, 32, 18, 26, 0.99998, 10)

  function run_noise_set(sigmas::Array{Float64}, gammas::Array{Float64}, KMAX::Int64=38, N::Int64=32, L::Int64=18, K::Int64=26, temperature_decay::Float64=0.99998, ppi::Int64=10)
    for sigma in sigmas
      for gamma in gammas
        histo_name = histogram_name("parallel/data_generation/SH_", ppi, N, KMAX, float(KMAX), 3276800000, "", gamma, sigma)
        run_determination("paper_noise_$(ppi)p_KMAX$(KMAX)_N$(N)_K$(K)_L$(L)_$(temperature_decay)/G$(gamma)_S$(sigma)", histograms=histo_name, initial_stepsize=pi/4.0, K=K, L=L, KMAX=KMAX, rmax=float(KMAX), optimizer="rotate_all_at_once", initial_temperature_factor=0.1, temperature_decay=temperature_decay, N=N, successive_jobs=1, measure="Bayes", range=1000:1019, postprocess=true, gpu=true, Ncores=8, stepsizefactor=1.01, run_denoise=true, sigma=sigma, force_repostprocess=true, reference_pdb_path="$(ENV["DETERMINATION_PATH"])/data/structures/crambin.pdb")
      end
    end
  end

  run_noise_set([0.5, 0.75, 1.125, 2.5], collect(0.1:0.1:0.5), 38, 32, 18, 26, 0.99998, 10)

end
