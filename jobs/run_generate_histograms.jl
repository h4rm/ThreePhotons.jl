include("runs.jl")

if ENV_name == "gwdg" || ENV_name == "owl"

  function generate_histogram_set_ppi(ppi::Int64=10; K::Int64=38, N::Int64=32, incident_10p::Int64=550)
    images = calculate_images_ppi(ppi)
    setsize = images[4]
    for img in images
      if img <= images[4]
        generate_histograms(; max_pictures = img, max_triplets = Integer(0), Ncores=8, N=N, photons_per_image=ppi, batchsize = Integer(img/8), successive_jobs=1, prefix="SH_", suffix="", use_cube=false, qcut_ratio=1.0, K=K, rmax=float(K), histogram_method="histogramCorrelationsInPicture_alltoall", structure_pdb_path="$(ENV["THREEPHOTONS_PATH"])/structures/crambin.pdb")

      else
        numbersets = ceil(Int64, img / setsize)
        for i = 1:numbersets
          generate_histograms(; max_pictures = setsize, max_triplets = Integer(0), Ncores=8, N=N, photons_per_image=ppi, batchsize = Integer(setsize/8), successive_jobs=1, prefix="parts/SH_", suffix="_$(i)", use_cube=false, qcut_ratio=1.0, K=K, rmax=float(K), histogram_method="histogramCorrelationsInPicture_alltoall", structure_pdb_path="$(ENV["THREEPHOTONS_PATH"])/structures/crambin.pdb")
        end
      end
    end
  end

  # generate_histogram_set_ppi(25)
  # generate_histogram_set_ppi(50)

  noise_photons = Dict(0.5=>14, 0.75=>14, 1.125=>14, 2.5=>25)
  let K = 38, N = 32
    # for sigma in [0.5, 0.75, 1.125]
    let sigma = 2.5
      for gamma in collect(0.1:0.1:0.5)
        setsize = Integer(2*2.048e7)
        numbersets = Integer(3.2768e9 / setsize)
        for i = 1:numbersets
          generate_histograms(; max_pictures = setsize, max_triplets = Integer(0), Ncores=8, N=N, photons_per_image=10, batchsize = Integer(2.56e5), successive_jobs=3, prefix="parts/SH_", suffix="_$(i)", use_cube=false, qcut_ratio=1.0, K=K, rmax=float(K), histogram_method="histogramCorrelationsInPicture_alltoall", gamma=gamma, sigma=sigma, noise_photons=round(Int64, noise_photons[sigma] * gamma), structure_pdb_path="$(ENV["THREEPHOTONS_PATH"])/structures/crambin.pdb")
        end
      end
    end
  end

end
