include("runs.jl")

function generate_histogram_image(img::Int64, ppi::Int64, K::Int64, N::Int64; setsize::Int64=Integer(2*2.048e7), name::String="", lambda::Float64=1.0)
    if img <= setsize
        generate_histograms(; max_pictures = img, max_triplets = Integer(0), Ncores=8, N=N, photons_per_image=ppi, batchsize = Integer(img/8), successive_jobs=1, prefix="SH_", suffix="", use_cube=false, qcut_ratio=1.0, K=K, rmax=float(K), histogram_method="histogramCorrelationsInPicture_alltoall", structure_pdb_path="$(ENV["DETERMINATION_DATA"])/structures/crambin.pdb", lambda=lambda)

    else
        numbersets = ceil(Int64, img / setsize)
        for i = 1:numbersets
            generate_histograms(; max_pictures = setsize, max_triplets = Integer(0), Ncores=8, N=N, photons_per_image=ppi, batchsize = Integer(setsize/8), successive_jobs=1, prefix="parts/$(name)SH_", suffix="_$(i)", use_cube=false, qcut_ratio=1.0, K=K, rmax=float(K), histogram_method="histogramCorrelationsInPicture_alltoall", structure_pdb_path="$(ENV["DETERMINATION_DATA"])/structures/crambin.pdb", lambda=lambda)
        end
    end
end

function generate_histogram_set_ppi(ppi::Int64=10; K::Int64=38, N::Int64=32)
    for img in calculate_images_ppi(ppi)
        generate_histogram_image(img, ppi, K, N)
    end
end

# generate_histogram_set_ppi(25)
# generate_histogram_set_ppi(50)

function generate_noisy_histograms()
    noise_photons = Dict(0.5=>14, 0.75=>14, 1.125=>14, 2.5=>25)
    let K = 38, N = 32
        # for sigma in [0.5, 0.75, 1.125]
        let sigma = 2.5
            for gamma in collect(0.1:0.1:0.5)
                setsize = Integer(2*2.048e7)
                numbersets = Integer(3.2768e9 / setsize)
                for i = 1:numbersets
                    generate_histograms(; max_pictures = setsize, max_triplets = Integer(0), Ncores=8, N=N, photons_per_image=10, batchsize = Integer(2.56e5), successive_jobs=3, prefix="parts/SH_", suffix="_$(i)", use_cube=false, qcut_ratio=1.0, K=K, rmax=float(K), histogram_method="histogramCorrelationsInPicture_alltoall", gamma=gamma, sigma=sigma, noise_photons=round(Int64, noise_photons[sigma] * gamma), structure_pdb_path="$(ENV["DETERMINATION_DATA"])/structures/crambin.pdb")
                end
            end
        end
    end
end

# generate_noisy_histograms()

function generate_single_multiparticle_histogram(number_images::Int64, setsize::Int64, number_particles::Int64=2, ppi::Int64=10; K::Int64=38, N::Int64=32, lambda::Float64=0.0)
    numbersets = ceil(Int64, number_images / setsize)
    for i = 1:numbersets
        generate_histograms(; max_pictures = setsize, max_triplets = Integer(0), Ncores=8, N=N, photons_per_image=ppi, batchsize = Integer(setsize/8), successive_jobs=1, prefix="parts/multi_$(number_particles)_SH_", suffix="_$(i)", use_cube=false, qcut_ratio=1.0, K=K, rmax=float(K), histogram_method="histogramCorrelationsInPicture_alltoall", structure_pdb_path="$(ENV["DETERMINATION_DATA"])/structures/crambin.pdb", number_particles=number_particles, lambda=lambda)
    end
end

# generate_single_multiparticle_histogram(Integer(3.2768e8), Integer(2*2.048e6), 2)
# combine_set([Integer(3.2768e8)], Integer(2*2.048e6), 10, 38, 32, environment_path("data_generation/parts/multi_2_SH_"))

# generate_histogram_image(Integer(3.2768e8), 10, 38, 32; setsize=Integer(2*2.048e6), name="Ewald_", lambda=1.0)
# combine_set([Integer(3.2768e8)], Integer(2*2.048e6), 10, 38, 32, environment_path("data_generation/parts/Ewald_SH_"))
