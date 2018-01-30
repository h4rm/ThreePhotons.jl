function histogram_name(prefix::String, ppi::Int64, N::Int64, K2::Int64, K3::Int64, rmax::Float64, max_pictures::Int64, suffix::String, gamma::Float64=0.0, sigma::Float64=1.0, lambda::Float64=0.0)
    return "$(prefix)$(ppi)p_N$(N)_K2_$(K2)_K3_$(K3)_R$(rmax)_P$(max_pictures)$(gamma > 0.0 ? "_G$(gamma)_S$(sigma)" : "")_lambda_$(lambda)$(suffix)"
end

"""Starts a cluster job for synthetic correlation generation"""
function generate_histograms(; max_pictures::Integer=Integer(0), N::Integer=32, number_incident_photons::Integer=500, incident_photon_variance::Integer = 0, gamma::Float64=0.0, sigma::Float64=1.0, noise_photons::Int64=0, Ncores::Integer=8, batchsize::Integer = Integer(1e4), successive_jobs::Integer=1, prefix::String="correlations_", suffix::String="", use_cube::Bool=true, qcut_ratio::Float64=1.0, K2::Int64=38, K3::Int64=26, rmax::Float64=35.0, histogram_method="histogramCorrelationsInPicture_alltoall", structure_pdb_path::String="", number_particles::Int64=1, lambda::Float64=0.0, beamstop_width::Float64=0.0, beamstop_only::Bool=false)
    name = histogram_name(prefix, calculate_ppi(number_incident_photons), N, K2, K3, rmax, max_pictures, suffix, gamma, sigma, lambda)

    julia_script = """
    using ThreePhotons

    volume = 0.0

    if $(use_cube)
        # _,_,volume = createCubicStructure("$(structure_pdb_path)", 4*$K2+1, 2.0*$rmax)
        volume = loadCube("$(ENV["DETERMINATION_DATA"])/histograms/synthetic_crambin/intensityCube_high.mrc")
    else
        _,_,intensity = createSphericalHarmonicsStructure("$(structure_pdb_path)", 35, $K2, $rmax)
        volume = getSurfaceVolume(intensity)
        volume.radial_interp = false
    end

    generateHistogram(volume; qcut=$(qcut_ratio)*volume.rmax, K2=$K2, K3=$K3, N=$N, max_pictures=$max_pictures, number_incident_photons=$number_incident_photons, incident_photon_variance=$incident_photon_variance, numprocesses=$(Ncores), file="histo.dat", noise=GaussianNoise($gamma, $sigma, $(noise_photons)), batchsize = $batchsize, histogramMethod=$histogram_method, number_particles=$(number_particles), lambda=$(lambda), beamstop_width=$(beamstop_width), beamstop_only=$(beamstop_only))
    """

    launch_job("data_generation/$name", Ncores, false, julia_script, successive_jobs, architecture="")
end

function generate_histogram_image(img::Int64, ppi::Int64, K2::Int64, K3::Int64, N::Int64, lambda::Float64=0.0; setsize::Int64=Integer(2*2.048e7), name::String="", beamstop_width::Float64=0.0)
    numbersets = ceil(Int64, img / setsize)
    for i = 1:numbersets
        generate_histograms(; max_pictures = setsize, Ncores=8, N=N, number_incident_photons=calculate_incident_photons(ppi), batchsize = Integer(setsize/8), successive_jobs=1, prefix="parts/$(name)SH_", suffix="_$(i)", use_cube=false, qcut_ratio=1.0, K2=K2, K3=K3, rmax=float(K2), histogram_method="histogramCorrelationsInPicture_alltoall", structure_pdb_path="$(ENV["DETERMINATION_DATA"])/structures/crambin.pdb", lambda=lambda, beamstop_width=beamstop_width)
    end
end

function generate_histogram_set_ppi(ppi::Int64=10; K2::Int64=38, K3::Int64=26, N::Int64=32, lambda::Float64=0.0)
    for img in calculate_images_ppi(ppi)
        generate_histogram_image(img, ppi, K2, K3, N, lambda)
    end
end

# generate_histogram_set_ppi(10, 38, 26, 32)
# generate_histogram_set_ppi(25)
# generate_histogram_set_ppi(50)

function generate_noisy_histograms()
    noise_photons = Dict(0.5=>14, 0.75=>14, 1.125=>14, 2.5=>25)
    let K2 = 38, K3 = 26, N = 32
        # for sigma in [0.5, 0.75, 1.125]
        let sigma = 2.5
            for gamma in collect(0.1:0.1:0.5)
                setsize = Integer(2*2.048e7)
                numbersets = Integer(3.2768e9 / setsize)
                for i = 1:numbersets
                    generate_histograms(; max_pictures = setsize, Ncores=8, N=N, number_incident_photons=calculate_incident_photons(10), batchsize = Integer(2.56e5), successive_jobs=3, prefix="parts/SH_", suffix="_$(i)", use_cube=false, qcut_ratio=1.0, K2=K2, K3=K3, rmax=float(K2), histogram_method="histogramCorrelationsInPicture_alltoall", gamma=gamma, sigma=sigma, noise_photons=round(Int64, noise_photons[sigma] * gamma), structure_pdb_path="$(ENV["DETERMINATION_DATA"])/structures/crambin.pdb")
                end
            end
        end
    end
end

# generate_noisy_histograms()

function generate_single_multiparticle_histogram(number_images::Int64, setsize::Int64, number_particles::Int64=2, ppi::Int64=10; K2::Int64=38, K3::Int64=26, N::Int64=32, lambda::Float64=0.0)
    numbersets = ceil(Int64, number_images / setsize)
    for i = 1:numbersets
        generate_histograms(; max_pictures = setsize, Ncores=8, N=N, number_incident_photons=calculate_incident_photons(ppi), batchsize = Integer(setsize/8), successive_jobs=1, prefix="parts/multi_$(number_particles)_SH_", suffix="_$(i)", use_cube=false, qcut_ratio=1.0, K2=K2, K3=K3, rmax=float(K2), histogram_method="histogramCorrelationsInPicture_alltoall", structure_pdb_path="$(ENV["DETERMINATION_DATA"])/structures/crambin.pdb", number_particles=number_particles, lambda=lambda)
    end
end

using HDF5
"""Distributes the calculation of correlations among many jobs"""
function run_calculate_correlation_from_images(particle_name::String, images_path::String, images_per_job::Int64, K2::Int64, K3::Int64, N::Int64; Ncores::Int64=8, symmetrize::Bool=false, sample_photons::Int64=Integer(1e6))

    file = h5open(images_path, "r")
    photonConverter = read(file["photonConverter"])
    _,_,images_in_file = Base.size(photonConverter["pnccdBack"]["photonCount"])
    number_runs = floor(Int64, images_in_file / images_per_job)

    for n in 1:number_runs
        julia_script = """
        using HDF5
        using ThreePhotons
        using Images

        K2 = $K2
        K3 = $K3
        N = $N

        photonConverter = h5open("$(images_path)", "r") do file
            read(file["photonConverter"])
        end

        #Calculate maximum intensity integral for scaling
        overall_maximum = maximum(Float64[sum(abs,photonConverter["pnccdBack"]["photonCount"][:,:,i]) for i = 1:500])

        resized_image_list = [ Images.imresize(convert(Array{Float64},photonConverter["pnccdBack"]["photonCount"][:,:,i]), (2*K2+1, 2*K2+1)) for i=$((n-1)*images_per_job+1):$(n*images_per_job)]
        calculate_correlations_in_image_using_single_photons(resized_image_list, K2, K3, N, "histo.dat", overall_maximum, Integer($(sample_photons)), $(symmetrize))
        """
        launch_job("exp_data/parts/$(particle_name)_$(n)", Ncores, false, julia_script, 1)#, memory="$(Ncores*1.5)G")
    end
end

"""Distributes the calculation of correlations among many jobs"""
function run_calculate_beamstop_correlation(jobname::String, images_path::String, K2::Int64, K3::Int64, N::Int64; Ncores::Int64=8)

    julia_script = """
    using HDF5
    using ThreePhotons
    using Images

    K2 = $K2
    K3 = $K3
    N = $N

    file = h5open("$(images_path)", "r")
    photonConverter = read(file["photonConverter"])

    sum = 0
    for i = 1:500
        sum += convert(Array{Float64},photonConverter["pnccdBack"]["photonCount"][:,:,i])
    end
    beamstop_map = (convert(Array{Float64},(abs(sum) .< eps())))
    beamstop_map_resized = Images.imresize(convert(Images.Image,beamstop_map), (2*K2+1, 2*K2+1)).data

    calculate_correlations_in_image([1.0 - beamstop_map_resized], K2, K3, N)
    """
    launch_job("exp_data/$(jobname)_K2_$(K2)_K3_$(K3)_N$(N)", Ncores, false, julia_script, 1)
end

# """Distributes the calculation of correlations among many jobs"""
# function run_calculate_generic_beamstop_correlation(jobname::String, K2::Int64, K3::Int64, N::Int64; Ncores::Int64=8)
#
#     julia_script = """
#     using ThreePhotons
#     using Images
#
#     K2 = $K2
#     K3 = $K3
#     N = $N
#
#     function beamstop(K::Vector{Float64}, beamstop_width::Float64)
#         (abs(K[1]) > beamstop_width && abs(K[2]) > beamstop_width) ? 1.0 : 0.0
#     end
#
#     qm = 3.141592653589793
#     delq = 0.08267349088394192
#     range = -qm:delq:qm
#     beamstop_width = qmax(K2,float(K2))/20.0
#     crambin_beamstop = [beamstop([k1,k2], beamstop_width) for k1 in range, k2 in range]
#
#     calculate_correlations_in_image([crambin_beamstop], K2, K3, N)
#     """
#     launch_job("exp_data/$(jobname)", Ncores, false, julia_script, 1)
# end

exp_filelist = String[
"amo86615_186_PR772_single.h5",
"amo86615_188_PR772_single.h5",
"amo86615_190_PR772_single.h5",
"amo86615_191_PR772_single.h5",
"amo86615_192_PR772_single.h5",
"amo86615_193_PR772_single.h5",
"amo86615_194_PR772_single.h5",
"amo86615_196_PR772_single.h5",
"amo86615_197_PR772_single.h5"
]

# for file in exp_filelist
#     run_calculate_correlation_from_images("coliphage_single_photons_K2_40_K3_30_N32/$file", environment_path("exp_data/Coliphage_PR772/$file"), 100, 40, 30, 32, symmetrize=false)
# end

#Calculate beamstop of Coliphage_PR772
# run_calculate_beamstop_correlation("coliphage_beamstop", environment_path("exp_data/Coliphage_PR772/amo86615_186_PR772_single.h5"), 40, 30, 32)

#Calculate beamstop of crambin
# run_calculate_generic_beamstop_correlation("crambin_beamstop", 38, 26, 32)

function combine_histograms(dir::String, num::Int64)
    println("Processing $dir")

    println("\tLoading file #1")
    params,c2_full,c3_full,c1_full = deserializeFromFile("$(dir)_1/histo.dat")
    part_images = params["num_pictures"]

    for i = 2:num
        println("\tLoading file #$i")
        histofile = "$(dir)_$(i)/histo.dat"
        if isfile(histofile)
            params_part, c2_part, c3_part, c1_part =  deserializeFromFile(histofile)
            part_triplets = sum(c3_part)
            if part_triplets < float(1.0e20)
                params["num_pictures"] += params_part["num_pictures"]
                c1_full += c1_part
                c2_full += c2_part
                c3_full += c3_part
                println("\tAdded $(part_triplets) triplets and $(params_part["num_pictures"]) images.")
            else
                println("\t------ERROR: $i file.")
            end
        end
    end

    finalpath = replace(dir, "parts/", "")
    finalpath = replace(finalpath, "P$(part_images)", "P$(params["num_pictures"])")

    println("Writing to $(finalpath)/histo.dat")
    try mkdir(finalpath) end
    serializeToFile("$(finalpath)/histo.dat", (params, c2_full, c3_full, c1_full))
end

function combine_histograms(filelist::Array{String}, output_dir::String)

    println("Processing $(filelist[1]).")
    params,c2_full,c3_full,c1_full = deserializeFromFile(filelist[1])
    part_images = params["num_pictures"]

    for i = 2:length(filelist)
        histofile = filelist[i]
        println("Processing $histofile.")
        if isfile("$histofile.jld")
            try
                params_part, c2_part, c3_part, c1_part =  deserializeFromFile(histofile)
                part_triplets = sum(c3_part)
                if part_triplets < float(1.0e20)
                    params["num_pictures"] += params_part["num_pictures"]
                    c1_full += c1_part
                    c2_full += c2_part
                    c3_full += c3_part
                    println("\tAdded $(part_triplets) triplets and $(params_part["num_pictures"]) images.")
                else
                    println("\t------ERROR: $histofile.")
                end
            catch err
                println("\t------ERROR: $histofile.\n$err")
            end
        end
    end

    println("Writing to $(output_dir)/histo.dat")
    try mkdir(output_dir) end
    serializeToFile("$(output_dir)/histo.dat", (params, c2_full, c3_full, c1_full))
end

function create_exp_filelist(name::String="coliphage")
    root = environment_path("exp_data/parts/$(name)")
    list = readdir(root)
    map((p)->"$root/$p/histo.dat", list)
end

function process_exp_data(name::String="coliphage", beamstop::String="coliphage_beamstop_K2_38_K3_30_N_32", Gauss_filter::Bool=false)
    combine_histograms(create_exp_filelist(name), environment_path("exp_data/$(name)"))
    p,c2,c3,c1 = deserializeFromFile(environment_path("exp_data/$(name)/histo.dat"))
    _,c2_beamstop,c3_beamstop,_ = deserializeFromFile(environment_path("exp_data/$(beamstop)/histo.dat"))
    c2_filtered,c3_filtered = postprocess_correlations(c2, c3, c2_beamstop, c3_beamstop, Gauss_filter)
    newname = "$(name)_processed$(Gauss_filter ? "_smoothed" : "_notsmoothed")"
    try mkdir(environment_path("exp_data/$(newname)")) end
    serializeToFile(environment_path("exp_data/$(newname)/histo.dat"), (p, c2_filtered, c3_filtered, c1))
end

function combine_set_noise(img::Int64, setsize::Int64, sigmavals::Vector{Float64}=Float64[0.5, 0.75, 1.125], gammavals::Vector{Float64}=[0.1, 0.2, 0.3, 0.4, 0.5], ppi::Int64=10, K::Int64=38, N::Int64=32, lambda::Float64=0.0)
    for sigma in sigmavals
        for gamma in gammavals
            name = histogram_name("parallel/data_generation/parts/SH_", ppi, N, K, float(K), setsize, "", gamma, sigma, lambda)
            combine_histograms(name, img, K, N, setsize)
        end
    end
end

#Combine noisy histograms
# combine_set_noise(Integer(3.2768e9), Integer(2*2.048e7), [2.5], [0.4, 0.5], 10, 38, 32)

function combine_set(images::Array{Int64}, setsize::Int64, ppi::Int64=10, K::Int64=38, N::Int64=32, lambda::Float64=0.0, prefix::String="$(ENV["DETERMINATION_DATA"])/data_generation/parts/SH_")
    for pic in images
        name = histogram_name(prefix, ppi, N, K, float(K), setsize, "", 0.0, 0.0, lambda)
        combine_histograms(name, pic, K, N, setsize)
    end
end

# generate_single_multiparticle_histogram(Integer(3.2768e9), Integer(2.048e7), 2)
# combine_histograms(environment_path("data_generation/parts/multi_2_SH_10p_N32_K2_38_K3_26_R38.0_P40960000"), 80)



# generate_histogram_image(Integer(3.2768e9), 10, 38, 26, 32, 2.5; setsize=Integer(2*2.048e7), name="Ewald_large_")
# generate_histogram_image(Integer(3.2768e8), 10, 38, 26, 32, 2.5; setsize=Integer(2*2.048e7), name="Ewald_medium_")
# combine_histograms(environment_path("data_generation/parts/Ewald_lambda_2.0_SH_10p_N32_K2_38_K3_26_R38.0_P40960000"), 80)

#With beamstop
# generate_histogram_image(Integer(3.2768e9), 10, 38, 26, 32; setsize=Integer(2*2.048e7), name="Ewald_lambda_0.0_beamstop_", lambda=0.0, beamstop_width=qmax(38,38.0)/20.0)
# combine_histograms(environment_path("data_generation/parts/Ewald_lambda_0.0_beamstop_SH_10p_N32_K2_38_K3_26_R38.0_P4096000"), 80)

#beamstop only
# generate_histograms(; max_pictures = Integer(2*2.048e7), Ncores=8, N=32, number_incident_photons=100, batchsize = Integer(2.048e6), successive_jobs=1, prefix="beamstop_only_3_", suffix="", use_cube=false, qcut_ratio=1.0, K2=38, K3=26, rmax=float(38), histogram_method="histogramCorrelationsInPicture_alltoall", structure_pdb_path="$(ENV["DETERMINATION_DATA"])/structures/crambin.pdb", lambda=0.0, beamstop_width=qmax(38,38.0)/20.0, beamstop_only=true)




#Coliphage
# combine_histograms(environment_path("exp_data/parts/coliphage_symmetric"), 42)
