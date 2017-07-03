using ThreePhotons
include("runs.jl")

function combine_histograms(dir::String, pic::Int64, K::Int64, N::Int64, batchsize::Int64=20480000)
  println("Processing $pic")
  c2_full = zeros(N,K,K)
  c3_full = zeros(N,N,K,K,K)
  params,_,_ = deserializeFromFile("$(dir)_1/histo.dat")

  for i = 1:ceil(Int64,pic/batchsize)
    println("\tLoading $i. file.")
    histofile = "$(dir)_$(i)/histo.dat"
    if isfile(histofile)
      _, c2_full_part, c3_full_part =  deserializeFromFile(histofile)
      part_triplets = sum(c3_full_part)
      if part_triplets < float(1.0e20)
        c2_full += c2_full_part
        c3_full += c3_full_part
        println("\tAdded $(part_triplets) triplets.")
      else
        println("\t------ERROR: $i file.")
      end
    end
  end

  params["num_pictures"] = pic
  finalpath = replace(replace(dir, "parts/", ""), "P$(batchsize)", "P$(pic)")

  println("Writing to $(finalpath).")
  try mkdir(finalpath) end
  serializeToFile("$(finalpath)/histo.dat", (params, c2_full, c3_full))
end

function combine_set_noise(img::Int64, setsize::Int64, sigmavals::Vector{Float64}=Float64[0.5, 0.75, 1.125], gammavals::Vector{Float64}=[0.1, 0.2, 0.3, 0.4, 0.5], ppi::Int64=10, K::Int64=38, N::Int64=32)
  for sigma in sigmavals
    for gamma in gammavals
      name = histogram_name("parallel/data_generation/parts/SH_", ppi, N, K, float(K), setsize, "", gamma, sigma)
      combine_histograms(name, img, K, N, setsize)
    end
  end
end

#Combine noisy histograms
combine_set_noise(Integer(3.2768e9), Integer(2*2.048e7), [2.5], [0.4, 0.5], 10, 38, 32)

function combine_set(images::Array{Int64}, setsize::Int64, ppi::Int64=10, K::Int64=38, N::Int64=32)
  for pic in images
    name = histogram_name("parallel/data_generation/parts/SH_", ppi, N, K, float(K), setsize, "")
    combine_histograms(name, pic, K, N, setsize)
  end
end

#Combine non noisy histograms
# combine_set(calculate_images_ppi(25)[5:6], calculate_images_ppi(25)[4], 25)
# combine_set(calculate_images_ppi(50)[5:6], calculate_images_ppi(50)[4], 50)
