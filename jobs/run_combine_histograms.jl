using ThreePhotons
include("runs.jl")

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

function combine_set_noise(img::Int64, setsize::Int64, sigmavals::Vector{Float64}=Float64[0.5, 0.75, 1.125], gammavals::Vector{Float64}=[0.1, 0.2, 0.3, 0.4, 0.5], ppi::Int64=10, K::Int64=38, N::Int64=32)
  for sigma in sigmavals
    for gamma in gammavals
      name = histogram_name("parallel/data_generation/parts/SH_", ppi, N, K, float(K), setsize, "", gamma, sigma)
      combine_histograms(name, img, K, N, setsize)
    end
  end
end

#Combine noisy histograms
# combine_set_noise(Integer(3.2768e9), Integer(2*2.048e7), [2.5], [0.4, 0.5], 10, 38, 32)

function combine_set(images::Array{Int64}, setsize::Int64, ppi::Int64=10, K::Int64=38, N::Int64=32, prefix::String="$(ENV["DETERMINATION_DATA"])/data_generation/parts/SH_")
  for pic in images
    name = histogram_name(prefix, ppi, N, K, float(K), setsize, "")
    combine_histograms(name, pic, K, N, setsize)
  end
end
