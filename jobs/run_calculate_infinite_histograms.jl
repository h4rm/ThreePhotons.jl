addprocs()
using ThreePhotons

# using ThreePhotons
# include("$(ENV["THREEPHOTONS_PATH"])/src/cuda.jl")
# CUDA_init()

LMAX = 25
for N in [32]
  for L in [16]
    for K2 in [38]
        for K3 in [20]
          density,fourier,intensity = createSphericalHarmonicsStructure("$(ENV["DETERMINATION_DATA"])/structures/crambin.pdb", LMAX, K2, float(K2))

          basis = calculate_basis(L, LMAX, N, K3, 4.0, dq(intensity))

          c3 = FullCorrelation_parallized(intensity, basis, true)
          c2 = twoPhotons(intensity, basis, K2, true)
          serializeToFile("$(ENV["DETERMINATION_DATA"])/output/data_generation/correlations_N$(N)_K2_$(K2)_K3_$(K3)_L$(L)_inf.dat", (Dict("num_pictures"=>"inf"),c2,c3))
      end
    end
  end
end
