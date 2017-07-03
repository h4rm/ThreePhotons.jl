addprocs()
@everywhere using ThreePhotons

# using ThreePhotons
# include("$(ENV["THREEPHOTONS_PATH"])/src/cuda.jl")
# CUDA_init()


for N in [48]
  for L in [20]
    basis = complexBasis(L,N,25)
    for K in [35, 38]
      density,fourier,intensity = createSphericalHarmonicsStructure("$(ENV["THREEPHOTONS_PATH"])/data/structures/crambin.pdb", 25, K, float(K))
      c3 = FullCorrelation_parallized_piecewise(intensity, basis, K, true)
      c2 = twoPhotons(intensity, basis, K, true)
      serializeToFile("$(ENV["THREEPHOTONS_PATH"])/expdata/correlations_N$(N)_K$(K)_L$(L)_inf.dat", (Dict("num_pictures"=>"inf"),c2,c3))
    end
  end
end
