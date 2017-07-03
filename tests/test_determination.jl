using Base.Test
using ThreePhotons

# Determination

params,state = rotation_search(Dict("stepsizefactor"=>1.02, "initial_stepsize" => Float64(pi), "iterations"=>1, "lcut"=>4, "kcut" =>4, "N"=>32, "histograms"=>"$(ENV["DETERMINATION_PATH"])/expdata/correlations_N32_K35_L20_inf.dat", "optimizer"=>rotate_hierarchical, "initial_temperature_factor"=>1.0, "measure"=>"Bayes", "temperature_decay"=>0.95, "lmax"=>25, "kmax"=>35, "rmax"=>35.0))
end

c2all,c2, c3all, c3 = loadHistograms(kcut, "$(ENV["DETERMINATION_PATH"])/expdata/correlations_N48_K25_P2048000.dat")
retrievedIntensity = retrieveSolution(c2all, lcut, intensity.lmax, intensity.kmax, intensity.rmax)
basis = complexBasis(lcut,N,intensity.lmax)
rotatedIntensity,FSC = checkRotationSearch(retrievedIntensity, deleteTerms(intensity, kcut, lcut), kcut, lcut, basis, iterations=Integer(1.0e4), reduce_stepsize=1000, energy = (volume)->energy(volume, basis, c3, kcut), save_structures=false)
@test FSC > 0.8
