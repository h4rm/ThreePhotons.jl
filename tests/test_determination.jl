using Base.Test
using ThreePhotons

# Determination

params,state = rotation_search(Dict("reference_pdb_path"=>"../data/structures/crambin.pdb","stepsizefactor"=>1.01, "initial_stepsize" => Float64(pi/100), "L"=>L, "K" =>K, "N"=>N, "histograms"=>"/tmp/correlations.dat", "optimizer"=>rotate_hierarchical, "initial_temperature_factor"=>1.0e-6, "measure"=>"Bayes", "temperature_decay"=>0.96, "LMAX"=>LMAX, "KMAX"=>KMAX, "rmax"=>float(KMAX)))
