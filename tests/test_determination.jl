using Base.Test
using ThreePhotons

# Determination

params,state = rotation_search(Dict("reference_pdb_path"=>"$(ENV["DETERMINATION_DATA"])/structures/crambin.pdb","stepsizefactor"=>1.01, "initial_stepsize" => Float64(pi/100), "L"=>L, "K3_range" =>1:K, "N"=>N, "histograms"=>"/tmp/correlations.dat", "optimizer"=>rotate_all_at_once, "initial_temperature_factor"=>1.0e-6, "measure"=>"Bayes", "temperature_decay"=>0.96, "LMAX"=>LMAX, "K2_range"=>1:KMAX, "qmax"=>float(KMAX), "lambda"=>0.0))
