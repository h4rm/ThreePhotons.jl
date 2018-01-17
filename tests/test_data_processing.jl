using Base.Test
using ThreePhotons

#Postprocess runs

params = Dict("L"=>L, "K3_range" =>1:K, "N"=>N, "LMAX"=>LMAX, "K2_range"=>1:KMAX, "rmax"=>float(KMAX), "qmax"=>float(KMAX))
state = Dict{Any,Any}("intensity"=>deleteTerms(intensity, K, L))
postprocess_run(params, state, "../data/structures/crambin.pdb", false)
