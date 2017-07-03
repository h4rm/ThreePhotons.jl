using Base.Test
using ThreePhotons

#Postprocess runs

params = Dict("L"=>L, "K" =>K, "N"=>N, "LMAX"=>LMAX, "KMAX"=>KMAX, "rmax"=>float(KMAX))
state = Dict{Any,Any}("intensity"=>deleteTerms(intensity, K, L))
postprocess_run(params, state, "../data/structures/crambin.pdb", false)
