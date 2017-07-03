using Base.Test
using ThreePhotons

#Postprocess runs

params = Dict("lcut"=>L, "kcut" =>K, "N"=>N, "lmax"=>LMAX, "kmax"=>KMAX, "rmax"=>float(KMAX))
state = Dict{Any,Any}("intensity"=>deleteTerms(intensity, K, L))
postprocess_run(params, state, false)
