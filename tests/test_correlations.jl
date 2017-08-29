#Correlations

#basis for three-photon correlation
basis = calculate_basis(L,LMAX,N,K)
#two-photon correlation
println("Calculate two photon correlation.")
c2 = twoPhotons(intensity, basis, K)
#dummy reference c3
c3ref = zeros(N,N,K,K,K)
println("calculate energy.")
#energy calculation including three-photon correlation
energy(intensity, basis, c3ref)
