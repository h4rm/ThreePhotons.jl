#Correlations

#basis for three-photon correlation
basis = complexBasis(L,N,LMAX)
#two-photon correlation
c2 = twoPhotons(intensity, basis, K)
#dummy reference c3
c3ref = zeros(N,N,K,K,K)
#energy calculation including three-photon correlation
energy(intensity, basis, c3ref, K)
