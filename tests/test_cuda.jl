CUDA_init()
basis = complexBasis(L,N,LMAX)
d_basis = complexBasis_CUDA(basis, L,N,LMAX)
c3ref = ones(N,N,K,K,K)

#Cuda
println("Calculating CUDA version")
c_cuda = FullCorrelation_parallized(intensity, d_basis, K)

#Normal
println("Calculating CPU version")
c_normal = FullCorrelation_parallized(intensity, basis, K)

#Are they the same?
@test sumabs(c_normal-c_cuda)/sumabs(c_normal) < 1e-5

#For performance tests
for i=1:5
  println("Cuda:")
  @time FullCorrelation_parallized(intensity, d_basis, K)
  println("Normal:")
  @time FullCorrelation_parallized(intensity, basis, K)
  println("-----")
end

println("CUDA energy:")
e_cuda = energy(intensity, d_basis, c3ref, K, "Bayes")
println(e_cuda)

println("CPU energy:")
e_cpu = energy(intensity, basis, c3ref, K, "Bayes")
println(e_cpu)

@test abs(e_cuda - e_cpu)/abs(e_cuda) < 1e-6

for i=1:5
  println("CUDA energy:")
  @time energy(intensity, d_basis, c3ref, K, "Bayes")
end
