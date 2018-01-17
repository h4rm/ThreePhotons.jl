CUDA_init()
basis = calculate_basis(L,LMAX,N,K, 2.5, dq(intensity))
d_basis = CUDA_store_basis(basis)
c3ref = ones(N,2*N,K,K,K)

#Cuda
println("Calculating CUDA version")
c_cuda = FullCorrelation_parallized(intensity, d_basis)

# #Normal
# println("Calculating CPU version")
# c_normal = FullCorrelation_parallized(intensity, basis)
#
# #Are they the same?
# @test sum(abs, c_normal-c_cuda)/sum(abs, c_normal) < 1e-5

# #For performance tests
# for i=1:5
#   println("Cuda:")
#   @time FullCorrelation_parallized(intensity, d_basis)
#   println("Normal:")
#   @time FullCorrelation_parallized(intensity, basis)
#   println("-----")
# end

println("CUDA energy:")
e_cuda = energy(intensity, d_basis, c3ref, "Bayes")
println(e_cuda)

println("CPU energy:")
e_cpu = energy(intensity, basis, c3ref, "Bayes")
println(e_cpu)

@test abs(e_cuda - e_cpu)/abs(e_cuda) < 1e-6

for i=1:5
    println("CUDA energy:")
    @time e_cuda = energy(intensity, d_basis, c3ref, "Bayes")

    println("CPU CUDA energy:")
    @time e_cpu = energy(intensity, basis, c3ref, "Bayes")
end
