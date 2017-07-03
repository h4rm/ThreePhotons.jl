using CUDArt
using CUBLAS
global CUDA_enabled = false

function CUDA_init()
  CUDArt.init(0)
  CUDArt.device(0)
  global md = CuModule("$(ENV["DETERMINATION_PATH"])/src/cuda_kernel.ptx", false)
  global calculate_coefficient_matrix_cuda = CuFunction(md, "calculate_coefficient_matrix")
  global calculate_coefficient_matrix_optimized_cuda = CuFunction(md, "calculate_coefficient_matrix_optimized")
  global CUDA_enabled = true
  println("Initialization of CUDA complete.")
end

"""Datatype for precalculated three photon correlation basis function"""
type BasisTypeCuda <: AbstractBasisType
  basis::CudaArray
  basisindices::CudaArray
  basislen::Int64
  N::Int64
  lcut::Int64
  lmax::Int64
  lrange::StepRange
  ctr::Dict
  rtc::Dict
end

"""Overrides the original complexBasis function"""
function complexBasis_CUDA(basis::BasisType, lcut::Int64, N::Int64, lmax::Int64, forIntensity=true)
  if CUDA_enabled
    d_basis = CudaArray(convert(Array{Float32}, sdata(basis.basis)))
    d_basisindices = CudaArray(convert(Array{Int32}, sdata(basis.basisindices)))
    cuda_basis =  BasisTypeCuda(d_basis, d_basisindices, basis.basislen, basis.N, basis.lcut, basis.lmax, basis.lrange, basis.ctr, basis.rtc)
    println("Initialized CUDA basis.")
    return cuda_basis
  else
    return basis
  end
end

"""Overrides the original complexBasis function"""
function complexBasis_CUDA(lcut::Int64, N::Int64, lmax::Int64, forIntensity=true)
  basis = complexBasis(lcut, N, lmax, forIntensity)
  if CUDA_enabled
    d_basis = CudaArray(convert(Array{Float32}, sdata(basis.basis)))
    d_basisindices = CudaArray(convert(Array{Int32}, sdata(basis.basisindices)))
    cuda_basis =  BasisTypeCuda(d_basis, d_basisindices, basis.basislen, basis.N, basis.lcut, basis.lmax, basis.lrange, basis.ctr, basis.rtc)
    println("Initialized CUDA basis.")
    return cuda_basis
  else
    return basis
  end
end

complexBasis_choice = complexBasis_CUDA

"""Calculates the full three photon correlation"""
function FullCorrelation_parallized(intensity::SphericalHarmonicsVolume, basis::BasisTypeCuda, kcut::Int64, minimal::Bool=true, normalize::Bool=false, return_raw::Bool=false)

    #Prepare all arrays on GPU
    numcoeff = num_coeff(intensity.lmax)
    kmax = intensity.kmax
    d_coeff = CudaArray(Complex{Float32}[intensity.coeff[k][i] for k = 1:kmax, i=1:numcoeff]')

    klength = round(Int64, kcut*(kcut+1)*(kcut+2)/6)
    d_faclist = CudaArray(Float32,(basis.basislen, klength))
    d_correlation = CudaArray(Float32,(basis.N^2,klength))

    # Set up grid and block, see below
    # kmapping = [(k1,k2,k3) for k1=1:kcut for k2=1:k1 for k3=1:k2]
    # d_kmapping = CudaArray(Int32[kmapping[i][j] for i=1:klength,j=1:3])
    # threadsperblock = (32,32)
    # blockspergrid = ( ceil( Int32, basis.basislen / threadsperblock[1]), ceil(Int32, klength / threadsperblock[2]) )
    # launch(calculate_coefficient_matrix_optimized_cuda, blockspergrid, threadsperblock, (d_coeff, d_kmapping, numcoeff, kcut, basis.basisindices, basis.basislen, d_faclist, klength))

    threadsperblock = 1024
    blockspergrid = ceil(Int32, basis.basislen / threadsperblock)
    launch(calculate_coefficient_matrix_cuda, blockspergrid, threadsperblock, (d_coeff, numcoeff, kcut, basis.basisindices, basis.basislen, d_faclist))

    CUBLAS.gemm!('N','N',Float32(1.0), basis.basis, d_faclist, Float32(0.0), d_correlation)

    #Transfer result to host
    res = to_host(d_correlation)

    #Enforce GC to avoid crashes
    CUDArt.free(d_faclist)
    CUDArt.free(d_correlation)
    CUDArt.free(d_coeff)
    if return_raw return res end

    #Translate matrix into 5dimensional correlation format
    t = zeros(Float64, basis.N, basis.N, kcut, kcut, kcut)
    i = 1
    for k1=1:kcut for k2=1:k1 for k3=1:k2
      t[:,:, k3, k2, k1] = reshape(res[:, i], basis.N, basis.N)
      i += 1
    end end end
    return normalize ? t / sumabs(t) : t
end

# "Calculates the three photon correlation from spherical harmonics coefficients in a paralllized way."
# function FullCorrelation_parallized_combined(intensity::SphericalHarmonicsVolume, basis_cpu::BasisType, basis_gpu::BasisTypeCuda, kcut::Int64=8, minimal::Bool=true, normalize::Bool=false)
#     kcombinations = Tuple{Int64,Int64,Int64}[(k1,k2,k3) for k1 = 1:kcut for k2 = 1:(minimal ? k1 : kcut) for k3 = 1:(minimal ? k2 : kcut)]
#
#     #Prepare CPU arrays
#     c = SharedArray(Float64,(basis_cpu.N,basis_cpu.N,kcut,kcut,kcut))
#
#     #Prepare all arrays on GPU
#     numcoeff = num_coeff(intensity.lmax)
#     kmax = intensity.kmax
#     d_coeff = CudaArray(convert(Array{Complex{Float32}}, Complex{Float32}[intensity.coeff[k][i] for k = 1:kmax, i=1:numcoeff]'))
#     d_faclist = CudaArray(zeros(Float32,basis_gpu.basislen))
#     d_correlation = CudaArray(zeros(Float32, kcut^3*basis_gpu.N^2))
#     stream = Stream()
#
#     np = nprocs()  # determine the number of processes available
#     n = length(kcombinations)
#
#     i = 1
#     # function to produce the next work item from the queue.
#     nextidx() = (idx=i; i+=1; idx)
#
#     @sync begin
#       #Go over all CPUs + GPU
#       for p=2:np+1
#         #Start infinite scheduling for CPU + GPU
#         @async begin
#           while true
#             #Fetch next available combination
#             idx = nextidx()
#             if idx > n
#               break
#             end
#             combo = kcombinations[idx]
#             k1,k2,k3 = combo[1], combo[2], combo[3]
#             #feed cpu job
#             if p >= 2 && p <= np
#               # println("Launching cpu job on $p ($k1,$k2,$k3)")
#               remotecall_wait(calculateCorrelationSlices!, p, c, intensity.coeff, basis_cpu, k1,k2,k3)
#             #feed gpu job
#             else
#               # println("Launching gpu job ($k1,$k2,$k3)")
#               wait(stream)
#               calculateCorrelationSlices(d_correlation, d_faclist, d_coeff, basis_gpu, numcoeff, kcut, kmax, k1,k2,k3)
#             end
#           end
#         end
#       end
#     end
#
#     #Combine CPU and GPU results
#     res_cpu = sdata(c)
#     res_gpu = Array{Float64}(reshape(to_host(d_correlation), basis_gpu.N, basis_gpu.N, kcut, kcut, kcut))
#     res = res_cpu + res_gpu
#     return normalize ? res / sumabs(res) : res
# end
