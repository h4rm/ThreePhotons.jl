using CUDArt
using CUBLAS
global CUDA_enabled = false

#cuda
export BasisTypeCuda, CUDA_store_basis, CUDA_init, FullCorrelation_parallized

function CUDA_init()
    CUDArt.init(0)
    CUDArt.device(0)
    global md = CuModule("$(ENV["THREEPHOTONS_PATH"])/src/cuda_kernel.ptx", false)
    global calculate_coefficient_matrix_cuda = CuFunction(md, "calculate_coefficient_matrix")
    global CUDA_enabled = true
    println("Initialization of CUDA complete.")
end

"""Datatype for precalculated three photon correlation basis function"""
type BasisTypeCuda <: AbstractBasisType
    d_wignerlist::CudaArray
    d_indices::CudaArray
    d_PAcombos::CudaArray
    B::CudaArray
    P::SharedArray{Float64,2}
    d_correlation::CudaArray

    basislen::Int64
    N::Int64
    L::Int64
    LMAX::Int64
    lrange::StepRange
    ctr::Dict
    rtc::Dict
    K::Int64
    lambda::Float64
    dq::Float64
end

"""Overrides the original complexBasis function"""
function CUDA_store_basis(basis::BasisType, L::Int64, LMAX::Int64, N::Int64, K::Int64, lambda::Float64, dq::Float64, forIntensity=true)
    if CUDA_enabled
        d_wignerlist = CudaArray(convert(Array{Float32}, sdata(basis.wignerlist)))
        d_indices = CudaArray(convert(Array{Int32}, sdata(basis.indices)))
        d_PAcombos = CudaArray(convert(Array{Int32}, sdata(basis.PAcombos)))
        d_B = CudaArray(convert(Array{Float32}, sdata(basis.B)))
        d_correlation = CudaArray(Float32,(basis.N^2,round(Int64, basis.K*(basis.K+1)*(basis.K+2)/6)))

        cuda_basis =  BasisTypeCuda(d_wignerlist, d_indices, d_PAcombos, d_B, basis.P, d_correlation, basis.basislen, basis.N, basis.L, basis.LMAX, basis.lrange, basis.ctr, basis.rtc, basis.K, basis.lambda, basis.dq)
        println("Initialized CUDA basis.")
        return cuda_basis
    else
        return basis
    end
end

function CUDA_store_basis( L::Int64, LMAX::Int64, N::Int64, K::Int64, lambda::Float64, dq::Float64, forIntensity=true)
    basis = calculate_basis(L, LMAX, N, K, lambda, dq, forIntensity)
    if CUDA_enabled
        return CUDA_store_basis(basis, L, LMAX, N, K, lambda, dq, forIntensity)
    else
        return basis
    end
end


complexBasis_choice = CUDA_store_basis

"""Calculates the full three photon correlation"""
function FullCorrelation_parallized(intensity::SphericalHarmonicsVolume, basis::BasisTypeCuda, minimal::Bool=true, normalize::Bool=false, return_raw::Bool=false)

    #Prepare all arrays on GPU
    klength = Integer(basis.K*(basis.K+1)*(basis.K+2)/6)
    numcoeff = num_coeff(basis.LMAX)
    d_coeff = CudaArray(Complex{Float32}[intensity.coeff[k][i] for k = 1:basis.K, i=1:numcoeff]')
    d_PA = CudaArray(convert(Array{Float32}, basis.P))

    threadsperblock = 1024
    blockspergrid = ceil(Int32, Base.size(basis.d_PAcombos)[2] / threadsperblock)
    launch(calculate_coefficient_matrix_cuda, blockspergrid, threadsperblock, (d_coeff, numcoeff, basis.d_wignerlist, basis.d_indices, Base.size(basis.d_indices)[2], basis.d_PAcombos, Base.size(basis.d_PAcombos)[2], d_PA, klength))
    println("Done with factor matrix")
    CUBLAS.gemm!('N','N',Float32(1.0), basis.B, d_PA, Float32(0.0), basis.d_correlation)

    #Transfer result to host
    res = to_host(basis.d_correlation)

    #Enforce GC to avoid crashes
    CUDArt.free(d_PA)
    CUDArt.free(d_coeff)
    if return_raw return res end

    #Translate matrix into 5dimensional correlation format
    t = zeros(Float64, basis.N, basis.N, K, K, K)
    i = 1
    for k1=1:K for k2=1:k1 for k3=1:k2
        t[:,:, k3, k2, k1] = reshape(res[:, i], basis.N, basis.N)
        i += 1
    end end end
    return normalize ? t / sumabs(t) : t
end

# "Calculates the three photon correlation from spherical harmonics coefficients in a paralllized way."
# function FullCorrelation_parallized_combined(intensity::SphericalHarmonicsVolume, basis_cpu::BasisType, basis_gpu::BasisTypeCuda, K::Int64=8, minimal::Bool=true, normalize::Bool=false)
#     kcombinations = Tuple{Int64,Int64,Int64}[(k1,k2,k3) for k1 = 1:K for k2 = 1:(minimal ? k1 : K) for k3 = 1:(minimal ? k2 : K)]
#
#     #Prepare CPU arrays
#     c = SharedArray(Float64,(basis_cpu.N,basis_cpu.N,K,K,K))
#
#     #Prepare all arrays on GPU
#     numcoeff = num_coeff(intensity.LMAX)
#     KMAX = intensity.KMAX
#     d_coeff = CudaArray(convert(Array{Complex{Float32}}, Complex{Float32}[intensity.coeff[k][i] for k = 1:KMAX, i=1:numcoeff]'))
#     d_faclist = CudaArray(zeros(Float32,basis_gpu.basislen))
#     d_correlation = CudaArray(zeros(Float32, K^3*basis_gpu.N^2))
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
#               calculateCorrelationSlices(d_correlation, d_faclist, d_coeff, basis_gpu, numcoeff, K, KMAX, k1,k2,k3)
#             end
#           end
#         end
#       end
#     end
#
#     #Combine CPU and GPU results
#     res_cpu = sdata(c)
#     res_gpu = Array{Float64}(reshape(to_host(d_correlation), basis_gpu.N, basis_gpu.N, K, K, K))
#     res = res_cpu + res_gpu
#     return normalize ? res / sumabs(res) : res
# end
