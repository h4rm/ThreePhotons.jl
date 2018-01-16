using CUDAdrv
using CUDArt
include("CUBLAS.jl-0.0.2/src/CUBLAS.jl")
# using "CUBLAS.jl-0.0.2/src/CUBLAS"
CUDA_enabled = false

#cuda
export BasisTypeCuda, CUDA_store_basis, CUDA_calculate_basis, complexBasis_choice, CUDA_init, FullCorrelation_parallized, new_energy

function CUDA_init()
    CUBLAS_init()
    CUDArt.init(0)
    CUDArt.device(0)
    global md = CUDAdrv.CuModuleFile("$(ENV["THREEPHOTONS_PATH"])/src/cuda_kernel.ptx")
    global calculate_coefficient_matrix_cuda = CUDAdrv.CuFunction(md, "calculate_coefficient_matrix")
    global energy_cuda = CUDAdrv.CuFunction(md, "energy")

    global CUDA_enabled = true
    CUDAdrv.gc()
    println("Initialization of CUDA complete.")
end

"""Datatype for precalculated three photon correlation basis function"""
type BasisTypeCuda <: AbstractBasisType
    d_wignerlist::CudaArray
    d_indices::CudaArray
    d_PAcombos::CudaArray
    d_B::CudaArray
    h_P::Array{Float32,2}
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
function CUDA_store_basis(basis::BasisType)
    d_wignerlist = CudaArray(convert(Array{Float32}, sdata(basis.wignerlist)))
    d_indices = CudaArray(convert(Array{Int32}, sdata(basis.indices)))
    d_PAcombos = CudaArray(convert(Array{Int32}, sdata(basis.PAcombos)))
    d_B = CudaArray(sdata(basis.B))
    d_correlation = CudaArray(Float32,(2*basis.N^2,round(Int64, basis.K*(basis.K+1)*(basis.K+2)/6)))

    cuda_basis =  BasisTypeCuda(d_wignerlist, d_indices, d_PAcombos, d_B, basis.h_P, d_correlation, basis.basislen, basis.N, basis.L, basis.LMAX, basis.lrange, basis.ctr, basis.rtc, basis.K, basis.lambda, basis.dq)
    println("Initialized CUDA basis.")
    return cuda_basis
end

function CUDA_calculate_basis( L::Int64, LMAX::Int64, N::Int64, K::Int64, lambda::Float64, dq::Float64, forIntensity=true)
    println("Intialize basis for CUDA.")
    basis = calculate_basis(L, LMAX, N, K, lambda, dq, forIntensity)
    if CUDA_enabled
        return CUDA_store_basis(basis)
    else
        return basis
    end
end

function CUDA_calculate_basis( L::Int64, LMAX::Int64, N::Int64, K::Int64, forIntensity=true)
    return CUDA_calculate_basis(L,LMAX,N,K, 0.0, 0.0, forIntensity)
end

global complexBasis_choice = CUDA_calculate_basis

"""Calculates the full three photon correlation"""
function FullCorrelation_parallized(intensity::SphericalHarmonicsVolume, basis::BasisTypeCuda, minimal::Bool=true, normalize::Bool=false, return_raw::Bool=false)
    # println("CUDA correlation.")
    #Prepare all arrays on GPU
    klength = Integer(basis.K*(basis.K+1)*(basis.K+2)/6)
    numcoeff = num_coeff(basis.LMAX)
    d_coeff = CudaArray(Complex{Float32}[intensity.coeff[k][i] for k = 1:basis.K, i=1:numcoeff]')
    d_PA = CudaArray(basis.h_P)

    threadsperblock = 1024
    blockspergrid = ceil(Int32, Base.size(basis.d_PAcombos)[2] / threadsperblock)
    CUDAdrv.cudacall(calculate_coefficient_matrix_cuda, blockspergrid, threadsperblock,
    (Ptr{Cfloat}, Cint, Ptr{Cfloat}, Ptr{Cint}, Cint, Ptr{Cint}, Cint, Ptr{Cfloat}, Cint),
    d_coeff, numcoeff, basis.d_wignerlist, basis.d_indices, Base.size(basis.d_indices)[2], basis.d_PAcombos, Base.size(basis.d_PAcombos)[2], d_PA, klength)

    gemm!('N','N',Float32(1.0), basis.d_B, d_PA, Float32(0.0), basis.d_correlation)

    #Wait for computation to finish and transfer result to host
    res = to_host(basis.d_correlation)

    #Enforce GC to avoid crashes
    CUDArt.free(d_PA)
    CUDArt.free(d_coeff)
    # CUDAdrv.gc()

    if return_raw return res end

    #Translate matrix into 5dimensional correlation format
    t = zeros(Float64, basis.N, 2*basis.N, basis.K, basis.K, basis.K)
    i = 1
    for k1=1:basis.K for k2=1:k1 for k3=1:k2
        t[:,:, k3, k2, k1] = reshape(res[:, i], basis.N, 2*basis.N)
        i += 1
    end end end
    return normalize ? t / sumabs(t) : t
end

# """Calculate the energy on the GPU to avoid data transfer of large correlation Matrix
# However, this function is not much faster than the original.
# Bottlenck is still the large matrix operation."""
# function new_energy(intensity::SphericalHarmonicsVolume, basis::BasisTypeCuda, c3ref::C3, K3_range::UnitRange{Int64}, measure::String="Bayes", negativity_factor::Float64=0.0)
#     # println("CUDA correlation.")
#     #Prepare all arrays on GPU
#     klength = Integer(basis.K*(basis.K+1)*(basis.K+2)/6)
#     numcoeff = num_coeff(basis.LMAX)
#     d_coeff = CudaArray(Complex{Float32}[intensity.coeff[k][i] for k = 1:basis.K, i=1:numcoeff]')
#     d_PA = CudaArray(basis.h_P)
#
#     threadsperblock = 1024
#     blockspergrid = ceil(Int32, Base.size(basis.d_PAcombos)[2] / threadsperblock)
#     CUDAdrv.cudacall(calculate_coefficient_matrix_cuda, blockspergrid, threadsperblock, (d_coeff, numcoeff, basis.d_wignerlist, basis.d_indices, Base.size(basis.d_indices)[2], basis.d_PAcombos, Base.size(basis.d_PAcombos)[2], d_PA, klength))
#
#     gemm!('N','N',Float32(1.0), basis.d_B, d_PA, Float32(0.0), basis.d_correlation)
#
#     Kcombos = [(k1,k2,k3) for k1 in K3_range for k2=minimum(K3_range):k1 for k3=minimum(K3_range):k2]
#     Kcombos_length = length(Kcombos)
#     Kcombos_matrix = Int32[Kcombos[i][j] for i = 1:length(Kcombos), j=1:3]
#
#     d_Kcombos = CudaArray(Kcombos_matrix')
#     d_result = CudaArray(Float32, Kcombos_length)
#     c3ref_compressed = Array(Float32, 2*basis.N^2, Kcombos_length)
#
#     for i=1:Kcombos_length
#         k1,k2,k3 = Kcombos[i]
#         c3ref_compressed[:,i] = reshape(c3ref[:,:, k3,k2,k1], 2*basis.N^2)
#     end
#     d_c3ref = CudaArray(c3ref_compressed')
#
#     #  __global__ void energy(const float *c3, const float *c3ref, const int* Kcombos, const int combolength, const int slab_size, float *result)
#     threadsperblock = 1024
#     blockspergrid = ceil(Int32, Kcombos_length/threadsperblock)
#
#     launch(energy_cuda, blockspergrid, threadsperblock, (basis.d_correlation, d_c3ref, d_Kcombos, Kcombos_length, 2*basis.N^2, d_result))
#
#     return sumabs(to_host(d_result))
# end
