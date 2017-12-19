# Julia wrapper for header: /usr/local/cuda/include/cublas_v2.h
# Automatically generated using Clang.jl wrap_c, version v0.0.1
# Manually copied from ../gen to this directory

function cublasCreate_v2(handle)
  statuscheck(ccall( (:cublasCreate_v2, libcublas), cublasStatus_t, (Ptr{cublasHandle_t},), handle))
end
function cublasDestroy_v2(handle)
  statuscheck(ccall( (:cublasDestroy_v2, libcublas), cublasStatus_t, (cublasHandle_t,), handle))
end
function cublasGetVersion_v2(handle, version)
  statuscheck(ccall( (:cublasGetVersion_v2, libcublas), cublasStatus_t, (cublasHandle_t, Ptr{Cint}), handle, version))
end
function cublasSetStream_v2(handle, streamId)
  statuscheck(ccall( (:cublasSetStream_v2, libcublas), cublasStatus_t, (cublasHandle_t, cudaStream_t), handle, streamId))
end
function cublasGetStream_v2(handle, streamId)
  statuscheck(ccall( (:cublasGetStream_v2, libcublas), cublasStatus_t, (cublasHandle_t, Ptr{cudaStream_t}), handle, streamId))
end
function cublasGetPointerMode_v2(handle, mode)
  statuscheck(ccall( (:cublasGetPointerMode_v2, libcublas), cublasStatus_t, (cublasHandle_t, Ptr{cublasPointerMode_t}), handle, mode))
end
function cublasSetPointerMode_v2(handle, mode)
  statuscheck(ccall( (:cublasSetPointerMode_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasPointerMode_t), handle, mode))
end
function cublasGetAtomicsMode(handle, mode)
  statuscheck(ccall( (:cublasGetAtomicsMode, libcublas), cublasStatus_t, (cublasHandle_t, Ptr{cublasAtomicsMode_t}), handle, mode))
end
function cublasSetAtomicsMode(handle, mode)
  statuscheck(ccall( (:cublasSetAtomicsMode, libcublas), cublasStatus_t, (cublasHandle_t, cublasAtomicsMode_t), handle, mode))
end
function cublasSetVector(n, elemSize, x, incx, devicePtr, incy)
  statuscheck(ccall( (:cublasSetVector, libcublas), cublasStatus_t, (Cint, Cint, Ptr{None}, Cint, Ptr{None}, Cint), n, elemSize, x, incx, devicePtr, incy))
end
function cublasGetVector(n, elemSize, x, incx, y, incy)
  statuscheck(ccall( (:cublasGetVector, libcublas), cublasStatus_t, (Cint, Cint, Ptr{None}, Cint, Ptr{None}, Cint), n, elemSize, x, incx, y, incy))
end
function cublasSetMatrix(rows, cols, elemSize, A, lda, B, ldb)
  statuscheck(ccall( (:cublasSetMatrix, libcublas), cublasStatus_t, (Cint, Cint, Cint, Ptr{None}, Cint, Ptr{None}, Cint), rows, cols, elemSize, A, lda, B, ldb))
end
function cublasGetMatrix(rows, cols, elemSize, A, lda, B, ldb)
  statuscheck(ccall( (:cublasGetMatrix, libcublas), cublasStatus_t, (Cint, Cint, Cint, Ptr{None}, Cint, Ptr{None}, Cint), rows, cols, elemSize, A, lda, B, ldb))
end
function cublasSetVectorAsync(n, elemSize, hostPtr, incx, devicePtr, incy, stream)
  statuscheck(ccall( (:cublasSetVectorAsync, libcublas), cublasStatus_t, (Cint, Cint, Ptr{None}, Cint, Ptr{None}, Cint, cudaStream_t), n, elemSize, hostPtr, incx, devicePtr, incy, stream))
end
function cublasGetVectorAsync(n, elemSize, devicePtr, incx, hostPtr, incy, stream)
  statuscheck(ccall( (:cublasGetVectorAsync, libcublas), cublasStatus_t, (Cint, Cint, Ptr{None}, Cint, Ptr{None}, Cint, cudaStream_t), n, elemSize, devicePtr, incx, hostPtr, incy, stream))
end
function cublasSetMatrixAsync(rows, cols, elemSize, A, lda, B, ldb, stream)
  statuscheck(ccall( (:cublasSetMatrixAsync, libcublas), cublasStatus_t, (Cint, Cint, Cint, Ptr{None}, Cint, Ptr{None}, Cint, cudaStream_t), rows, cols, elemSize, A, lda, B, ldb, stream))
end
function cublasGetMatrixAsync(rows, cols, elemSize, A, lda, B, ldb, stream)
  statuscheck(ccall( (:cublasGetMatrixAsync, libcublas), cublasStatus_t, (Cint, Cint, Cint, Ptr{None}, Cint, Ptr{None}, Cint, cudaStream_t), rows, cols, elemSize, A, lda, B, ldb, stream))
end
function cublasXerbla(srName, info)
  ccall( (:cublasXerbla, libcublas), None, (Ptr{UInt8}, Cint), srName, info)
end
function cublasSnrm2_v2(handle, n, x, incx, result)
  statuscheck(ccall( (:cublasSnrm2_v2, libcublas), cublasStatus_t, (cublasHandle_t, Cint, Ptr{Cfloat}, Cint, Ptr{Cfloat}), handle, n, x, incx, result))
end
function cublasDnrm2_v2(handle, n, x, incx, result)
  statuscheck(ccall( (:cublasDnrm2_v2, libcublas), cublasStatus_t, (cublasHandle_t, Cint, Ptr{Cdouble}, Cint, Ptr{Cdouble}), handle, n, x, incx, result))
end
function cublasScnrm2_v2(handle, n, x, incx, result)
  statuscheck(ccall( (:cublasScnrm2_v2, libcublas), cublasStatus_t, (cublasHandle_t, Cint, Ptr{cuComplex}, Cint, Ptr{Cfloat}), handle, n, x, incx, result))
end
function cublasDznrm2_v2(handle, n, x, incx, result)
  statuscheck(ccall( (:cublasDznrm2_v2, libcublas), cublasStatus_t, (cublasHandle_t, Cint, Ptr{cuDoubleComplex}, Cint, Ptr{Cdouble}), handle, n, x, incx, result))
end
function cublasSdot_v2(handle, n, x, incx, y, incy, result)
  statuscheck(ccall( (:cublasSdot_v2, libcublas), cublasStatus_t, (cublasHandle_t, Cint, Ptr{Cfloat}, Cint, Ptr{Cfloat}, Cint, Ptr{Cfloat}), handle, n, x, incx, y, incy, result))
end
function cublasDdot_v2(handle, n, x, incx, y, incy, result)
  statuscheck(ccall( (:cublasDdot_v2, libcublas), cublasStatus_t, (cublasHandle_t, Cint, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint, Ptr{Cdouble}), handle, n, x, incx, y, incy, result))
end
function cublasCdotu_v2(handle, n, x, incx, y, incy, result)
  statuscheck(ccall( (:cublasCdotu_v2, libcublas), cublasStatus_t, (cublasHandle_t, Cint, Ptr{cuComplex}, Cint, Ptr{cuComplex}, Cint, Ptr{cuComplex}), handle, n, x, incx, y, incy, result))
end
function cublasCdotc_v2(handle, n, x, incx, y, incy, result)
  statuscheck(ccall( (:cublasCdotc_v2, libcublas), cublasStatus_t, (cublasHandle_t, Cint, Ptr{cuComplex}, Cint, Ptr{cuComplex}, Cint, Ptr{cuComplex}), handle, n, x, incx, y, incy, result))
end
function cublasZdotu_v2(handle, n, x, incx, y, incy, result)
  statuscheck(ccall( (:cublasZdotu_v2, libcublas), cublasStatus_t, (cublasHandle_t, Cint, Ptr{cuDoubleComplex}, Cint, Ptr{cuDoubleComplex}, Cint, Ptr{cuDoubleComplex}), handle, n, x, incx, y, incy, result))
end
function cublasZdotc_v2(handle, n, x, incx, y, incy, result)
  statuscheck(ccall( (:cublasZdotc_v2, libcublas), cublasStatus_t, (cublasHandle_t, Cint, Ptr{cuDoubleComplex}, Cint, Ptr{cuDoubleComplex}, Cint, Ptr{cuDoubleComplex}), handle, n, x, incx, y, incy, result))
end
function cublasSscal_v2(handle, n, alpha, x, incx)
  statuscheck(ccall( (:cublasSscal_v2, libcublas), cublasStatus_t, (cublasHandle_t, Cint, Ptr{Cfloat}, Ptr{Cfloat}, Cint), handle, n, alpha, x, incx))
end
function cublasDscal_v2(handle, n, alpha, x, incx)
  statuscheck(ccall( (:cublasDscal_v2, libcublas), cublasStatus_t, (cublasHandle_t, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cint), handle, n, alpha, x, incx))
end
function cublasCscal_v2(handle, n, alpha, x, incx)
  statuscheck(ccall( (:cublasCscal_v2, libcublas), cublasStatus_t, (cublasHandle_t, Cint, Ptr{cuComplex}, Ptr{cuComplex}, Cint), handle, n, alpha, x, incx))
end
function cublasCsscal_v2(handle, n, alpha, x, incx)
  statuscheck(ccall( (:cublasCsscal_v2, libcublas), cublasStatus_t, (cublasHandle_t, Cint, Ptr{Cfloat}, Ptr{cuComplex}, Cint), handle, n, alpha, x, incx))
end
function cublasZscal_v2(handle, n, alpha, x, incx)
  statuscheck(ccall( (:cublasZscal_v2, libcublas), cublasStatus_t, (cublasHandle_t, Cint, Ptr{cuDoubleComplex}, Ptr{cuDoubleComplex}, Cint), handle, n, alpha, x, incx))
end
function cublasZdscal_v2(handle, n, alpha, x, incx)
  statuscheck(ccall( (:cublasZdscal_v2, libcublas), cublasStatus_t, (cublasHandle_t, Cint, Ptr{Cdouble}, Ptr{cuDoubleComplex}, Cint), handle, n, alpha, x, incx))
end
function cublasSaxpy_v2(handle, n, alpha, x, incx, y, incy)
  statuscheck(ccall( (:cublasSaxpy_v2, libcublas), cublasStatus_t, (cublasHandle_t, Cint, Ptr{Cfloat}, Ptr{Cfloat}, Cint, Ptr{Cfloat}, Cint), handle, n, alpha, x, incx, y, incy))
end
function cublasDaxpy_v2(handle, n, alpha, x, incx, y, incy)
  statuscheck(ccall( (:cublasDaxpy_v2, libcublas), cublasStatus_t, (cublasHandle_t, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint), handle, n, alpha, x, incx, y, incy))
end
function cublasCaxpy_v2(handle, n, alpha, x, incx, y, incy)
  statuscheck(ccall( (:cublasCaxpy_v2, libcublas), cublasStatus_t, (cublasHandle_t, Cint, Ptr{cuComplex}, Ptr{cuComplex}, Cint, Ptr{cuComplex}, Cint), handle, n, alpha, x, incx, y, incy))
end
function cublasZaxpy_v2(handle, n, alpha, x, incx, y, incy)
  statuscheck(ccall( (:cublasZaxpy_v2, libcublas), cublasStatus_t, (cublasHandle_t, Cint, Ptr{cuDoubleComplex}, Ptr{cuDoubleComplex}, Cint, Ptr{cuDoubleComplex}, Cint), handle, n, alpha, x, incx, y, incy))
end
function cublasScopy_v2(handle, n, x, incx, y, incy)
  statuscheck(ccall( (:cublasScopy_v2, libcublas), cublasStatus_t, (cublasHandle_t, Cint, Ptr{Cfloat}, Cint, Ptr{Cfloat}, Cint), handle, n, x, incx, y, incy))
end
function cublasDcopy_v2(handle, n, x, incx, y, incy)
  statuscheck(ccall( (:cublasDcopy_v2, libcublas), cublasStatus_t, (cublasHandle_t, Cint, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint), handle, n, x, incx, y, incy))
end
function cublasCcopy_v2(handle, n, x, incx, y, incy)
  statuscheck(ccall( (:cublasCcopy_v2, libcublas), cublasStatus_t, (cublasHandle_t, Cint, Ptr{cuComplex}, Cint, Ptr{cuComplex}, Cint), handle, n, x, incx, y, incy))
end
function cublasZcopy_v2(handle, n, x, incx, y, incy)
  statuscheck(ccall( (:cublasZcopy_v2, libcublas), cublasStatus_t, (cublasHandle_t, Cint, Ptr{cuDoubleComplex}, Cint, Ptr{cuDoubleComplex}, Cint), handle, n, x, incx, y, incy))
end
function cublasSswap_v2(handle, n, x, incx, y, incy)
  statuscheck(ccall( (:cublasSswap_v2, libcublas), cublasStatus_t, (cublasHandle_t, Cint, Ptr{Cfloat}, Cint, Ptr{Cfloat}, Cint), handle, n, x, incx, y, incy))
end
function cublasDswap_v2(handle, n, x, incx, y, incy)
  statuscheck(ccall( (:cublasDswap_v2, libcublas), cublasStatus_t, (cublasHandle_t, Cint, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint), handle, n, x, incx, y, incy))
end
function cublasCswap_v2(handle, n, x, incx, y, incy)
  statuscheck(ccall( (:cublasCswap_v2, libcublas), cublasStatus_t, (cublasHandle_t, Cint, Ptr{cuComplex}, Cint, Ptr{cuComplex}, Cint), handle, n, x, incx, y, incy))
end
function cublasZswap_v2(handle, n, x, incx, y, incy)
  statuscheck(ccall( (:cublasZswap_v2, libcublas), cublasStatus_t, (cublasHandle_t, Cint, Ptr{cuDoubleComplex}, Cint, Ptr{cuDoubleComplex}, Cint), handle, n, x, incx, y, incy))
end
function cublasIsamax_v2(handle, n, x, incx, result)
  statuscheck(ccall( (:cublasIsamax_v2, libcublas), cublasStatus_t, (cublasHandle_t, Cint, Ptr{Cfloat}, Cint, Ptr{Cint}), handle, n, x, incx, result))
end
function cublasIdamax_v2(handle, n, x, incx, result)
  statuscheck(ccall( (:cublasIdamax_v2, libcublas), cublasStatus_t, (cublasHandle_t, Cint, Ptr{Cdouble}, Cint, Ptr{Cint}), handle, n, x, incx, result))
end
function cublasIcamax_v2(handle, n, x, incx, result)
  statuscheck(ccall( (:cublasIcamax_v2, libcublas), cublasStatus_t, (cublasHandle_t, Cint, Ptr{cuComplex}, Cint, Ptr{Cint}), handle, n, x, incx, result))
end
function cublasIzamax_v2(handle, n, x, incx, result)
  statuscheck(ccall( (:cublasIzamax_v2, libcublas), cublasStatus_t, (cublasHandle_t, Cint, Ptr{cuDoubleComplex}, Cint, Ptr{Cint}), handle, n, x, incx, result))
end
function cublasIsamin_v2(handle, n, x, incx, result)
  statuscheck(ccall( (:cublasIsamin_v2, libcublas), cublasStatus_t, (cublasHandle_t, Cint, Ptr{Cfloat}, Cint, Ptr{Cint}), handle, n, x, incx, result))
end
function cublasIdamin_v2(handle, n, x, incx, result)
  statuscheck(ccall( (:cublasIdamin_v2, libcublas), cublasStatus_t, (cublasHandle_t, Cint, Ptr{Cdouble}, Cint, Ptr{Cint}), handle, n, x, incx, result))
end
function cublasIcamin_v2(handle, n, x, incx, result)
  statuscheck(ccall( (:cublasIcamin_v2, libcublas), cublasStatus_t, (cublasHandle_t, Cint, Ptr{cuComplex}, Cint, Ptr{Cint}), handle, n, x, incx, result))
end
function cublasIzamin_v2(handle, n, x, incx, result)
  statuscheck(ccall( (:cublasIzamin_v2, libcublas), cublasStatus_t, (cublasHandle_t, Cint, Ptr{cuDoubleComplex}, Cint, Ptr{Cint}), handle, n, x, incx, result))
end
function cublasSasum_v2(handle, n, x, incx, result)
  statuscheck(ccall( (:cublasSasum_v2, libcublas), cublasStatus_t, (cublasHandle_t, Cint, Ptr{Cfloat}, Cint, Ptr{Cfloat}), handle, n, x, incx, result))
end
function cublasDasum_v2(handle, n, x, incx, result)
  statuscheck(ccall( (:cublasDasum_v2, libcublas), cublasStatus_t, (cublasHandle_t, Cint, Ptr{Cdouble}, Cint, Ptr{Cdouble}), handle, n, x, incx, result))
end
function cublasScasum_v2(handle, n, x, incx, result)
  statuscheck(ccall( (:cublasScasum_v2, libcublas), cublasStatus_t, (cublasHandle_t, Cint, Ptr{cuComplex}, Cint, Ptr{Cfloat}), handle, n, x, incx, result))
end
function cublasDzasum_v2(handle, n, x, incx, result)
  statuscheck(ccall( (:cublasDzasum_v2, libcublas), cublasStatus_t, (cublasHandle_t, Cint, Ptr{cuDoubleComplex}, Cint, Ptr{Cdouble}), handle, n, x, incx, result))
end
function cublasSrot_v2(handle, n, x, incx, y, incy, c, s)
  statuscheck(ccall( (:cublasSrot_v2, libcublas), cublasStatus_t, (cublasHandle_t, Cint, Ptr{Cfloat}, Cint, Ptr{Cfloat}, Cint, Ptr{Cfloat}, Ptr{Cfloat}), handle, n, x, incx, y, incy, c, s))
end
function cublasDrot_v2(handle, n, x, incx, y, incy, c, s)
  statuscheck(ccall( (:cublasDrot_v2, libcublas), cublasStatus_t, (cublasHandle_t, Cint, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Ptr{Cdouble}), handle, n, x, incx, y, incy, c, s))
end
function cublasCrot_v2(handle, n, x, incx, y, incy, c, s)
  statuscheck(ccall( (:cublasCrot_v2, libcublas), cublasStatus_t, (cublasHandle_t, Cint, Ptr{cuComplex}, Cint, Ptr{cuComplex}, Cint, Ptr{Cfloat}, Ptr{cuComplex}), handle, n, x, incx, y, incy, c, s))
end
function cublasCsrot_v2(handle, n, x, incx, y, incy, c, s)
  statuscheck(ccall( (:cublasCsrot_v2, libcublas), cublasStatus_t, (cublasHandle_t, Cint, Ptr{cuComplex}, Cint, Ptr{cuComplex}, Cint, Ptr{Cfloat}, Ptr{Cfloat}), handle, n, x, incx, y, incy, c, s))
end
function cublasZrot_v2(handle, n, x, incx, y, incy, c, s)
  statuscheck(ccall( (:cublasZrot_v2, libcublas), cublasStatus_t, (cublasHandle_t, Cint, Ptr{cuDoubleComplex}, Cint, Ptr{cuDoubleComplex}, Cint, Ptr{Cdouble}, Ptr{cuDoubleComplex}), handle, n, x, incx, y, incy, c, s))
end
function cublasZdrot_v2(handle, n, x, incx, y, incy, c, s)
  statuscheck(ccall( (:cublasZdrot_v2, libcublas), cublasStatus_t, (cublasHandle_t, Cint, Ptr{cuDoubleComplex}, Cint, Ptr{cuDoubleComplex}, Cint, Ptr{Cdouble}, Ptr{Cdouble}), handle, n, x, incx, y, incy, c, s))
end
function cublasSrotg_v2(handle, a, b, c, s)
  statuscheck(ccall( (:cublasSrotg_v2, libcublas), cublasStatus_t, (cublasHandle_t, Ptr{Cfloat}, Ptr{Cfloat}, Ptr{Cfloat}, Ptr{Cfloat}), handle, a, b, c, s))
end
function cublasDrotg_v2(handle, a, b, c, s)
  statuscheck(ccall( (:cublasDrotg_v2, libcublas), cublasStatus_t, (cublasHandle_t, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), handle, a, b, c, s))
end
function cublasCrotg_v2(handle, a, b, c, s)
  statuscheck(ccall( (:cublasCrotg_v2, libcublas), cublasStatus_t, (cublasHandle_t, Ptr{cuComplex}, Ptr{cuComplex}, Ptr{Cfloat}, Ptr{cuComplex}), handle, a, b, c, s))
end
function cublasZrotg_v2(handle, a, b, c, s)
  statuscheck(ccall( (:cublasZrotg_v2, libcublas), cublasStatus_t, (cublasHandle_t, Ptr{cuDoubleComplex}, Ptr{cuDoubleComplex}, Ptr{Cdouble}, Ptr{cuDoubleComplex}), handle, a, b, c, s))
end
function cublasSrotm_v2(handle, n, x, incx, y, incy, param)
  statuscheck(ccall( (:cublasSrotm_v2, libcublas), cublasStatus_t, (cublasHandle_t, Cint, Ptr{Cfloat}, Cint, Ptr{Cfloat}, Cint, Ptr{Cfloat}), handle, n, x, incx, y, incy, param))
end
function cublasDrotm_v2(handle, n, x, incx, y, incy, param)
  statuscheck(ccall( (:cublasDrotm_v2, libcublas), cublasStatus_t, (cublasHandle_t, Cint, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint, Ptr{Cdouble}), handle, n, x, incx, y, incy, param))
end
function cublasSrotmg_v2(handle, d1, d2, x1, y1, param)
  statuscheck(ccall( (:cublasSrotmg_v2, libcublas), cublasStatus_t, (cublasHandle_t, Ptr{Cfloat}, Ptr{Cfloat}, Ptr{Cfloat}, Ptr{Cfloat}, Ptr{Cfloat}), handle, d1, d2, x1, y1, param))
end
function cublasDrotmg_v2(handle, d1, d2, x1, y1, param)
  statuscheck(ccall( (:cublasDrotmg_v2, libcublas), cublasStatus_t, (cublasHandle_t, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), handle, d1, d2, x1, y1, param))
end
function cublasSgemv_v2(handle, trans, m, n, alpha, A, lda, x, incx, beta, y, incy)
  statuscheck(ccall( (:cublasSgemv_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasOperation_t, Cint, Cint, Ptr{Cfloat}, Ptr{Cfloat}, Cint, Ptr{Cfloat}, Cint, Ptr{Cfloat}, Ptr{Cfloat}, Cint), handle, trans, m, n, alpha, A, lda, x, incx, beta, y, incy))
end
function cublasDgemv_v2(handle, trans, m, n, alpha, A, lda, x, incx, beta, y, incy)
  statuscheck(ccall( (:cublasDgemv_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasOperation_t, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cint), handle, trans, m, n, alpha, A, lda, x, incx, beta, y, incy))
end
function cublasCgemv_v2(handle, trans, m, n, alpha, A, lda, x, incx, beta, y, incy)
  statuscheck(ccall( (:cublasCgemv_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasOperation_t, Cint, Cint, Ptr{cuComplex}, Ptr{cuComplex}, Cint, Ptr{cuComplex}, Cint, Ptr{cuComplex}, Ptr{cuComplex}, Cint), handle, trans, m, n, alpha, A, lda, x, incx, beta, y, incy))
end
function cublasZgemv_v2(handle, trans, m, n, alpha, A, lda, x, incx, beta, y, incy)
  statuscheck(ccall( (:cublasZgemv_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasOperation_t, Cint, Cint, Ptr{cuDoubleComplex}, Ptr{cuDoubleComplex}, Cint, Ptr{cuDoubleComplex}, Cint, Ptr{cuDoubleComplex}, Ptr{cuDoubleComplex}, Cint), handle, trans, m, n, alpha, A, lda, x, incx, beta, y, incy))
end
function cublasSgbmv_v2(handle, trans, m, n, kl, ku, alpha, A, lda, x, incx, beta, y, incy)
  statuscheck(ccall( (:cublasSgbmv_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasOperation_t, Cint, Cint, Cint, Cint, Ptr{Cfloat}, Ptr{Cfloat}, Cint, Ptr{Cfloat}, Cint, Ptr{Cfloat}, Ptr{Cfloat}, Cint), handle, trans, m, n, kl, ku, alpha, A, lda, x, incx, beta, y, incy))
end
function cublasDgbmv_v2(handle, trans, m, n, kl, ku, alpha, A, lda, x, incx, beta, y, incy)
  statuscheck(ccall( (:cublasDgbmv_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasOperation_t, Cint, Cint, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cint), handle, trans, m, n, kl, ku, alpha, A, lda, x, incx, beta, y, incy))
end
function cublasCgbmv_v2(handle, trans, m, n, kl, ku, alpha, A, lda, x, incx, beta, y, incy)
  statuscheck(ccall( (:cublasCgbmv_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasOperation_t, Cint, Cint, Cint, Cint, Ptr{cuComplex}, Ptr{cuComplex}, Cint, Ptr{cuComplex}, Cint, Ptr{cuComplex}, Ptr{cuComplex}, Cint), handle, trans, m, n, kl, ku, alpha, A, lda, x, incx, beta, y, incy))
end
function cublasZgbmv_v2(handle, trans, m, n, kl, ku, alpha, A, lda, x, incx, beta, y, incy)
  statuscheck(ccall( (:cublasZgbmv_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasOperation_t, Cint, Cint, Cint, Cint, Ptr{cuDoubleComplex}, Ptr{cuDoubleComplex}, Cint, Ptr{cuDoubleComplex}, Cint, Ptr{cuDoubleComplex}, Ptr{cuDoubleComplex}, Cint), handle, trans, m, n, kl, ku, alpha, A, lda, x, incx, beta, y, incy))
end
function cublasStrmv_v2(handle, uplo, trans, diag, n, A, lda, x, incx)
  statuscheck(ccall( (:cublasStrmv_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasFillMode_t, cublasOperation_t, cublasDiagType_t, Cint, Ptr{Cfloat}, Cint, Ptr{Cfloat}, Cint), handle, uplo, trans, diag, n, A, lda, x, incx))
end
function cublasDtrmv_v2(handle, uplo, trans, diag, n, A, lda, x, incx)
  statuscheck(ccall( (:cublasDtrmv_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasFillMode_t, cublasOperation_t, cublasDiagType_t, Cint, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint), handle, uplo, trans, diag, n, A, lda, x, incx))
end
function cublasCtrmv_v2(handle, uplo, trans, diag, n, A, lda, x, incx)
  statuscheck(ccall( (:cublasCtrmv_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasFillMode_t, cublasOperation_t, cublasDiagType_t, Cint, Ptr{cuComplex}, Cint, Ptr{cuComplex}, Cint), handle, uplo, trans, diag, n, A, lda, x, incx))
end
function cublasZtrmv_v2(handle, uplo, trans, diag, n, A, lda, x, incx)
  statuscheck(ccall( (:cublasZtrmv_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasFillMode_t, cublasOperation_t, cublasDiagType_t, Cint, Ptr{cuDoubleComplex}, Cint, Ptr{cuDoubleComplex}, Cint), handle, uplo, trans, diag, n, A, lda, x, incx))
end
function cublasStbmv_v2(handle, uplo, trans, diag, n, k, A, lda, x, incx)
  statuscheck(ccall( (:cublasStbmv_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasFillMode_t, cublasOperation_t, cublasDiagType_t, Cint, Cint, Ptr{Cfloat}, Cint, Ptr{Cfloat}, Cint), handle, uplo, trans, diag, n, k, A, lda, x, incx))
end
function cublasDtbmv_v2(handle, uplo, trans, diag, n, k, A, lda, x, incx)
  statuscheck(ccall( (:cublasDtbmv_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasFillMode_t, cublasOperation_t, cublasDiagType_t, Cint, Cint, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint), handle, uplo, trans, diag, n, k, A, lda, x, incx))
end
function cublasCtbmv_v2(handle, uplo, trans, diag, n, k, A, lda, x, incx)
  statuscheck(ccall( (:cublasCtbmv_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasFillMode_t, cublasOperation_t, cublasDiagType_t, Cint, Cint, Ptr{cuComplex}, Cint, Ptr{cuComplex}, Cint), handle, uplo, trans, diag, n, k, A, lda, x, incx))
end
function cublasZtbmv_v2(handle, uplo, trans, diag, n, k, A, lda, x, incx)
  statuscheck(ccall( (:cublasZtbmv_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasFillMode_t, cublasOperation_t, cublasDiagType_t, Cint, Cint, Ptr{cuDoubleComplex}, Cint, Ptr{cuDoubleComplex}, Cint), handle, uplo, trans, diag, n, k, A, lda, x, incx))
end
function cublasStpmv_v2(handle, uplo, trans, diag, n, AP, x, incx)
  statuscheck(ccall( (:cublasStpmv_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasFillMode_t, cublasOperation_t, cublasDiagType_t, Cint, Ptr{Cfloat}, Ptr{Cfloat}, Cint), handle, uplo, trans, diag, n, AP, x, incx))
end
function cublasDtpmv_v2(handle, uplo, trans, diag, n, AP, x, incx)
  statuscheck(ccall( (:cublasDtpmv_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasFillMode_t, cublasOperation_t, cublasDiagType_t, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cint), handle, uplo, trans, diag, n, AP, x, incx))
end
function cublasCtpmv_v2(handle, uplo, trans, diag, n, AP, x, incx)
  statuscheck(ccall( (:cublasCtpmv_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasFillMode_t, cublasOperation_t, cublasDiagType_t, Cint, Ptr{cuComplex}, Ptr{cuComplex}, Cint), handle, uplo, trans, diag, n, AP, x, incx))
end
function cublasZtpmv_v2(handle, uplo, trans, diag, n, AP, x, incx)
  statuscheck(ccall( (:cublasZtpmv_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasFillMode_t, cublasOperation_t, cublasDiagType_t, Cint, Ptr{cuDoubleComplex}, Ptr{cuDoubleComplex}, Cint), handle, uplo, trans, diag, n, AP, x, incx))
end
function cublasStrsv_v2(handle, uplo, trans, diag, n, A, lda, x, incx)
  statuscheck(ccall( (:cublasStrsv_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasFillMode_t, cublasOperation_t, cublasDiagType_t, Cint, Ptr{Cfloat}, Cint, Ptr{Cfloat}, Cint), handle, uplo, trans, diag, n, A, lda, x, incx))
end
function cublasDtrsv_v2(handle, uplo, trans, diag, n, A, lda, x, incx)
  statuscheck(ccall( (:cublasDtrsv_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasFillMode_t, cublasOperation_t, cublasDiagType_t, Cint, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint), handle, uplo, trans, diag, n, A, lda, x, incx))
end
function cublasCtrsv_v2(handle, uplo, trans, diag, n, A, lda, x, incx)
  statuscheck(ccall( (:cublasCtrsv_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasFillMode_t, cublasOperation_t, cublasDiagType_t, Cint, Ptr{cuComplex}, Cint, Ptr{cuComplex}, Cint), handle, uplo, trans, diag, n, A, lda, x, incx))
end
function cublasZtrsv_v2(handle, uplo, trans, diag, n, A, lda, x, incx)
  statuscheck(ccall( (:cublasZtrsv_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasFillMode_t, cublasOperation_t, cublasDiagType_t, Cint, Ptr{cuDoubleComplex}, Cint, Ptr{cuDoubleComplex}, Cint), handle, uplo, trans, diag, n, A, lda, x, incx))
end
function cublasStpsv_v2(handle, uplo, trans, diag, n, AP, x, incx)
  statuscheck(ccall( (:cublasStpsv_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasFillMode_t, cublasOperation_t, cublasDiagType_t, Cint, Ptr{Cfloat}, Ptr{Cfloat}, Cint), handle, uplo, trans, diag, n, AP, x, incx))
end
function cublasDtpsv_v2(handle, uplo, trans, diag, n, AP, x, incx)
  statuscheck(ccall( (:cublasDtpsv_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasFillMode_t, cublasOperation_t, cublasDiagType_t, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cint), handle, uplo, trans, diag, n, AP, x, incx))
end
function cublasCtpsv_v2(handle, uplo, trans, diag, n, AP, x, incx)
  statuscheck(ccall( (:cublasCtpsv_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasFillMode_t, cublasOperation_t, cublasDiagType_t, Cint, Ptr{cuComplex}, Ptr{cuComplex}, Cint), handle, uplo, trans, diag, n, AP, x, incx))
end
function cublasZtpsv_v2(handle, uplo, trans, diag, n, AP, x, incx)
  statuscheck(ccall( (:cublasZtpsv_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasFillMode_t, cublasOperation_t, cublasDiagType_t, Cint, Ptr{cuDoubleComplex}, Ptr{cuDoubleComplex}, Cint), handle, uplo, trans, diag, n, AP, x, incx))
end
function cublasStbsv_v2(handle, uplo, trans, diag, n, k, A, lda, x, incx)
  statuscheck(ccall( (:cublasStbsv_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasFillMode_t, cublasOperation_t, cublasDiagType_t, Cint, Cint, Ptr{Cfloat}, Cint, Ptr{Cfloat}, Cint), handle, uplo, trans, diag, n, k, A, lda, x, incx))
end
function cublasDtbsv_v2(handle, uplo, trans, diag, n, k, A, lda, x, incx)
  statuscheck(ccall( (:cublasDtbsv_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasFillMode_t, cublasOperation_t, cublasDiagType_t, Cint, Cint, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint), handle, uplo, trans, diag, n, k, A, lda, x, incx))
end
function cublasCtbsv_v2(handle, uplo, trans, diag, n, k, A, lda, x, incx)
  statuscheck(ccall( (:cublasCtbsv_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasFillMode_t, cublasOperation_t, cublasDiagType_t, Cint, Cint, Ptr{cuComplex}, Cint, Ptr{cuComplex}, Cint), handle, uplo, trans, diag, n, k, A, lda, x, incx))
end
function cublasZtbsv_v2(handle, uplo, trans, diag, n, k, A, lda, x, incx)
  statuscheck(ccall( (:cublasZtbsv_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasFillMode_t, cublasOperation_t, cublasDiagType_t, Cint, Cint, Ptr{cuDoubleComplex}, Cint, Ptr{cuDoubleComplex}, Cint), handle, uplo, trans, diag, n, k, A, lda, x, incx))
end
function cublasSsymv_v2(handle, uplo, n, alpha, A, lda, x, incx, beta, y, incy)
  statuscheck(ccall( (:cublasSsymv_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasFillMode_t, Cint, Ptr{Cfloat}, Ptr{Cfloat}, Cint, Ptr{Cfloat}, Cint, Ptr{Cfloat}, Ptr{Cfloat}, Cint), handle, uplo, n, alpha, A, lda, x, incx, beta, y, incy))
end
function cublasDsymv_v2(handle, uplo, n, alpha, A, lda, x, incx, beta, y, incy)
  statuscheck(ccall( (:cublasDsymv_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasFillMode_t, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cint), handle, uplo, n, alpha, A, lda, x, incx, beta, y, incy))
end
function cublasCsymv_v2(handle, uplo, n, alpha, A, lda, x, incx, beta, y, incy)
  statuscheck(ccall( (:cublasCsymv_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasFillMode_t, Cint, Ptr{cuComplex}, Ptr{cuComplex}, Cint, Ptr{cuComplex}, Cint, Ptr{cuComplex}, Ptr{cuComplex}, Cint), handle, uplo, n, alpha, A, lda, x, incx, beta, y, incy))
end
function cublasZsymv_v2(handle, uplo, n, alpha, A, lda, x, incx, beta, y, incy)
  statuscheck(ccall( (:cublasZsymv_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasFillMode_t, Cint, Ptr{cuDoubleComplex}, Ptr{cuDoubleComplex}, Cint, Ptr{cuDoubleComplex}, Cint, Ptr{cuDoubleComplex}, Ptr{cuDoubleComplex}, Cint), handle, uplo, n, alpha, A, lda, x, incx, beta, y, incy))
end
function cublasChemv_v2(handle, uplo, n, alpha, A, lda, x, incx, beta, y, incy)
  statuscheck(ccall( (:cublasChemv_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasFillMode_t, Cint, Ptr{cuComplex}, Ptr{cuComplex}, Cint, Ptr{cuComplex}, Cint, Ptr{cuComplex}, Ptr{cuComplex}, Cint), handle, uplo, n, alpha, A, lda, x, incx, beta, y, incy))
end
function cublasZhemv_v2(handle, uplo, n, alpha, A, lda, x, incx, beta, y, incy)
  statuscheck(ccall( (:cublasZhemv_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasFillMode_t, Cint, Ptr{cuDoubleComplex}, Ptr{cuDoubleComplex}, Cint, Ptr{cuDoubleComplex}, Cint, Ptr{cuDoubleComplex}, Ptr{cuDoubleComplex}, Cint), handle, uplo, n, alpha, A, lda, x, incx, beta, y, incy))
end
function cublasSsbmv_v2(handle, uplo, n, k, alpha, A, lda, x, incx, beta, y, incy)
  statuscheck(ccall( (:cublasSsbmv_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasFillMode_t, Cint, Cint, Ptr{Cfloat}, Ptr{Cfloat}, Cint, Ptr{Cfloat}, Cint, Ptr{Cfloat}, Ptr{Cfloat}, Cint), handle, uplo, n, k, alpha, A, lda, x, incx, beta, y, incy))
end
function cublasDsbmv_v2(handle, uplo, n, k, alpha, A, lda, x, incx, beta, y, incy)
  statuscheck(ccall( (:cublasDsbmv_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasFillMode_t, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cint), handle, uplo, n, k, alpha, A, lda, x, incx, beta, y, incy))
end
function cublasChbmv_v2(handle, uplo, n, k, alpha, A, lda, x, incx, beta, y, incy)
  statuscheck(ccall( (:cublasChbmv_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasFillMode_t, Cint, Cint, Ptr{cuComplex}, Ptr{cuComplex}, Cint, Ptr{cuComplex}, Cint, Ptr{cuComplex}, Ptr{cuComplex}, Cint), handle, uplo, n, k, alpha, A, lda, x, incx, beta, y, incy))
end
function cublasZhbmv_v2(handle, uplo, n, k, alpha, A, lda, x, incx, beta, y, incy)
  statuscheck(ccall( (:cublasZhbmv_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasFillMode_t, Cint, Cint, Ptr{cuDoubleComplex}, Ptr{cuDoubleComplex}, Cint, Ptr{cuDoubleComplex}, Cint, Ptr{cuDoubleComplex}, Ptr{cuDoubleComplex}, Cint), handle, uplo, n, k, alpha, A, lda, x, incx, beta, y, incy))
end
function cublasSspmv_v2(handle, uplo, n, alpha, AP, x, incx, beta, y, incy)
  statuscheck(ccall( (:cublasSspmv_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasFillMode_t, Cint, Ptr{Cfloat}, Ptr{Cfloat}, Ptr{Cfloat}, Cint, Ptr{Cfloat}, Ptr{Cfloat}, Cint), handle, uplo, n, alpha, AP, x, incx, beta, y, incy))
end
function cublasDspmv_v2(handle, uplo, n, alpha, AP, x, incx, beta, y, incy)
  statuscheck(ccall( (:cublasDspmv_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasFillMode_t, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cint), handle, uplo, n, alpha, AP, x, incx, beta, y, incy))
end
function cublasChpmv_v2(handle, uplo, n, alpha, AP, x, incx, beta, y, incy)
  statuscheck(ccall( (:cublasChpmv_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasFillMode_t, Cint, Ptr{cuComplex}, Ptr{cuComplex}, Ptr{cuComplex}, Cint, Ptr{cuComplex}, Ptr{cuComplex}, Cint), handle, uplo, n, alpha, AP, x, incx, beta, y, incy))
end
function cublasZhpmv_v2(handle, uplo, n, alpha, AP, x, incx, beta, y, incy)
  statuscheck(ccall( (:cublasZhpmv_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasFillMode_t, Cint, Ptr{cuDoubleComplex}, Ptr{cuDoubleComplex}, Ptr{cuDoubleComplex}, Cint, Ptr{cuDoubleComplex}, Ptr{cuDoubleComplex}, Cint), handle, uplo, n, alpha, AP, x, incx, beta, y, incy))
end
function cublasSger_v2(handle, m, n, alpha, x, incx, y, incy, A, lda)
  statuscheck(ccall( (:cublasSger_v2, libcublas), cublasStatus_t, (cublasHandle_t, Cint, Cint, Ptr{Cfloat}, Ptr{Cfloat}, Cint, Ptr{Cfloat}, Cint, Ptr{Cfloat}, Cint), handle, m, n, alpha, x, incx, y, incy, A, lda))
end
function cublasDger_v2(handle, m, n, alpha, x, incx, y, incy, A, lda)
  statuscheck(ccall( (:cublasDger_v2, libcublas), cublasStatus_t, (cublasHandle_t, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint), handle, m, n, alpha, x, incx, y, incy, A, lda))
end
function cublasCgeru_v2(handle, m, n, alpha, x, incx, y, incy, A, lda)
  statuscheck(ccall( (:cublasCgeru_v2, libcublas), cublasStatus_t, (cublasHandle_t, Cint, Cint, Ptr{cuComplex}, Ptr{cuComplex}, Cint, Ptr{cuComplex}, Cint, Ptr{cuComplex}, Cint), handle, m, n, alpha, x, incx, y, incy, A, lda))
end
function cublasCgerc_v2(handle, m, n, alpha, x, incx, y, incy, A, lda)
  statuscheck(ccall( (:cublasCgerc_v2, libcublas), cublasStatus_t, (cublasHandle_t, Cint, Cint, Ptr{cuComplex}, Ptr{cuComplex}, Cint, Ptr{cuComplex}, Cint, Ptr{cuComplex}, Cint), handle, m, n, alpha, x, incx, y, incy, A, lda))
end
function cublasZgeru_v2(handle, m, n, alpha, x, incx, y, incy, A, lda)
  statuscheck(ccall( (:cublasZgeru_v2, libcublas), cublasStatus_t, (cublasHandle_t, Cint, Cint, Ptr{cuDoubleComplex}, Ptr{cuDoubleComplex}, Cint, Ptr{cuDoubleComplex}, Cint, Ptr{cuDoubleComplex}, Cint), handle, m, n, alpha, x, incx, y, incy, A, lda))
end
function cublasZgerc_v2(handle, m, n, alpha, x, incx, y, incy, A, lda)
  statuscheck(ccall( (:cublasZgerc_v2, libcublas), cublasStatus_t, (cublasHandle_t, Cint, Cint, Ptr{cuDoubleComplex}, Ptr{cuDoubleComplex}, Cint, Ptr{cuDoubleComplex}, Cint, Ptr{cuDoubleComplex}, Cint), handle, m, n, alpha, x, incx, y, incy, A, lda))
end
function cublasSsyr_v2(handle, uplo, n, alpha, x, incx, A, lda)
  statuscheck(ccall( (:cublasSsyr_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasFillMode_t, Cint, Ptr{Cfloat}, Ptr{Cfloat}, Cint, Ptr{Cfloat}, Cint), handle, uplo, n, alpha, x, incx, A, lda))
end
function cublasDsyr_v2(handle, uplo, n, alpha, x, incx, A, lda)
  statuscheck(ccall( (:cublasDsyr_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasFillMode_t, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint), handle, uplo, n, alpha, x, incx, A, lda))
end
function cublasCsyr_v2(handle, uplo, n, alpha, x, incx, A, lda)
  statuscheck(ccall( (:cublasCsyr_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasFillMode_t, Cint, Ptr{cuComplex}, Ptr{cuComplex}, Cint, Ptr{cuComplex}, Cint), handle, uplo, n, alpha, x, incx, A, lda))
end
function cublasZsyr_v2(handle, uplo, n, alpha, x, incx, A, lda)
  statuscheck(ccall( (:cublasZsyr_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasFillMode_t, Cint, Ptr{cuDoubleComplex}, Ptr{cuDoubleComplex}, Cint, Ptr{cuDoubleComplex}, Cint), handle, uplo, n, alpha, x, incx, A, lda))
end
function cublasCher_v2(handle, uplo, n, alpha, x, incx, A, lda)
  statuscheck(ccall( (:cublasCher_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasFillMode_t, Cint, Ptr{Cfloat}, Ptr{cuComplex}, Cint, Ptr{cuComplex}, Cint), handle, uplo, n, alpha, x, incx, A, lda))
end
function cublasZher_v2(handle, uplo, n, alpha, x, incx, A, lda)
  statuscheck(ccall( (:cublasZher_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasFillMode_t, Cint, Ptr{Cdouble}, Ptr{cuDoubleComplex}, Cint, Ptr{cuDoubleComplex}, Cint), handle, uplo, n, alpha, x, incx, A, lda))
end
function cublasSspr_v2(handle, uplo, n, alpha, x, incx, AP)
  statuscheck(ccall( (:cublasSspr_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasFillMode_t, Cint, Ptr{Cfloat}, Ptr{Cfloat}, Cint, Ptr{Cfloat}), handle, uplo, n, alpha, x, incx, AP))
end
function cublasDspr_v2(handle, uplo, n, alpha, x, incx, AP)
  statuscheck(ccall( (:cublasDspr_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasFillMode_t, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cint, Ptr{Cdouble}), handle, uplo, n, alpha, x, incx, AP))
end
function cublasChpr_v2(handle, uplo, n, alpha, x, incx, AP)
  statuscheck(ccall( (:cublasChpr_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasFillMode_t, Cint, Ptr{Cfloat}, Ptr{cuComplex}, Cint, Ptr{cuComplex}), handle, uplo, n, alpha, x, incx, AP))
end
function cublasZhpr_v2(handle, uplo, n, alpha, x, incx, AP)
  statuscheck(ccall( (:cublasZhpr_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasFillMode_t, Cint, Ptr{Cdouble}, Ptr{cuDoubleComplex}, Cint, Ptr{cuDoubleComplex}), handle, uplo, n, alpha, x, incx, AP))
end
function cublasSsyr2_v2(handle, uplo, n, alpha, x, incx, y, incy, A, lda)
  statuscheck(ccall( (:cublasSsyr2_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasFillMode_t, Cint, Ptr{Cfloat}, Ptr{Cfloat}, Cint, Ptr{Cfloat}, Cint, Ptr{Cfloat}, Cint), handle, uplo, n, alpha, x, incx, y, incy, A, lda))
end
function cublasDsyr2_v2(handle, uplo, n, alpha, x, incx, y, incy, A, lda)
  statuscheck(ccall( (:cublasDsyr2_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasFillMode_t, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint), handle, uplo, n, alpha, x, incx, y, incy, A, lda))
end
function cublasCsyr2_v2(handle, uplo, n, alpha, x, incx, y, incy, A, lda)
  statuscheck(ccall( (:cublasCsyr2_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasFillMode_t, Cint, Ptr{cuComplex}, Ptr{cuComplex}, Cint, Ptr{cuComplex}, Cint, Ptr{cuComplex}, Cint), handle, uplo, n, alpha, x, incx, y, incy, A, lda))
end
function cublasZsyr2_v2(handle, uplo, n, alpha, x, incx, y, incy, A, lda)
  statuscheck(ccall( (:cublasZsyr2_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasFillMode_t, Cint, Ptr{cuDoubleComplex}, Ptr{cuDoubleComplex}, Cint, Ptr{cuDoubleComplex}, Cint, Ptr{cuDoubleComplex}, Cint), handle, uplo, n, alpha, x, incx, y, incy, A, lda))
end
function cublasCher2_v2(handle, uplo, n, alpha, x, incx, y, incy, A, lda)
  statuscheck(ccall( (:cublasCher2_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasFillMode_t, Cint, Ptr{cuComplex}, Ptr{cuComplex}, Cint, Ptr{cuComplex}, Cint, Ptr{cuComplex}, Cint), handle, uplo, n, alpha, x, incx, y, incy, A, lda))
end
function cublasZher2_v2(handle, uplo, n, alpha, x, incx, y, incy, A, lda)
  statuscheck(ccall( (:cublasZher2_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasFillMode_t, Cint, Ptr{cuDoubleComplex}, Ptr{cuDoubleComplex}, Cint, Ptr{cuDoubleComplex}, Cint, Ptr{cuDoubleComplex}, Cint), handle, uplo, n, alpha, x, incx, y, incy, A, lda))
end
function cublasSspr2_v2(handle, uplo, n, alpha, x, incx, y, incy, AP)
  statuscheck(ccall( (:cublasSspr2_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasFillMode_t, Cint, Ptr{Cfloat}, Ptr{Cfloat}, Cint, Ptr{Cfloat}, Cint, Ptr{Cfloat}), handle, uplo, n, alpha, x, incx, y, incy, AP))
end
function cublasDspr2_v2(handle, uplo, n, alpha, x, incx, y, incy, AP)
  statuscheck(ccall( (:cublasDspr2_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasFillMode_t, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint, Ptr{Cdouble}), handle, uplo, n, alpha, x, incx, y, incy, AP))
end
function cublasChpr2_v2(handle, uplo, n, alpha, x, incx, y, incy, AP)
  statuscheck(ccall( (:cublasChpr2_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasFillMode_t, Cint, Ptr{cuComplex}, Ptr{cuComplex}, Cint, Ptr{cuComplex}, Cint, Ptr{cuComplex}), handle, uplo, n, alpha, x, incx, y, incy, AP))
end
function cublasZhpr2_v2(handle, uplo, n, alpha, x, incx, y, incy, AP)
  statuscheck(ccall( (:cublasZhpr2_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasFillMode_t, Cint, Ptr{cuDoubleComplex}, Ptr{cuDoubleComplex}, Cint, Ptr{cuDoubleComplex}, Cint, Ptr{cuDoubleComplex}), handle, uplo, n, alpha, x, incx, y, incy, AP))
end
function cublasSgemm_v2(handle, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
  statuscheck(ccall( (:cublasSgemm_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasOperation_t, cublasOperation_t, Cint, Cint, Cint, Ptr{Cfloat}, Ptr{Cfloat}, Cint, Ptr{Cfloat}, Cint, Ptr{Cfloat}, Ptr{Cfloat}, Cint), handle, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc))
end
function cublasDgemm_v2(handle, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
  statuscheck(ccall( (:cublasDgemm_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasOperation_t, cublasOperation_t, Cint, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cint), handle, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc))
end
function cublasCgemm_v2(handle, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
  statuscheck(ccall( (:cublasCgemm_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasOperation_t, cublasOperation_t, Cint, Cint, Cint, Ptr{cuComplex}, Ptr{cuComplex}, Cint, Ptr{cuComplex}, Cint, Ptr{cuComplex}, Ptr{cuComplex}, Cint), handle, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc))
end
function cublasZgemm_v2(handle, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
  statuscheck(ccall( (:cublasZgemm_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasOperation_t, cublasOperation_t, Cint, Cint, Cint, Ptr{cuDoubleComplex}, Ptr{cuDoubleComplex}, Cint, Ptr{cuDoubleComplex}, Cint, Ptr{cuDoubleComplex}, Ptr{cuDoubleComplex}, Cint), handle, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc))
end
function cublasSsyrk_v2(handle, uplo, trans, n, k, alpha, A, lda, beta, C, ldc)
  statuscheck(ccall( (:cublasSsyrk_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasFillMode_t, cublasOperation_t, Cint, Cint, Ptr{Cfloat}, Ptr{Cfloat}, Cint, Ptr{Cfloat}, Ptr{Cfloat}, Cint), handle, uplo, trans, n, k, alpha, A, lda, beta, C, ldc))
end
function cublasDsyrk_v2(handle, uplo, trans, n, k, alpha, A, lda, beta, C, ldc)
  statuscheck(ccall( (:cublasDsyrk_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasFillMode_t, cublasOperation_t, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cint), handle, uplo, trans, n, k, alpha, A, lda, beta, C, ldc))
end
function cublasCsyrk_v2(handle, uplo, trans, n, k, alpha, A, lda, beta, C, ldc)
  statuscheck(ccall( (:cublasCsyrk_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasFillMode_t, cublasOperation_t, Cint, Cint, Ptr{cuComplex}, Ptr{cuComplex}, Cint, Ptr{cuComplex}, Ptr{cuComplex}, Cint), handle, uplo, trans, n, k, alpha, A, lda, beta, C, ldc))
end
function cublasZsyrk_v2(handle, uplo, trans, n, k, alpha, A, lda, beta, C, ldc)
  statuscheck(ccall( (:cublasZsyrk_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasFillMode_t, cublasOperation_t, Cint, Cint, Ptr{cuDoubleComplex}, Ptr{cuDoubleComplex}, Cint, Ptr{cuDoubleComplex}, Ptr{cuDoubleComplex}, Cint), handle, uplo, trans, n, k, alpha, A, lda, beta, C, ldc))
end
function cublasCherk_v2(handle, uplo, trans, n, k, alpha, A, lda, beta, C, ldc)
  statuscheck(ccall( (:cublasCherk_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasFillMode_t, cublasOperation_t, Cint, Cint, Ptr{Cfloat}, Ptr{cuComplex}, Cint, Ptr{Cfloat}, Ptr{cuComplex}, Cint), handle, uplo, trans, n, k, alpha, A, lda, beta, C, ldc))
end
function cublasZherk_v2(handle, uplo, trans, n, k, alpha, A, lda, beta, C, ldc)
  statuscheck(ccall( (:cublasZherk_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasFillMode_t, cublasOperation_t, Cint, Cint, Ptr{Cdouble}, Ptr{cuDoubleComplex}, Cint, Ptr{Cdouble}, Ptr{cuDoubleComplex}, Cint), handle, uplo, trans, n, k, alpha, A, lda, beta, C, ldc))
end
function cublasSsyr2k_v2(handle, uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
  statuscheck(ccall( (:cublasSsyr2k_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasFillMode_t, cublasOperation_t, Cint, Cint, Ptr{Cfloat}, Ptr{Cfloat}, Cint, Ptr{Cfloat}, Cint, Ptr{Cfloat}, Ptr{Cfloat}, Cint), handle, uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc))
end
function cublasDsyr2k_v2(handle, uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
  statuscheck(ccall( (:cublasDsyr2k_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasFillMode_t, cublasOperation_t, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cint), handle, uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc))
end
function cublasCsyr2k_v2(handle, uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
  statuscheck(ccall( (:cublasCsyr2k_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasFillMode_t, cublasOperation_t, Cint, Cint, Ptr{cuComplex}, Ptr{cuComplex}, Cint, Ptr{cuComplex}, Cint, Ptr{cuComplex}, Ptr{cuComplex}, Cint), handle, uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc))
end
function cublasZsyr2k_v2(handle, uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
  statuscheck(ccall( (:cublasZsyr2k_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasFillMode_t, cublasOperation_t, Cint, Cint, Ptr{cuDoubleComplex}, Ptr{cuDoubleComplex}, Cint, Ptr{cuDoubleComplex}, Cint, Ptr{cuDoubleComplex}, Ptr{cuDoubleComplex}, Cint), handle, uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc))
end
function cublasCher2k_v2(handle, uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
  statuscheck(ccall( (:cublasCher2k_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasFillMode_t, cublasOperation_t, Cint, Cint, Ptr{cuComplex}, Ptr{cuComplex}, Cint, Ptr{cuComplex}, Cint, Ptr{Cfloat}, Ptr{cuComplex}, Cint), handle, uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc))
end
function cublasZher2k_v2(handle, uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
  statuscheck(ccall( (:cublasZher2k_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasFillMode_t, cublasOperation_t, Cint, Cint, Ptr{cuDoubleComplex}, Ptr{cuDoubleComplex}, Cint, Ptr{cuDoubleComplex}, Cint, Ptr{Cdouble}, Ptr{cuDoubleComplex}, Cint), handle, uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc))
end
function cublasSsyrkx(handle, uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
  statuscheck(ccall( (:cublasSsyrkx, libcublas), cublasStatus_t, (cublasHandle_t, cublasFillMode_t, cublasOperation_t, Cint, Cint, Ptr{Cfloat}, Ptr{Cfloat}, Cint, Ptr{Cfloat}, Cint, Ptr{Cfloat}, Ptr{Cfloat}, Cint), handle, uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc))
end
function cublasDsyrkx(handle, uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
  statuscheck(ccall( (:cublasDsyrkx, libcublas), cublasStatus_t, (cublasHandle_t, cublasFillMode_t, cublasOperation_t, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cint), handle, uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc))
end
function cublasCsyrkx(handle, uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
  statuscheck(ccall( (:cublasCsyrkx, libcublas), cublasStatus_t, (cublasHandle_t, cublasFillMode_t, cublasOperation_t, Cint, Cint, Ptr{cuComplex}, Ptr{cuComplex}, Cint, Ptr{cuComplex}, Cint, Ptr{cuComplex}, Ptr{cuComplex}, Cint), handle, uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc))
end
function cublasZsyrkx(handle, uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
  statuscheck(ccall( (:cublasZsyrkx, libcublas), cublasStatus_t, (cublasHandle_t, cublasFillMode_t, cublasOperation_t, Cint, Cint, Ptr{cuDoubleComplex}, Ptr{cuDoubleComplex}, Cint, Ptr{cuDoubleComplex}, Cint, Ptr{cuDoubleComplex}, Ptr{cuDoubleComplex}, Cint), handle, uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc))
end
function cublasCherkx(handle, uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
  statuscheck(ccall( (:cublasCherkx, libcublas), cublasStatus_t, (cublasHandle_t, cublasFillMode_t, cublasOperation_t, Cint, Cint, Ptr{cuComplex}, Ptr{cuComplex}, Cint, Ptr{cuComplex}, Cint, Ptr{Cfloat}, Ptr{cuComplex}, Cint), handle, uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc))
end
function cublasZherkx(handle, uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
  statuscheck(ccall( (:cublasZherkx, libcublas), cublasStatus_t, (cublasHandle_t, cublasFillMode_t, cublasOperation_t, Cint, Cint, Ptr{cuDoubleComplex}, Ptr{cuDoubleComplex}, Cint, Ptr{cuDoubleComplex}, Cint, Ptr{Cdouble}, Ptr{cuDoubleComplex}, Cint), handle, uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc))
end
function cublasSsymm_v2(handle, side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc)
  statuscheck(ccall( (:cublasSsymm_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasSideMode_t, cublasFillMode_t, Cint, Cint, Ptr{Cfloat}, Ptr{Cfloat}, Cint, Ptr{Cfloat}, Cint, Ptr{Cfloat}, Ptr{Cfloat}, Cint), handle, side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc))
end
function cublasDsymm_v2(handle, side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc)
  statuscheck(ccall( (:cublasDsymm_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasSideMode_t, cublasFillMode_t, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cint), handle, side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc))
end
function cublasCsymm_v2(handle, side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc)
  statuscheck(ccall( (:cublasCsymm_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasSideMode_t, cublasFillMode_t, Cint, Cint, Ptr{cuComplex}, Ptr{cuComplex}, Cint, Ptr{cuComplex}, Cint, Ptr{cuComplex}, Ptr{cuComplex}, Cint), handle, side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc))
end
function cublasZsymm_v2(handle, side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc)
  statuscheck(ccall( (:cublasZsymm_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasSideMode_t, cublasFillMode_t, Cint, Cint, Ptr{cuDoubleComplex}, Ptr{cuDoubleComplex}, Cint, Ptr{cuDoubleComplex}, Cint, Ptr{cuDoubleComplex}, Ptr{cuDoubleComplex}, Cint), handle, side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc))
end
function cublasChemm_v2(handle, side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc)
  statuscheck(ccall( (:cublasChemm_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasSideMode_t, cublasFillMode_t, Cint, Cint, Ptr{cuComplex}, Ptr{cuComplex}, Cint, Ptr{cuComplex}, Cint, Ptr{cuComplex}, Ptr{cuComplex}, Cint), handle, side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc))
end
function cublasZhemm_v2(handle, side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc)
  statuscheck(ccall( (:cublasZhemm_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasSideMode_t, cublasFillMode_t, Cint, Cint, Ptr{cuDoubleComplex}, Ptr{cuDoubleComplex}, Cint, Ptr{cuDoubleComplex}, Cint, Ptr{cuDoubleComplex}, Ptr{cuDoubleComplex}, Cint), handle, side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc))
end
function cublasStrsm_v2(handle, side, uplo, trans, diag, m, n, alpha, A, lda, B, ldb)
  statuscheck(ccall( (:cublasStrsm_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasSideMode_t, cublasFillMode_t, cublasOperation_t, cublasDiagType_t, Cint, Cint, Ptr{Cfloat}, Ptr{Cfloat}, Cint, Ptr{Cfloat}, Cint), handle, side, uplo, trans, diag, m, n, alpha, A, lda, B, ldb))
end
function cublasDtrsm_v2(handle, side, uplo, trans, diag, m, n, alpha, A, lda, B, ldb)
  statuscheck(ccall( (:cublasDtrsm_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasSideMode_t, cublasFillMode_t, cublasOperation_t, cublasDiagType_t, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint), handle, side, uplo, trans, diag, m, n, alpha, A, lda, B, ldb))
end
function cublasCtrsm_v2(handle, side, uplo, trans, diag, m, n, alpha, A, lda, B, ldb)
  statuscheck(ccall( (:cublasCtrsm_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasSideMode_t, cublasFillMode_t, cublasOperation_t, cublasDiagType_t, Cint, Cint, Ptr{cuComplex}, Ptr{cuComplex}, Cint, Ptr{cuComplex}, Cint), handle, side, uplo, trans, diag, m, n, alpha, A, lda, B, ldb))
end
function cublasZtrsm_v2(handle, side, uplo, trans, diag, m, n, alpha, A, lda, B, ldb)
  statuscheck(ccall( (:cublasZtrsm_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasSideMode_t, cublasFillMode_t, cublasOperation_t, cublasDiagType_t, Cint, Cint, Ptr{cuDoubleComplex}, Ptr{cuDoubleComplex}, Cint, Ptr{cuDoubleComplex}, Cint), handle, side, uplo, trans, diag, m, n, alpha, A, lda, B, ldb))
end
function cublasStrmm_v2(handle, side, uplo, trans, diag, m, n, alpha, A, lda, B, ldb, C, ldc)
  statuscheck(ccall( (:cublasStrmm_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasSideMode_t, cublasFillMode_t, cublasOperation_t, cublasDiagType_t, Cint, Cint, Ptr{Cfloat}, Ptr{Cfloat}, Cint, Ptr{Cfloat}, Cint, Ptr{Cfloat}, Cint), handle, side, uplo, trans, diag, m, n, alpha, A, lda, B, ldb, C, ldc))
end
function cublasDtrmm_v2(handle, side, uplo, trans, diag, m, n, alpha, A, lda, B, ldb, C, ldc)
  statuscheck(ccall( (:cublasDtrmm_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasSideMode_t, cublasFillMode_t, cublasOperation_t, cublasDiagType_t, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint), handle, side, uplo, trans, diag, m, n, alpha, A, lda, B, ldb, C, ldc))
end
function cublasCtrmm_v2(handle, side, uplo, trans, diag, m, n, alpha, A, lda, B, ldb, C, ldc)
  statuscheck(ccall( (:cublasCtrmm_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasSideMode_t, cublasFillMode_t, cublasOperation_t, cublasDiagType_t, Cint, Cint, Ptr{cuComplex}, Ptr{cuComplex}, Cint, Ptr{cuComplex}, Cint, Ptr{cuComplex}, Cint), handle, side, uplo, trans, diag, m, n, alpha, A, lda, B, ldb, C, ldc))
end
function cublasZtrmm_v2(handle, side, uplo, trans, diag, m, n, alpha, A, lda, B, ldb, C, ldc)
  statuscheck(ccall( (:cublasZtrmm_v2, libcublas), cublasStatus_t, (cublasHandle_t, cublasSideMode_t, cublasFillMode_t, cublasOperation_t, cublasDiagType_t, Cint, Cint, Ptr{cuDoubleComplex}, Ptr{cuDoubleComplex}, Cint, Ptr{cuDoubleComplex}, Cint, Ptr{cuDoubleComplex}, Cint), handle, side, uplo, trans, diag, m, n, alpha, A, lda, B, ldb, C, ldc))
end
function cublasSgemmBatched(handle, transa, transb, m, n, k, alpha, Aarray, lda, Barray, ldb, beta, Carray, ldc, batchCount)
  statuscheck(ccall( (:cublasSgemmBatched, libcublas), cublasStatus_t, (cublasHandle_t, cublasOperation_t, cublasOperation_t, Cint, Cint, Cint, Ptr{Cfloat}, Ptr{Ptr{Cfloat}}, Cint, Ptr{Ptr{Cfloat}}, Cint, Ptr{Cfloat}, Ptr{Ptr{Cfloat}}, Cint, Cint), handle, transa, transb, m, n, k, alpha, Aarray, lda, Barray, ldb, beta, Carray, ldc, batchCount))
end
function cublasDgemmBatched(handle, transa, transb, m, n, k, alpha, Aarray, lda, Barray, ldb, beta, Carray, ldc, batchCount)
  statuscheck(ccall( (:cublasDgemmBatched, libcublas), cublasStatus_t, (cublasHandle_t, cublasOperation_t, cublasOperation_t, Cint, Cint, Cint, Ptr{Cdouble}, Ptr{Ptr{Cdouble}}, Cint, Ptr{Ptr{Cdouble}}, Cint, Ptr{Cdouble}, Ptr{Ptr{Cdouble}}, Cint, Cint), handle, transa, transb, m, n, k, alpha, Aarray, lda, Barray, ldb, beta, Carray, ldc, batchCount))
end
function cublasCgemmBatched(handle, transa, transb, m, n, k, alpha, Aarray, lda, Barray, ldb, beta, Carray, ldc, batchCount)
  statuscheck(ccall( (:cublasCgemmBatched, libcublas), cublasStatus_t, (cublasHandle_t, cublasOperation_t, cublasOperation_t, Cint, Cint, Cint, Ptr{cuComplex}, Ptr{Ptr{cuComplex}}, Cint, Ptr{Ptr{cuComplex}}, Cint, Ptr{cuComplex}, Ptr{Ptr{cuComplex}}, Cint, Cint), handle, transa, transb, m, n, k, alpha, Aarray, lda, Barray, ldb, beta, Carray, ldc, batchCount))
end
function cublasZgemmBatched(handle, transa, transb, m, n, k, alpha, Aarray, lda, Barray, ldb, beta, Carray, ldc, batchCount)
  statuscheck(ccall( (:cublasZgemmBatched, libcublas), cublasStatus_t, (cublasHandle_t, cublasOperation_t, cublasOperation_t, Cint, Cint, Cint, Ptr{cuDoubleComplex}, Ptr{Ptr{cuDoubleComplex}}, Cint, Ptr{Ptr{cuDoubleComplex}}, Cint, Ptr{cuDoubleComplex}, Ptr{Ptr{cuDoubleComplex}}, Cint, Cint), handle, transa, transb, m, n, k, alpha, Aarray, lda, Barray, ldb, beta, Carray, ldc, batchCount))
end
function cublasSgeam(handle, transa, transb, m, n, alpha, A, lda, beta, B, ldb, C, ldc)
  statuscheck(ccall( (:cublasSgeam, libcublas), cublasStatus_t, (cublasHandle_t, cublasOperation_t, cublasOperation_t, Cint, Cint, Ptr{Cfloat}, Ptr{Cfloat}, Cint, Ptr{Cfloat}, Ptr{Cfloat}, Cint, Ptr{Cfloat}, Cint), handle, transa, transb, m, n, alpha, A, lda, beta, B, ldb, C, ldc))
end
function cublasDgeam(handle, transa, transb, m, n, alpha, A, lda, beta, B, ldb, C, ldc)
  statuscheck(ccall( (:cublasDgeam, libcublas), cublasStatus_t, (cublasHandle_t, cublasOperation_t, cublasOperation_t, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint), handle, transa, transb, m, n, alpha, A, lda, beta, B, ldb, C, ldc))
end
function cublasCgeam(handle, transa, transb, m, n, alpha, A, lda, beta, B, ldb, C, ldc)
  statuscheck(ccall( (:cublasCgeam, libcublas), cublasStatus_t, (cublasHandle_t, cublasOperation_t, cublasOperation_t, Cint, Cint, Ptr{cuComplex}, Ptr{cuComplex}, Cint, Ptr{cuComplex}, Ptr{cuComplex}, Cint, Ptr{cuComplex}, Cint), handle, transa, transb, m, n, alpha, A, lda, beta, B, ldb, C, ldc))
end
function cublasZgeam(handle, transa, transb, m, n, alpha, A, lda, beta, B, ldb, C, ldc)
  statuscheck(ccall( (:cublasZgeam, libcublas), cublasStatus_t, (cublasHandle_t, cublasOperation_t, cublasOperation_t, Cint, Cint, Ptr{cuDoubleComplex}, Ptr{cuDoubleComplex}, Cint, Ptr{cuDoubleComplex}, Ptr{cuDoubleComplex}, Cint, Ptr{cuDoubleComplex}, Cint), handle, transa, transb, m, n, alpha, A, lda, beta, B, ldb, C, ldc))
end
function cublasSgetrfBatched(handle, n, A, lda, P, info, batchSize)
  statuscheck(ccall( (:cublasSgetrfBatched, libcublas), cublasStatus_t, (cublasHandle_t, Cint, Ptr{Ptr{Cfloat}}, Cint, Ptr{Cint}, Ptr{Cint}, Cint), handle, n, A, lda, P, info, batchSize))
end
function cublasDgetrfBatched(handle, n, A, lda, P, info, batchSize)
  statuscheck(ccall( (:cublasDgetrfBatched, libcublas), cublasStatus_t, (cublasHandle_t, Cint, Ptr{Ptr{Cdouble}}, Cint, Ptr{Cint}, Ptr{Cint}, Cint), handle, n, A, lda, P, info, batchSize))
end
function cublasCgetrfBatched(handle, n, A, lda, P, info, batchSize)
  statuscheck(ccall( (:cublasCgetrfBatched, libcublas), cublasStatus_t, (cublasHandle_t, Cint, Ptr{Ptr{cuComplex}}, Cint, Ptr{Cint}, Ptr{Cint}, Cint), handle, n, A, lda, P, info, batchSize))
end
function cublasZgetrfBatched(handle, n, A, lda, P, info, batchSize)
  statuscheck(ccall( (:cublasZgetrfBatched, libcublas), cublasStatus_t, (cublasHandle_t, Cint, Ptr{Ptr{cuDoubleComplex}}, Cint, Ptr{Cint}, Ptr{Cint}, Cint), handle, n, A, lda, P, info, batchSize))
end
function cublasSgetriBatched(handle, n, A, lda, P, C, ldc, info, batchSize)
  statuscheck(ccall( (:cublasSgetriBatched, libcublas), cublasStatus_t, (cublasHandle_t, Cint, Ptr{Ptr{Cfloat}}, Cint, Ptr{Cint}, Ptr{Ptr{Cfloat}}, Cint, Ptr{Cint}, Cint), handle, n, A, lda, P, C, ldc, info, batchSize))
end
function cublasDgetriBatched(handle, n, A, lda, P, C, ldc, info, batchSize)
  statuscheck(ccall( (:cublasDgetriBatched, libcublas), cublasStatus_t, (cublasHandle_t, Cint, Ptr{Ptr{Cdouble}}, Cint, Ptr{Cint}, Ptr{Ptr{Cdouble}}, Cint, Ptr{Cint}, Cint), handle, n, A, lda, P, C, ldc, info, batchSize))
end
function cublasCgetriBatched(handle, n, A, lda, P, C, ldc, info, batchSize)
  statuscheck(ccall( (:cublasCgetriBatched, libcublas), cublasStatus_t, (cublasHandle_t, Cint, Ptr{Ptr{cuComplex}}, Cint, Ptr{Cint}, Ptr{Ptr{cuComplex}}, Cint, Ptr{Cint}, Cint), handle, n, A, lda, P, C, ldc, info, batchSize))
end
function cublasZgetriBatched(handle, n, A, lda, P, C, ldc, info, batchSize)
  statuscheck(ccall( (:cublasZgetriBatched, libcublas), cublasStatus_t, (cublasHandle_t, Cint, Ptr{Ptr{cuDoubleComplex}}, Cint, Ptr{Cint}, Ptr{Ptr{cuDoubleComplex}}, Cint, Ptr{Cint}, Cint), handle, n, A, lda, P, C, ldc, info, batchSize))
end
function cublasStrsmBatched(handle, side, uplo, trans, diag, m, n, alpha, A, lda, B, ldb, batchCount)
  statuscheck(ccall( (:cublasStrsmBatched, libcublas), cublasStatus_t, (cublasHandle_t, cublasSideMode_t, cublasFillMode_t, cublasOperation_t, cublasDiagType_t, Cint, Cint, Ptr{Cfloat}, Ptr{Ptr{Cfloat}}, Cint, Ptr{Ptr{Cfloat}}, Cint, Cint), handle, side, uplo, trans, diag, m, n, alpha, A, lda, B, ldb, batchCount))
end
function cublasDtrsmBatched(handle, side, uplo, trans, diag, m, n, alpha, A, lda, B, ldb, batchCount)
  statuscheck(ccall( (:cublasDtrsmBatched, libcublas), cublasStatus_t, (cublasHandle_t, cublasSideMode_t, cublasFillMode_t, cublasOperation_t, cublasDiagType_t, Cint, Cint, Ptr{Cdouble}, Ptr{Ptr{Cdouble}}, Cint, Ptr{Ptr{Cdouble}}, Cint, Cint), handle, side, uplo, trans, diag, m, n, alpha, A, lda, B, ldb, batchCount))
end
function cublasCtrsmBatched(handle, side, uplo, trans, diag, m, n, alpha, A, lda, B, ldb, batchCount)
  statuscheck(ccall( (:cublasCtrsmBatched, libcublas), cublasStatus_t, (cublasHandle_t, cublasSideMode_t, cublasFillMode_t, cublasOperation_t, cublasDiagType_t, Cint, Cint, Ptr{cuComplex}, Ptr{Ptr{cuComplex}}, Cint, Ptr{Ptr{cuComplex}}, Cint, Cint), handle, side, uplo, trans, diag, m, n, alpha, A, lda, B, ldb, batchCount))
end
function cublasZtrsmBatched(handle, side, uplo, trans, diag, m, n, alpha, A, lda, B, ldb, batchCount)
  statuscheck(ccall( (:cublasZtrsmBatched, libcublas), cublasStatus_t, (cublasHandle_t, cublasSideMode_t, cublasFillMode_t, cublasOperation_t, cublasDiagType_t, Cint, Cint, Ptr{cuDoubleComplex}, Ptr{Ptr{cuDoubleComplex}}, Cint, Ptr{Ptr{cuDoubleComplex}}, Cint, Cint), handle, side, uplo, trans, diag, m, n, alpha, A, lda, B, ldb, batchCount))
end
function cublasSmatinvBatched(handle, n, A, lda, Ainv, lda_inv, info, batchSize)
  statuscheck(ccall( (:cublasSmatinvBatched, libcublas), cublasStatus_t, (cublasHandle_t, Cint, Ptr{Ptr{Cfloat}}, Cint, Ptr{Ptr{Cfloat}}, Cint, Ptr{Cint}, Cint), handle, n, A, lda, Ainv, lda_inv, info, batchSize))
end
function cublasDmatinvBatched(handle, n, A, lda, Ainv, lda_inv, info, batchSize)
  statuscheck(ccall( (:cublasDmatinvBatched, libcublas), cublasStatus_t, (cublasHandle_t, Cint, Ptr{Ptr{Cdouble}}, Cint, Ptr{Ptr{Cdouble}}, Cint, Ptr{Cint}, Cint), handle, n, A, lda, Ainv, lda_inv, info, batchSize))
end
function cublasCmatinvBatched(handle, n, A, lda, Ainv, lda_inv, info, batchSize)
  statuscheck(ccall( (:cublasCmatinvBatched, libcublas), cublasStatus_t, (cublasHandle_t, Cint, Ptr{Ptr{cuComplex}}, Cint, Ptr{Ptr{cuComplex}}, Cint, Ptr{Cint}, Cint), handle, n, A, lda, Ainv, lda_inv, info, batchSize))
end
function cublasZmatinvBatched(handle, n, A, lda, Ainv, lda_inv, info, batchSize)
  statuscheck(ccall( (:cublasZmatinvBatched, libcublas), cublasStatus_t, (cublasHandle_t, Cint, Ptr{Ptr{cuDoubleComplex}}, Cint, Ptr{Ptr{cuDoubleComplex}}, Cint, Ptr{Cint}, Cint), handle, n, A, lda, Ainv, lda_inv, info, batchSize))
end
function cublasSgeqrfBatched(handle, m, n, Aarray, lda, TauArray, info, batchSize)
  statuscheck(ccall( (:cublasSgeqrfBatched, libcublas), cublasStatus_t, (cublasHandle_t, Cint, Cint, Ptr{Ptr{Cfloat}}, Cint, Ptr{Ptr{Cfloat}}, Ptr{Cint}, Cint), handle, m, n, Aarray, lda, TauArray, info, batchSize))
end
function cublasDgeqrfBatched(handle, m, n, Aarray, lda, TauArray, info, batchSize)
  statuscheck(ccall( (:cublasDgeqrfBatched, libcublas), cublasStatus_t, (cublasHandle_t, Cint, Cint, Ptr{Ptr{Cdouble}}, Cint, Ptr{Ptr{Cdouble}}, Ptr{Cint}, Cint), handle, m, n, Aarray, lda, TauArray, info, batchSize))
end
function cublasCgeqrfBatched(handle, m, n, Aarray, lda, TauArray, info, batchSize)
  statuscheck(ccall( (:cublasCgeqrfBatched, libcublas), cublasStatus_t, (cublasHandle_t, Cint, Cint, Ptr{Ptr{cuComplex}}, Cint, Ptr{Ptr{cuComplex}}, Ptr{Cint}, Cint), handle, m, n, Aarray, lda, TauArray, info, batchSize))
end
function cublasZgeqrfBatched(handle, m, n, Aarray, lda, TauArray, info, batchSize)
  statuscheck(ccall( (:cublasZgeqrfBatched, libcublas), cublasStatus_t, (cublasHandle_t, Cint, Cint, Ptr{Ptr{cuDoubleComplex}}, Cint, Ptr{Ptr{cuDoubleComplex}}, Ptr{Cint}, Cint), handle, m, n, Aarray, lda, TauArray, info, batchSize))
end
function cublasSgelsBatched(handle, trans, m, n, nrhs, Aarray, lda, Carray, ldc, info, devInfoArray, batchSize)
  statuscheck(ccall( (:cublasSgelsBatched, libcublas), cublasStatus_t, (cublasHandle_t, cublasOperation_t, Cint, Cint, Cint, Ptr{Ptr{Cfloat}}, Cint, Ptr{Ptr{Cfloat}}, Cint, Ptr{Cint}, Ptr{Cint}, Cint), handle, trans, m, n, nrhs, Aarray, lda, Carray, ldc, info, devInfoArray, batchSize))
end
function cublasDgelsBatched(handle, trans, m, n, nrhs, Aarray, lda, Carray, ldc, info, devInfoArray, batchSize)
  statuscheck(ccall( (:cublasDgelsBatched, libcublas), cublasStatus_t, (cublasHandle_t, cublasOperation_t, Cint, Cint, Cint, Ptr{Ptr{Cdouble}}, Cint, Ptr{Ptr{Cdouble}}, Cint, Ptr{Cint}, Ptr{Cint}, Cint), handle, trans, m, n, nrhs, Aarray, lda, Carray, ldc, info, devInfoArray, batchSize))
end
function cublasCgelsBatched(handle, trans, m, n, nrhs, Aarray, lda, Carray, ldc, info, devInfoArray, batchSize)
  statuscheck(ccall( (:cublasCgelsBatched, libcublas), cublasStatus_t, (cublasHandle_t, cublasOperation_t, Cint, Cint, Cint, Ptr{Ptr{cuComplex}}, Cint, Ptr{Ptr{cuComplex}}, Cint, Ptr{Cint}, Ptr{Cint}, Cint), handle, trans, m, n, nrhs, Aarray, lda, Carray, ldc, info, devInfoArray, batchSize))
end
function cublasZgelsBatched(handle, trans, m, n, nrhs, Aarray, lda, Carray, ldc, info, devInfoArray, batchSize)
  statuscheck(ccall( (:cublasZgelsBatched, libcublas), cublasStatus_t, (cublasHandle_t, cublasOperation_t, Cint, Cint, Cint, Ptr{Ptr{cuDoubleComplex}}, Cint, Ptr{Ptr{cuDoubleComplex}}, Cint, Ptr{Cint}, Ptr{Cint}, Cint), handle, trans, m, n, nrhs, Aarray, lda, Carray, ldc, info, devInfoArray, batchSize))
end
function cublasSdgmm(handle, mode, m, n, A, lda, x, incx, C, ldc)
  statuscheck(ccall( (:cublasSdgmm, libcublas), cublasStatus_t, (cublasHandle_t, cublasSideMode_t, Cint, Cint, Ptr{Cfloat}, Cint, Ptr{Cfloat}, Cint, Ptr{Cfloat}, Cint), handle, mode, m, n, A, lda, x, incx, C, ldc))
end
function cublasDdgmm(handle, mode, m, n, A, lda, x, incx, C, ldc)
  statuscheck(ccall( (:cublasDdgmm, libcublas), cublasStatus_t, (cublasHandle_t, cublasSideMode_t, Cint, Cint, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint), handle, mode, m, n, A, lda, x, incx, C, ldc))
end
function cublasCdgmm(handle, mode, m, n, A, lda, x, incx, C, ldc)
  statuscheck(ccall( (:cublasCdgmm, libcublas), cublasStatus_t, (cublasHandle_t, cublasSideMode_t, Cint, Cint, Ptr{cuComplex}, Cint, Ptr{cuComplex}, Cint, Ptr{cuComplex}, Cint), handle, mode, m, n, A, lda, x, incx, C, ldc))
end
function cublasZdgmm(handle, mode, m, n, A, lda, x, incx, C, ldc)
  statuscheck(ccall( (:cublasZdgmm, libcublas), cublasStatus_t, (cublasHandle_t, cublasSideMode_t, Cint, Cint, Ptr{cuDoubleComplex}, Cint, Ptr{cuDoubleComplex}, Cint, Ptr{cuDoubleComplex}, Cint), handle, mode, m, n, A, lda, x, incx, C, ldc))
end
function cublasStpttr(handle, uplo, n, AP, A, lda)
  statuscheck(ccall( (:cublasStpttr, libcublas), cublasStatus_t, (cublasHandle_t, cublasFillMode_t, Cint, Ptr{Cfloat}, Ptr{Cfloat}, Cint), handle, uplo, n, AP, A, lda))
end
function cublasDtpttr(handle, uplo, n, AP, A, lda)
  statuscheck(ccall( (:cublasDtpttr, libcublas), cublasStatus_t, (cublasHandle_t, cublasFillMode_t, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cint), handle, uplo, n, AP, A, lda))
end
function cublasCtpttr(handle, uplo, n, AP, A, lda)
  statuscheck(ccall( (:cublasCtpttr, libcublas), cublasStatus_t, (cublasHandle_t, cublasFillMode_t, Cint, Ptr{cuComplex}, Ptr{cuComplex}, Cint), handle, uplo, n, AP, A, lda))
end
function cublasZtpttr(handle, uplo, n, AP, A, lda)
  statuscheck(ccall( (:cublasZtpttr, libcublas), cublasStatus_t, (cublasHandle_t, cublasFillMode_t, Cint, Ptr{cuDoubleComplex}, Ptr{cuDoubleComplex}, Cint), handle, uplo, n, AP, A, lda))
end
function cublasStrttp(handle, uplo, n, A, lda, AP)
  statuscheck(ccall( (:cublasStrttp, libcublas), cublasStatus_t, (cublasHandle_t, cublasFillMode_t, Cint, Ptr{Cfloat}, Cint, Ptr{Cfloat}), handle, uplo, n, A, lda, AP))
end
function cublasDtrttp(handle, uplo, n, A, lda, AP)
  statuscheck(ccall( (:cublasDtrttp, libcublas), cublasStatus_t, (cublasHandle_t, cublasFillMode_t, Cint, Ptr{Cdouble}, Cint, Ptr{Cdouble}), handle, uplo, n, A, lda, AP))
end
function cublasCtrttp(handle, uplo, n, A, lda, AP)
  statuscheck(ccall( (:cublasCtrttp, libcublas), cublasStatus_t, (cublasHandle_t, cublasFillMode_t, Cint, Ptr{cuComplex}, Cint, Ptr{cuComplex}), handle, uplo, n, A, lda, AP))
end
function cublasZtrttp(handle, uplo, n, A, lda, AP)
  statuscheck(ccall( (:cublasZtrttp, libcublas), cublasStatus_t, (cublasHandle_t, cublasFillMode_t, Cint, Ptr{cuDoubleComplex}, Cint, Ptr{cuDoubleComplex}), handle, uplo, n, A, lda, AP))
end
