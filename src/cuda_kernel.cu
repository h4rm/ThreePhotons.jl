//nvcc -ptx cuda_kernel.cu --gpu-architecture=compute_30 --gpu-code=compute_30
//In case of compilation error: https://github.com/arrayfire/arrayfire/issues/1384
#include "cuComplex.h"
// CUDA runtime
#include "cuda_runtime.h"
#include "stdint.h"

#define IDX2F(i,j,ld) ((((j)-1)*(ld))+((i)-1))
#define IDX2C(i,j,ld) (((j)*(ld))+(i))
// #define IDX3(k1,k2,k3,kcut) (((((k1)-1)*(kcut))+((k2)-1))*(kcut) + ((k3)-1))

extern "C"
{

  //Multiplies the coefficients according to basisindices list
  __global__ void calculate_coefficient_matrix(const cuFloatComplex* coeff, const int numcoeff, const int kcut, const int* basisindices, const int basislen, float *faclist)
  {
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    if(i<basislen)
    {
      int k1i = basisindices[IDX2F(i+1, 1, basislen)];
      int k2i = basisindices[IDX2F(i+1, 2, basislen)];
      int k3i = basisindices[IDX2F(i+1, 3, basislen)];

      int j = 0;
      for(int k1=1; k1<=kcut; k1++){
        const cuFloatComplex* ck1 = &coeff[(k1-1)*numcoeff];
        for(int k2=1; k2<=k1; k2++){
          const cuFloatComplex* ck2 = &coeff[(k2-1)*numcoeff];
          for(int k3=1; k3<=k2; k3++){
            const cuFloatComplex* ck3 = &coeff[(k3-1)*numcoeff];

            // int j = IDX3(k1,k2,k3,kcut);
            faclist[IDX2F(i+1, j+1, basislen)] = cuCrealf( cuCmulf(ck1[k1i-1], cuCmulf(ck2[k2i-1],ck3[k3i-1])) );
            j += 1;
          }
        }
      }
    }
  }

  //Multiplies the coefficients according to basisindices list
  __global__ void calculate_coefficient_matrix_optimized(const cuFloatComplex* coeff, const int* kmapping, const int numcoeff, const int kcut, const int* basisindices, const int basislen, float *faclist, const int klength)
  {
    const int i = threadIdx.x + blockIdx.x * blockDim.x;
    const int j = threadIdx.y + blockIdx.y * blockDim.y;

    if(i<basislen && j < klength)
    {
      const int k1i = basisindices[IDX2F(i+1, 1, basislen)];
      const int k2i = basisindices[IDX2F(i+1, 2, basislen)];
      const int k3i = basisindices[IDX2F(i+1, 3, basislen)];

      const int k1 = kmapping[IDX2F(j+1, 1, klength)];
      const int k2 = kmapping[IDX2F(j+1, 2, klength)];
      const int k3 = kmapping[IDX2F(j+1, 3, klength)];

      const cuFloatComplex* ck1 = &coeff[(k1-1)*numcoeff];
      const cuFloatComplex* ck2 = &coeff[(k2-1)*numcoeff];
      const cuFloatComplex* ck3 = &coeff[(k3-1)*numcoeff];

      faclist[IDX2F(i+1, j+1, basislen)] = cuCrealf( cuCmulf(ck1[k1i-1], cuCmulf(ck2[k2i-1],ck3[k3i-1])) );
    }
  }

}
