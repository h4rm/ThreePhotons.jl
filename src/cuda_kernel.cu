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
  //Multiplies the coefficients according to indices list
  __global__ void calculate_coefficient_matrix(const cuFloatComplex* coeff, const int numcoeff, const float* wignerlist, const int* indices, const int indiceslength, const int* PAcombos, const int combolength, float *PA, const int klength)
  {
    int i = threadIdx.x + blockIdx.x * blockDim.x + 1;
    if(i<=combolength)
    {
      const int k1 = PAcombos[IDX2F(1, i, 9)];
      const int k2 = PAcombos[IDX2F(2, i, 9)];
      const int k3 = PAcombos[IDX2F(3, i, 9)];
      const int ki = PAcombos[IDX2F(4, i, 9)];
      const int jstart = PAcombos[IDX2F(8, i, 9)];
      const int mcombos = PAcombos[IDX2F(9, i, 9)];

      const cuFloatComplex* ck1 = &coeff[(k1-1)*numcoeff];
      const cuFloatComplex* ck2 = &coeff[(k2-1)*numcoeff];
      const cuFloatComplex* ck3 = &coeff[(k3-1)*numcoeff];

      float As = 0.0f;
      for(int n=0; n < mcombos; n++){
        const int j = jstart + n;
        const int k1i = indices[IDX2F(1, j, 9)];
        const int k2i = indices[IDX2F(2, j, 9)];
        const int k3i = indices[IDX2F(3, j, 9)];
        As += wignerlist[j-1]*cuCrealf(cuCmulf(ck1[k1i-1], cuCmulf(ck2[k2i-1],ck3[k3i-1])) );
      }

      for(int n=0; n <mcombos; n++){
        PA[IDX2F(jstart+n,ki, indiceslength)] *= As;
      }
    }
  }
}
