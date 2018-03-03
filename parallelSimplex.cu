#include <iostream>
using namespace std;

__global__ void ratioCalculator(double** SimplexTableau, double* theta,
                                                double* Columnk, int k)
{
  int index = blockDim.x*blockIdx.x + threadIdx.x;
  double w = SimplexTableau[index][k];
  Columnk[index] = w;
  theta[index] = SimplexTableau[index][1]/w;
}

__global__ void normalizePivotRow( double** SimplexTableau, double* k, int k, int r)
{
  int index = blockDim.x*blockIdx.x + threadIdx.x;
  __shared__ double w;
  if(threadIdx.x == 0){w = Columnk[r];}
  __syncthreads();
  SimplexTableau[r][index] = SimplexTableau[r][index]/w;
}

__global__ void updateAllRows(double** SimplexTableau, double* Columnk, int k, int r)
{
  int xIndex = blockDim.x*blockIdx.x + threadIdx.x;
  int yIndex = blockDim.y*blockIdx.y + threadIdx.y;
  __shared__ double w[16];
  if(threadIdx.y == 0 && threadIdx.x <16)
  {
    w[threadIdx.x] = Columnk[blockIdx.y*blockDim.y+threadIdx.x];
  }
  __syncthreads();
  if(yIndex == r) return;
  SimplexTableau[yIndex][xIndex] = SimplexTableau[yIndex][xIndex] - w[threadIdx.y]
                                  *SimplexTableau[r][xIndex];
}

__global__ void updateSimplexTableau(double** SimplexTableau, double* Columnk, int k, int r)
{
  int index = blockDim.x*blockIdx.x + threadIdx.x;
  __shared__ double w;
  if(threadIdx.x == 0){w = Columnk[r]}
  __syncthreads();
  SimplexTableau[index][k] = -Columnk[index]/w;
  if(index == r){SimplexTableau[index][k]=1/w;}
}
