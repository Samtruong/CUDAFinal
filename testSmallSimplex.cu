#include <stdio.h>
#include <stdlib.h>
#include <curand.h>
#include <curand_kernel.h>
#include <iostream>

using namespace std;

__global__ void smallSimplex(double* job, int width, int height)
{
  extern __shared__ double shared[];
  double* pivotColumn = &shared[0];
  double* pivotRow = &shared[height];
  double* ratioColumn = &shared[width];


  __shared__ double smallestObj;
  __shared__ int pivotColIndx;
  __shared__ int pivotRowIndx;

  smallestObj = 10^20;
  int index = blockDim.x * blockIdx.x + threadIdx.x;
  int stride = blockDim.x * gridDim.x;
  while(true)
  {
    //  STEP 1
    if(threadIdx.x == 0)
    {
      for(int i = 0; i < width; i++)
      {
        if(smallestObj > job[(height-1)*width+i])
        {
          smallestObj = job[(height-1)*width+i];
          pivotColIndx = i;
        }
      }
    }
    __syncthreads();
    if(smallestObj >= 0){break;}
    // printf("Smallest collumn is %d with value %f \n",pivotColIndx,smallestObj);

    for(int i = index; i < height; i+=stride)
    {
        pivotColumn[i] = job[pivotColIndx + i*width];
        ratioColumn[i] = job[(i+1)*width-1] / pivotColumn[i];
    }

    // printf("The ratio column is:\n");
    // for(int i = 0; i < height; i++)
    // {
    //   printf("%f\n",ratioColumn[i]);
    // }
  //  STEP 2
    smallestObj = 10^20;
    if(threadIdx.x == 0)
    {
      for(int i = 0; i < height-1; i++)
      {
        if(smallestObj > ratioColumn[i])
        {
          smallestObj = ratioColumn[i];
          pivotRowIndx = i;
        }
      }
    }
    __syncthreads();
    // printf("Smallest row is %d with value %f \n",pivotRowIndx,smallestObj);
    for(int i = index; i < width; i+=stride)
    {
      job[width*pivotRowIndx + i] = job[width*pivotRowIndx + i]/pivotColumn[pivotRowIndx];
      pivotRow[i] = job[width*pivotRowIndx + i];
    }
  //  STEP 3
    for (int i = index; i < height; i+= stride)
    {
      if(i == pivotRowIndx){continue;}
      for (int j = 0; j < width; j++)
      {
        job[width * i + j] = job[width * i + j] - pivotColumn[i]* pivotRow[j];
      }
    }
  }

}

int main(int argc, char const *argv[])
{
  double* simplexTable;
  cudaMallocManaged(&simplexTable, sizeof(double)*18);
  // simplexTable[0] = 3.;
  // simplexTable[1] = 4.;
  // simplexTable[2] = 1.;
  // simplexTable[3] = 0.;
  // simplexTable[4] = 0.;
  // simplexTable[5] = 24.;
  //
  // simplexTable[6] = 7.;
  // simplexTable[7] = -4.;
  // simplexTable[8] = 0.;
  // simplexTable[9] = 1.;
  // simplexTable[10] = 0.;
  // simplexTable[11] = 16.;
  //
  // simplexTable[12] = -2.;
  // simplexTable[13] = 3.;
  // simplexTable[14] = 0.;
  // simplexTable[15] = 0.;
  // simplexTable[16] = 1.;
  // simplexTable[17] = 0.;

  simplexTable[0] = 10.;
  simplexTable[1] = 7.;
  simplexTable[2] = 1.;
  simplexTable[3] = 0.;
  simplexTable[4] = 0.;
  simplexTable[5] = 40.;

  simplexTable[6] = 1.;
  simplexTable[7] = 1.;
  simplexTable[8] = 0.;
  simplexTable[9] = 1.;
  simplexTable[10] = 0.;
  simplexTable[11] = 5.;

  simplexTable[12] = -17.;
  simplexTable[13] = -12.;
  simplexTable[14] = 0.;
  simplexTable[15] = 0.;
  simplexTable[16] = 1.;
  simplexTable[17] = 0.;
  int height = 3;
  int width = 6;
  int sharedMemory = sizeof(double)*(height*2 + width);
  smallSimplex<<<1,256,sharedMemory>>>(simplexTable,width,height);
  cudaDeviceSynchronize();
  for(int i = 0; i <18; i++)
  {
    cout<<simplexTable[i]<<" "<<"|"<<" ";
    if(i%6 == 5){cout<<endl;}
  }
  return 0;
}
