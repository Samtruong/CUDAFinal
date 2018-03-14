#include <stdio.h>
#include <stdlib.h>
#include <curand.h>
#include <curand_kernel.h>
#include <iostream>
using namespace std;

__device__ int policy(double currentJob, double * leftChild, double * rightChild)
{
  unsigned long long seed = threadIdx.x;
  unsigned long long sequence = threadIdx.x;
  unsigned long long offset = 0;
  curandState_t state;
  curand_init(seed,sequence,offset,&state);
  if (currentJob == 1)
  {
    *leftChild = 2.;
    *rightChild = 3.;
  }
  else if (currentJob == 2)
  {
    *leftChild = 4.;
    *rightChild = 5.;
  }
  else if (currentJob == 3)
  {
    *leftChild = 6.;
    *rightChild = 7.;
  }
  int randomInt = curand(&state) % 10;
  // if ((randomInt) <= 9)
  //   return 0;
  // else
  if (currentJob == 1 || currentJob == 2 || currentJob == 3)
    return 1;
  else
    return 0;
}


__global__ void branchAndBound(int numblocks,int upperBound, int* doneCPU, int* lengthCPU, int* otherLengthCPU,
  double* tableListCPU, double* othertableListCPU)
{
  __shared__ int* length;
  __shared__ int* otherLength;
  __shared__ int* done;
  __shared__ double* tableList;
  __shared__ double* othertableList;
  __shared__ double currentJob;
  __shared__ int action;
  __shared__ bool localDone;
  __shared__ int index;
  __shared__ int stride;

  index = blockIdx.x;
  stride = gridDim.x;

if (threadIdx.x == 0)
{

  length = lengthCPU;
  otherLength = otherLengthCPU;
  done = doneCPU;
  tableList = tableListCPU;
  othertableList = othertableListCPU;
  printf("Thread 0 done copying.\n");
}
  // while(!(((*done2 == numblocks) && (*length1 == 0)) || ((*done1 == numblocks) && (*length2 == 0))))
  // {//if there are jobs in either array, we will keep the kernel alive.
  //
  //   localDone = false;
  //   /*The following if else statements are used to determined which array to process */
  //   if (*done1 == numblocks){
  //     done = done2;
  //     tableList = tableList2;
  //     othertableList = tableList1;
  //     length = length2;
  //     otherLength = length1;
  //   }
  //   else{
  //     done = done1;
  //     tableList = tableList1;
  //     othertableList = tableList2;
  //     length = length1;
  //     otherLength = length2;
  //   }
  //   *otherLength = 0;
  //   *done1 = 0;
  //   *done2 = 0;

    while (*done != numblocks)
    {//if not everyone is done with the current array.
      if(localDone == false)
      {//if I am not done with the current array
        for (int i = index; i < *length; i+= stride)
        {//fetch all jobs I am suppose to process
          if (threadIdx.x == 0){currentJob = tableList[i];}
          __syncthreads();

          double * leftChild, * rightChild;
          leftChild = (double*) malloc(sizeof(double));
          rightChild = (double*) malloc(sizeof (double));
          action = policy(currentJob, leftChild, rightChild);

          if (threadIdx.x == 0)
          {
            if(action == 0) // BOUND
              continue;
            else // BRANCH
            {
              atomicAdd(otherLength, 2);
              othertableList[(*otherLength) - 2] = *leftChild;
              othertableList[(*otherLength) - 1] = *rightChild;
            }
          }
          __syncthreads();
        }//after for loop, I will run out of job, I will wait for other to be done
        if (threadIdx.x == 0){atomicAdd(done, 1);localDone = true; printf("atomic add performed on done %d\n", *done);}
        __syncthreads();
      }
    }
  //}
}


int main ()
{
  double *tableList1, * tableList2;
  cudaMallocManaged(&tableList1, sizeof(double) * 10);
  cudaMallocManaged(&tableList2, sizeof(double) * 10);
  *tableList1 = 1.0;
  int * done1;
  cudaMallocManaged(&done1, sizeof(int));
  *done1 = 0;
  int* done2;
  cudaMallocManaged(&done2, sizeof(int));
  *done2 = 1;
  int * length1;
  cudaMallocManaged(&length1, sizeof(int));
  *length1 = 1;
  int * length2;
  cudaMallocManaged(&length2, sizeof(int));
  *length2 = 0;
  int numblocks = 1;

  //CPU Logic
  int *done;
  cudaMallocManaged(&done, sizeof(int));
  double *tableList;
  cudaMallocManaged(&tableList, sizeof(double) * 10);
  double *othertableList;
  cudaMallocManaged(&othertableList, sizeof(double) * 10);
  int *length;
  cudaMallocManaged(&length, sizeof(int));
  int *otherLength;
  cudaMallocManaged(&otherLength, sizeof(int));

  while(!(((*done2 == numblocks) && (*length1 == 0)) || ((*done1 == numblocks) && (*length2 == 0))))
  {
    if (*done1 == numblocks){
      done = done2;
      tableList = tableList2;
      othertableList = tableList1;
      length = length2;
      otherLength = length1;
    }
    else{
      done = done1;
      tableList = tableList1;
      othertableList = tableList2;
      length = length1;
      otherLength = length2;
    }
    *otherLength = 0;
    *done1 = 0;
    *done2 = 0;
     printf("Kernel called\n");
    branchAndBound<<<numblocks, 1>>> (numblocks,10, done, length, otherLength, tableList, othertableList);
    cudaDeviceSynchronize();
  }
  cudaDeviceSynchronize();
  // branchAndBound<<<numblocks, 1>>> (numblocks,10, done1, done2, length1, length2, tableList1, tableList2);
  // cudaDeviceSynchronize();
  for (int i = 0; i < 10; i++)
  {
    cout << tableList1[i] << " ";
  }
    cout << endl;

  for (int i = 0; i < 10; i++)
    cout << tableList2[i] << " ";
    cout << endl;

  cudaFree(done1);
  cudaFree(done2);
  cudaFree(length1);
  cudaFree(length2);
  cudaFree(tableList1);
  cudaFree(tableList2);
  cudaFree(done);
  cudaFree(tableList);
  cudaFree(othertableList);
  cudaFree(length);
  cudaFree(otherLength);
  return 0;
}
