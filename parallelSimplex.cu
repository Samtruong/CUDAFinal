#include <stdio.h>
#include <stdlib.h>
#include <curand.h>
#include <curand_kernel.h>
#include <iostream>

#define TOTAL_BLOCKS 256;

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

__global__ void branchAndBound(int upperBound, int* done1, int* done2,
int* length1, int* length2, double* TableauList1, double* TableauList2)
{
  // printf("hello1\n");
  //__shared__ double* SimplexTableau;
  __shared__ int* length; //why pointer?
  __shared__ int* otherLength;
  __shared__ int* done;
  __shared__ double* tablelist;
  __shared__ double* otherTablelist;
  __shared__ double currentJob;
  __shared__ int action;
  __shared__ bool localDone;
  // printf("hello2\n");
//Assign shared memory to passed in values.

int numblocks = 2;
  // printf("done2 %d, length1 %d, done1 %d, length2 %d\n", *done2, *length1, *done1, *length2);
   while(!(((*done2 == numblocks) && (*length1 == 0)) || ((*done1 == numblocks) && (*length2 == 0))))//when there are still jobs
  {
     // printf("About to enter while loop\n");
     localDone = 0;
    if (*done1 == numblocks) //start with 2
    {
      done = done2;
      tablelist = TableauList2;
      otherTablelist = TableauList1;
      length = length2;
      otherLength = length1;
      *done1 = 0;
    }

    else //start with 1
    {
      done = done1;
      tablelist = TableauList1;
      otherTablelist = TableauList2;
      length = length1;
      otherLength = length2;
      *done2 = 0;
    }
    while (*done != numblocks)
    {
      if(!localDone)
      {
        // printf("inside if\n");
        int index = blockIdx.x;
        int stride = gridDim.x;
        for (int i = index; i < *length; i+= stride)
        {
          if (threadIdx.x == 0) //Load a job for current block.
          {
            currentJob = tablelist[i];
          }
          __syncthreads();
          double * leftChild, * rightChild;
          leftChild = (double*) malloc(sizeof(double));
          rightChild = (double*) malloc(sizeof (double));
          action = policy(currentJob, leftChild, rightChild);
          // printf("otherLength %d\n", *otherLength);
          // printf("current job %f\n leftChild %f\n rightChild %f\n", currentJob, *leftChild, *rightChild);
          // printf("otherLength %d\n", *otherLength);
          if (threadIdx.x == 0)
          {
            if(action == 0) // BOUND
              continue;
            else // BRANCH
            {
              atomicAdd(otherLength, 2);
              // printf("otherLength %d\n", *otherLength);
              otherTablelist[(*otherLength) - 2] = *leftChild;
              otherTablelist[(*otherLength) - 1] = *rightChild;
              // printf("leftChild %f\n", *leftChild);
              // printf("rightChild %f\n", *rightChild);
            }
          }
        }
        if (threadIdx.x == 0)
        {
          atomicAdd(done, 1);
          localDone = 1;
        }
        __syncthreads();
      }//do operation
      *length = 0;
    }
  }
}


int main ()
{
  double *TableauList1, * TableauList2;
  cudaMallocManaged(&TableauList1, sizeof(double) * 10);
  cudaMallocManaged(&TableauList2, sizeof(double) * 10);
  *TableauList1 = 1.0;
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
  branchAndBound<<<2, 1>>> (10, done1, done2, length1, length2, TableauList1, TableauList2);

  cudaDeviceSynchronize();
  for (int i = 0; i < 10; i++)
  {
    cout << TableauList1[i] << " ";
  }
    cout << endl;

  for (int i = 0; i < 10; i++)
    cout << TableauList2[i] << " ";
    cout << endl;

  cudaFree(done1);
  cudaFree(done2);
  cudaFree(length1);
  cudaFree(length2);
  cudaFree(TableauList1);
  cudaFree(TableauList2);
  return 0;
}
