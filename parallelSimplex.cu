#include <curand.h>
#include <curand_kernel.h>

#define TOTAL_BLOCKS 256;


__device__ int policy(double *currentJob, double * leftChild, double * rightChild)
{
  unsigned long long seed = threadIdx.x;
  unsigned long long sequence = threadIdx.x;
  unsigned long long offset = 0;
  curandState_t state;
  curand_init(seed,sequence,offset,&state);
  *leftChild = *currentJob + 1;
  *rightChild = *currentJob + 2;
  int randomInt = curand(&state) % 10;
  // if ((randomInt) <= 9)
  //   return 0;
  // else
     return 1;
}

__global__ void branchAndBound(int upperBound, int* done1, int* done2,
int* length1, int* length2, double** TableauList1, double** TableauList2)
{
  //__shared__ double* SimplexTableau;
  __shared__ int* length; //why pointer?
  __shared__ int* otherLength;
  __shared__ int* done;
  __shared__ double** tablelist;
  __shared__ double** otherTablelist;
  __shared__ double* currentJob;
  __shared__ int action;
  __shared__ bool localDone;
//Assign shared memory to passed in values.


  //while (*done2 == totalBlocks) //&& (length1 == 0))
   while(!(((*done2 == 256) && (length1 == 0)) || ((*done1 == 256) && (length2 == 0))))//when there are still jobs
  {
    if (*done1 == 256)
    {
      done = done2;
      tablelist = TableauList2;
      otherTablelist = TableauList1;
      length = length2;
      otherLength = length1;
      *done1 = 0;
    }

    else
    {
      done = done1;
      tablelist = TableauList1;
      otherTablelist = TableauList2;
      length = length1;
      otherLength = length2;
      *done2 = 0;
    }

    while (*done != 256)
    {
      if(!localDone)
      {
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
          action = policy(currentJob, leftChild, rightChild);


          if (threadIdx.x == 0)
          {
            if(action == 0) // STOP
              continue;
            else // BRANCH
            {
              atomicAdd(otherLength, 2);
              *otherTablelist[*otherLength - 2] = *leftChild;
              *otherTablelist[*otherLength - 1] = *rightChild;
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
  double **TableauList1, ** TableauList2;
  int * done1 = new int(0);
  int* done2 = new int(1);
  int * length1 = new int(1);
  int * length2 = new int(0);
  branchAndBound<<<256, 1>>> (10, done1, done2, length1, length2, TableauList1, TableauList2);

  return 0;
}
