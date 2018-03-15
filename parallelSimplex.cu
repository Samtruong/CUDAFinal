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

__device__ void smallSimplexSolver(double* job, int width, int height,
   double* pivotColumn, double* pivotRow, double* ratioColumn)
{
  __shared__ double smallestObj;
  __shared__ int pivotColIndx;
  __shared__ int pivotRowIndx;
  smallestObj = 10^20;
  int index = threadIdx.x;
  int stride = blockDim.x;
  // while(true)
  // {
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
    // if(smallestObj >= 0){break;}
    printf("Smallest collumn is %d with value %f \n",pivotColIndx,smallestObj);

    for(int i = index; i < height; i+=stride)
    {
        pivotColumn[i] = job[pivotColIndx + i*width];
        ratioColumn[i] = job[(i+1)*width-1] / pivotColumn[i];
    }

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
  //}
}


__global__ void branchAndBound(int numblocks, int* doneCPU, int* lengthCPU, int* otherLengthCPU,
  double* tableListCPU, double* othertableListCPU, int curWidth, int curHeight)
{
  __shared__ int* length;
  __shared__ double* tableList;
  __shared__ double* othertableList;
  __shared__ double* currentJob;
  __shared__ int action;
  __shared__ bool localDone;
  __shared__ int index;
  __shared__ int stride;
  __shared__ double* original;

  //simplex variables:

  index = blockIdx.x;
  stride = gridDim.x;

  if(threadIdx.x == 0)
  {
    length = lengthCPU;
    tableList = tableListCPU;
    othertableList = othertableListCPU;
  }
  __syncthreads();
  localDone = false;
  while (*doneCPU != numblocks)
  {//if not everyone is done with the current array.
    if(localDone == false)
    {//if I am not done with the current array
      for (int i = index; i < *length; i+= stride)
      {//fetch all jobs I am suppose to process
        __shared__ double* pivotRow;
        __shared__ double* pivotCol;
        __shared__ double* ratioCol;
        // c
        if(threadIdx.x == 0)
        {
          currentJob = &tableList[i*curWidth*curHeight];
          pivotCol = (double*) malloc(sizeof(double) * curHeight);
          pivotRow = (double*) malloc(sizeof(double) * curWidth);
          ratioCol = (double*) malloc(sizeof(double) * curHeight);
        }
        __syncthreads();
        original = (double*) malloc(sizeof(double) * curWidth * curHeight);
        for(int i = threadIdx.x; i < curWidth * curHeight; i+=stride)
        {
          original[i] = currentJob[i];
        }
        smallSimplexSolver(currentJob,curWidth,curHeight,pivotCol,pivotRow,ratioCol);
        for(int i = 0; i <curWidth*curHeight; i++)
        {
          printf("%f | ",(currentJob[i]));
          if(i%curWidth == curWidth-1){printf("\n");}
        }

        double  breakpoint;
        int leavingVariable;
        double* leftChild;
        double* rightChild;
        if(threadIdx.x == 0)
        {
          free(pivotCol);
          free(pivotRow);
          free(ratioCol);
          if(currentJob[curWidth*curHeight - 1] > 68.3)
          {
            action = true;
            for(int i = 0; i < curHeight - 1; i++)
            {
              if(currentJob[i + i*curWidth] == 1.)
              {
                breakpoint = currentJob[i*curHeight + curWidth - 1];
                leavingVariable = i;
                break;
              }
            }
          }
          else{action = false;}
        }
        __syncthreads();
        if(!action){continue;} //BOUND
        else //BRANCH
        {
          leftChild = (double*) malloc(sizeof(double) * (curHeight + 1) * (curWidth + 1));
          rightChild = (double*) malloc(sizeof(double) * (curHeight + 1) * (curWidth + 1));
          for(int i = threadIdx.x; i <(curHeight + 1) * (curWidth + 1); i+= blockDim.x)
          {
            leftChild[i] = 0.;
            rightChild[i] = 0.;
          }
          for(int i = threadIdx.x; i <curHeight; i+= blockDim.x)
          {
            for(int j = 0; j < curWidth-1; j++)
            {
              leftChild[i*(curWidth+1) + j] = original[i*curWidth + j];
              rightChild[i*(curWidth+1) + j] = original[i*curWidth + j];
            }
          }
          for(int i = threadIdx.x; i < curHeight; i+=blockDim.x)
          {
            leftChild[(curWidth) + i*(curWidth+1)] = original[(curWidth-1) + i*curWidth];
            rightChild[(curWidth) + i*(curWidth+1)] = original[(curWidth-1) + i*curWidth];
          }
          if(threadIdx.x == 0)
          {
            leftChild[leavingVariable+(curWidth+1)*curHeight] = 1.;
            leftChild[curWidth-1+(curWidth+1)*curHeight] = 1.;
            leftChild[(curWidth+1)*(curHeight+1)-1] = floor(breakpoint);
            rightChild[leavingVariable+(curWidth+1)*curHeight] = 1.;
            rightChild[curWidth-1+(curWidth+1)*curHeight] = -1.;
            rightChild[(curWidth+1)*(curHeight+1)-1] = -1*ceil(breakpoint);
          }
          double swap;
          for(int i = threadIdx.x; i < (curWidth+1); i+= blockDim.x)
          {
            swap = leftChild[(curWidth+1) * (curHeight+1) - 2*(curWidth+1) + i];
            leftChild[(curWidth+1) * (curHeight+1) - 2*(curWidth+1) + i] = leftChild[(curWidth+1) * (curHeight+1) - (curWidth+1) + i];
            leftChild[(curWidth+1) * (curHeight+1) - (curWidth+1) + i] = swap;

            swap = rightChild[(curWidth+1) * (curHeight+1) - 2*(curWidth+1) + i];
            rightChild[(curWidth+1) * (curHeight+1) - 2*(curWidth+1) + i] = rightChild[(curWidth+1) * (curHeight+1) - (curWidth+1) + i];
            rightChild[(curWidth+1) * (curHeight+1) - (curWidth+1) + i] = swap;
          }
          // for(int i = 0; i <28; i++)
          // {
          //   printf("%f | ",leftChild[i]);
          //   if(i%7 == 6){printf("\n");}
          // }
          // printf("\n");
          // for(int i = 0; i <28; i++)
          // {
          //   printf("%f | ",rightChild[i]);
          //   if(i%7 == 6){printf("\n");}
          // }
          atomicAdd(otherLengthCPU,2);
          for(int i = threadIdx.x; i < (curWidth+1)*(curHeight+1);i+=blockDim.x)
          {
            othertableList[(*otherLengthCPU - 2)*(curWidth+1)*(curHeight+1) + i] = leftChild[i];
            othertableList[(*otherLengthCPU - 1)*(curWidth+1)*(curHeight+1) + i] = rightChild[i];
          }
        }

      }//after for loop, I will run out of job, I will wait for other to be done
      if (threadIdx.x == 0){atomicAdd(doneCPU, 1);localDone = true;}
      __syncthreads();
    }
    __syncthreads();
  }
  //}
}


int main ()
{
  int numblocks = 1;
  int safeCount = 0;
  int height = 3;
  int width = 6;

  double *tableList1, * tableList2;
  int * done1;
  int* done2;
  int * length1;
  int * length2;


  cudaMallocManaged(&tableList1, sizeof(double) * 20 * width * height);
  cudaMallocManaged(&tableList2, sizeof(double) * 20 * (width + 1) * (height + 1));
  cudaMallocManaged(&done1, sizeof(int));
  cudaMallocManaged(&done2, sizeof(int));
  cudaMallocManaged(&length1, sizeof(int));
  cudaMallocManaged(&length2, sizeof(int));


  tableList1[0] = 10.;
  tableList1[1] = 7.;
  tableList1[2] = 1.;
  tableList1[3] = 0.;
  tableList1[4] = 0.;
  tableList1[5] = 40.;

  tableList1[6] = 1.;
  tableList1[7] = 1.;
  tableList1[8] = 0.;
  tableList1[9] = 1.;
  tableList1[10] = 0.;
  tableList1[11] = 5.;

  tableList1[12] = -17.;
  tableList1[13] = -12.;
  tableList1[14] = 0.;
  tableList1[15] = 0.;
  tableList1[16] = 1.;
  tableList1[17] = 0.;
  *length1 = 1;
  *length2 = 0;
  *done1 = 0;
  *done2 = numblocks;
  //CPU Logic
  int *done;
  double *tableList;
  double *othertableList;
  int *length;
  int *otherLength;

  cudaMallocManaged(&done, sizeof(int));
  cudaMallocManaged(&length, sizeof(int));
  cudaMallocManaged(&otherLength, sizeof(int));

  while(!(((*done2 == numblocks) && (*length1 == 0)) || ((*done1 == numblocks) && (*length2 == 0)) || safeCount >= 10))
  {
    cudaMallocManaged(&tableList, sizeof(double) * 20 * width * height);
    cudaMallocManaged(&othertableList, sizeof(double) * 20 * (width + 1) * (height + 1));
    safeCount++;
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

    // printf("Kernel called\n");
    branchAndBound<<<numblocks,1>>> (numblocks, done, length, otherLength, tableList, othertableList,
      width,height);
    cudaDeviceSynchronize();
    width++;
    height++;
    cudaFree(tableList);
    cudaFree(othertableList);
  }
  // for(int i = 0; i <28; i++)
  // {
  //   printf("%f | ",(othertableList[i]));
  //   if(i%7 == 6){printf("\n");}
  // }
  // printf("\n");
  // for(int i = 0; i <28; i++)
  // {
  //   printf("%f | ",othertableList[28+i]);
  //   if(i%7 == 6){printf("\n");}
  // }
  // cout <<endl;
  // for(int i = 0; i <18; i++)
  // {
  //   printf("%f | ",tableList1[i]);
  //   if(i%6 == 5){printf("\n");}
  // }

  // branchAndBound<<<numblocks, 1>>> (numblocks,10, done1, done2, length1, length2, tableList1, tableList2);
  // cudaDeviceSynchronize();
  // for (int i = 0; i < 10; i++)
  // {
  //   cout << tableList1[i] << " ";
  // }
  //   cout << endl;
  //
  // for (int i = 0; i < 10; i++)
  //   cout << tableList2[i] << " ";
  //   cout << endl;

  cudaFree(done1);
  cudaFree(done2);
  cudaFree(length1);
  cudaFree(length2);
  cudaFree(tableList1);
  cudaFree(tableList2);
  cudaFree(done);
  cudaFree(length);
  cudaFree(otherLength);
  return 0;
}
