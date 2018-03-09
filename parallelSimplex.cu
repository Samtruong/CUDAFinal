#define TOTAL_BLOCKS 256;
__global__ void branchAndBound(int upperBound,int* done1, int* done2,
int* length1, int* length2, double** TableauList1, double** TableauList2)
{
  __shared__ double* SimplexTableau;
  int* length;
  __shared__ int* done;
  __shared__ double** tablelist;
  while(*done1 != TOTAL_BLOCKS || *done2 != TOTAL_BLOCKS)//when there are still jobs
  {
    /*
      FIGURE IT OUT WHICH ARRAY TO USE HERE!!!!
    */
    if(threadIdx.x == 0){atomicAdd(*done,-1);} //signify this block is busy
    for(i = blockIdx.x; i < length; i+=gridDim.x)
    {
      if(threadIdx.x == 0)
      {
        SimplexTableau = TableauList[i];
      }
    }
    if(threadIdx.x == 0){atomicAdd(*done,1);}
  }
}
