#include <stdio.h>
#include <stdlib.h>
#include <set>
#include <sstream>
#include <string>
#include <fstream>
#include <iostream>
#include <cstring>
#include <thrust/scan.h>
#include <thrust/execution_policy.h>
#include <curand.h>
#include <curand_kernel.h>

using namespace std;


__global__ void GraphGenerator(int* matrix,int* dimension, int* address, int* h_graph, int V)
{
  int index = threadIdx.x + blockDim.x * blockIdx.x;
  int stride = blockDim.x * gridDim.x;
  for(int i = index; i < V; i += stride)
  {
    int a = address[i];
    int j = 0;
    for (int k = 0; k < V; k++)
    {
      if (matrix[i*V + k])
      {
        h_graph[a + j] = k;
        j++;
      }
    }
  }
}

__global__ void DimensionGenerator(int* matrix, int* dimension, int* address, int V)
{
  int index = threadIdx.x + blockDim.x * blockIdx.x;
  int stride = blockDim.x * gridDim.x;
  for(int i = index; i < V; i += stride)
  {
    for (int j = 0; j < V; j++)
    {
      if(matrix[i*V + j])
      {
        dimension[i]++;
      }
    }
  }
  __syncthreads();
}


//================================Utility Functions=======================================
void CountColors(int V,int length, int* color, int &minColors, int &minIndex)
{
	//int minColors = INT_MAX;
	//int minIndex;
   int *num_colors;
	num_colors = (int*) malloc(sizeof(int) * length);
	for (int i = 0; i < length; i++)
	{
		num_colors[i] = 0;
	}
   set<int> seen_colors;

   for (int i = 0; i < length; i++) {
      if (seen_colors.find(color[i]) == seen_colors.end())
      {
         seen_colors.insert(color[i]);
         num_colors[i/V]++;
      }
      if(i%V==V-1)
      {
        //cout<<num_colors[i/V]<<endl;
	if (num_colors[i/V] < minColors)
	{
		minColors = num_colors[i/V];
		minIndex = i / V;
	}
        seen_colors.clear();
        //num_colors = 0;
      }
   }
}

bool IsValidColoring(int* graph, int V, int* color)
{
   for (int i = 0; i < V; i++) {
      for (int j = 0; j < V; j++) {
         if (graph[i * V + j]) {
            if (i != j && color[i] == color[j]) {
               printf("Vertex %d and Vertex %d are connected and have the same color %d\n", i, j, color[i]);
               return false;
            }
            if (color[i] < 1) {
               printf("Vertex %d has invalid color %d\n", i, color[i]);

            }
         }
      }
   }

   return true;
}

//Load raw .co data
void getDimension(const char filename[], int* V)
{
   string line;
   ifstream infile(filename);
   if (infile.fail()) {
      printf("Failed to open %s\n", filename);
      return;
   }

   int num_rows;

   while (getline(infile, line))
   {
      istringstream iss(line);
      string s;
      iss >> s;
      if (s == "p") {
         iss >> s; // read string "edge"
         iss >> num_rows;
         *V = num_rows;
         break;
      }
   }
   infile.close();
}

void ReadColFile(const char filename[], int* graph, int V)
{
   string line;
   ifstream infile(filename);
   if (infile.fail()) {
      printf("Failed to open %s\n", filename);
      return;
   }

   while (getline(infile, line)) {
      istringstream iss(line);
      string s;
      int node1, node2;
      iss >> s;
      if (s != "e")
         continue;

      iss >> node1 >> node2;

      // Assume node numbering starts at 1
      (graph)[(node1 - 1) * V + (node2 - 1)] = 1;
      (graph)[(node2 - 1) * V + (node1 - 1)] = 1;
   }
   infile.close();
}

//print graph Matrix
void PrintMatrix(int* matrix, int M, int N) {
   for (int row=0; row<M; row++)
   {
      for(int columns=0; columns<N; columns++)
      {
         printf("%i", matrix[row * N + columns]);
      }
      printf("\n");
   }
}


// Read MatrixMarket graphs
// Assumes input nodes are numbered starting from 1
void ReadMMFile(const char filename[], bool** graph, int* V)
{
   string line;
   ifstream infile(filename);
   if (infile.fail()) {
      printf("Failed to open %s\n", filename);
      return;
   }

   // Reading comments
   while (getline(infile, line)) {
      istringstream iss(line);
      if (line.find('%') == string::npos)
         break;
   }

   // Reading metadata
   istringstream iss(line);
   int num_rows, num_cols, num_edges;
   iss >> num_rows >> num_cols >> num_edges;

   *graph = new bool[num_rows * num_rows];
   memset(*graph, 0, num_rows * num_rows * sizeof(bool));
   *V = num_rows;

   // Reading nodes
   while (getline(infile, line)) {
      istringstream iss(line);
      int node1, node2, weight;
      iss >> node1 >> node2 >> weight;

      // Assume node numbering starts at 1
      (*graph)[(node1 - 1) * num_rows + (node2 - 1)] = true;
      (*graph)[(node2 - 1) * num_rows + (node1 - 1)] = true;
   }
   infile.close();
}


//Constraints
int* const2(int numVertices)
{
  int *a;
  a = (int*) malloc(sizeof(int) * numVertices);
  for (int i = 0; i < numVertices; i++)
    a[i] = 1;
  return a;
}
//===================================Main=======================================
void GraphColoringGPU(const char filename[])
{
  int * matrix;
  int * h_graph;
  int * sequence;
  int * dimension;
  int * address;
  int * result;
  int V;



  if (string(filename).find(".col") != string::npos)
  {
    getDimension(filename, &V);
    cudaError_t result = cudaMallocManaged(&matrix,sizeof(int)*V*V);
    ReadColFile(filename,matrix,V);
  }
  /*
  else if (string(filename).find(".mm") != string::npos)
     ReadMMFile(filename, matrix, V);*/


  cudaMallocManaged(&sequence, sizeof(int) * V );
  cudaMallocManaged(&dimension,sizeof(int)*V);
  cudaMallocManaged(&address,sizeof(int)*V);
  cudaMallocManaged(&result, sizeof(int) *V);


  DimensionGenerator<<<256,1024>>>(matrix,dimension,address,V);
  cudaDeviceSynchronize();
  thrust::exclusive_scan(thrust::host,dimension,&dimension[V],address);
  cudaMallocManaged(&h_graph,sizeof(int)* (dimension[V-1]+address[V-1]));

  GraphGenerator<<<256,1024>>>(matrix,dimension,address,h_graph,V);
  cudaDeviceSynchronize();


  cudaFree(h_graph);
  cudaFree(dimension);
  cudaFree(sequence);
  cudaFree(address);
  cudaFree(matrix);
  cudaFree(result);

}

int main(int argc, char const *argv[]) {

  GraphColoringGPU(argv[1]);

  return 0;
}
