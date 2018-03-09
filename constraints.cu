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

//
// __global__ void GraphGenerator(int* matrix,int* dimension, int* address, int* h_graph, int V)
// {
//   int index = threadIdx.x + blockDim.x * blockIdx.x;
//   int stride = blockDim.x * gridDim.x;
//   for(int i = index; i < V; i += stride)
//   {
//     int a = address[i];
//     int j = 0;
//     for (int k = 0; k < V; k++)
//     {
//       if (matrix[i*V + k])
//       {
//         h_graph[a + j] = k;
//         j++;
//       }
//     }
//   }
// }
//
// __global__ void DimensionGenerator(int* matrix, int* dimension, int* address, int V)
// {
//   int index = threadIdx.x + blockDim.x * blockIdx.x;
//   int stride = blockDim.x * gridDim.x;
//   for(int i = index; i < V; i += stride)
//   {
//     for (int j = 0; j < V; j++)
//     {
//       if(matrix[i*V + j])
//       {
//         dimension[i]++;
//       }
//     }
//   }
//   __syncthreads();
// }
//
//
// //================================Utility Functions=======================================
// void CountColors(int V,int length, int* color, int &minColors, int &minIndex)
// {
// 	//int minColors = INT_MAX;
// 	//int minIndex;
//    int *num_colors;
// 	num_colors = (int*) malloc(sizeof(int) * length);
// 	for (int i = 0; i < length; i++)
// 	{
// 		num_colors[i] = 0;
// 	}
//    set<int> seen_colors;
//
//    for (int i = 0; i < length; i++) {
//       if (seen_colors.find(color[i]) == seen_colors.end())
//       {
//          seen_colors.insert(color[i]);
//          num_colors[i/V]++;
//       }
//       if(i%V==V-1)
//       {
//         //cout<<num_colors[i/V]<<endl;
// 	if (num_colors[i/V] < minColors)
// 	{
// 		minColors = num_colors[i/V];
// 		minIndex = i / V;
// 	}
//         seen_colors.clear();
//         //num_colors = 0;
//       }
//    }
// }
//
// bool IsValidColoring(int* graph, int V, int* color)
// {
//    for (int i = 0; i < V; i++) {
//       for (int j = 0; j < V; j++) {
//          if (graph[i * V + j]) {
//             if (i != j && color[i] == color[j]) {
//                printf("Vertex %d and Vertex %d are connected and have the same color %d\n", i, j, color[i]);
//                return false;
//             }
//             if (color[i] < 1) {
//                printf("Vertex %d has invalid color %d\n", i, color[i]);
//
//             }
//          }
//       }
//    }
//
//    return true;
// }
//
// //Load raw .co data
// void getDimension(const char filename[], int* V)
// {
//    string line;
//    ifstream infile(filename);
//    if (infile.fail()) {
//       printf("Failed to open %s\n", filename);
//       return;
//    }
//
//    int num_rows;
//
//    while (getline(infile, line))
//    {
//       istringstream iss(line);
//       string s;
//       iss >> s;
//       if (s == "p") {
//          iss >> s; // read string "edge"
//          iss >> num_rows;
//          *V = num_rows;
//          break;
//       }
//    }
//    infile.close();
// }


// //print graph Matrix
// void PrintMatrix(int* matrix, int M, int N) {
//    for (int row=0; row<M; row++)
//    {
//       for(int columns=0; columns<N; columns++)
//       {
//          printf("%i", matrix[row * N + columns]);
//       }
//       printf("\n");
//    }
// }


// Read MatrixMarket graphs
// // Assumes input nodes are numbered starting from 1
// void ReadMMFile(const char filename[], bool** graph, int* V)
// {
//    string line;
//    ifstream infile(filename);
//    if (infile.fail()) {
//       printf("Failed to open %s\n", filename);
//       return;
//    }
//
//    // Reading comments
//    while (getline(infile, line)) {
//       istringstream iss(line);
//       if (line.find('%') == string::npos)
//          break;
//    }
//
//    // Reading metadata
//    istringstream iss(line);
//    int num_rows, num_cols, num_edges;
//    iss >> num_rows >> num_cols >> num_edges;
//
//    *graph = new bool[num_rows * num_rows];
//    memset(*graph, 0, num_rows * num_rows * sizeof(bool));
//    *V = num_rows;
//
//    // Reading nodes
//    while (getline(infile, line)) {
//       istringstream iss(line);
//       int node1, node2, weight;
//       iss >> node1 >> node2 >> weight;
//
//       // Assume node numbering starts at 1
//       (*graph)[(node1 - 1) * num_rows + (node2 - 1)] = true;
//       (*graph)[(node2 - 1) * num_rows + (node1 - 1)] = true;
//    }
//    infile.close();
// }


//Constraints=======================================================================================
//
// int* const2(int numVertices) //This requires transpose.
// {
//   int *toRet;
//   toRet = (int*) malloc(sizeof(int) * numVertices);
//   for (int i = 0; i < numVertices; i++)
//     toRet[i] = 1;
//   return toRet;
// }

// int* allConsts(int * matrix, int numVertices)
// {
//   int *toRet;
//   int numEdges = 0;
//   for (int i = 0; i < numVertices * numVertices; i++)
//   {
//     if(matrix[i])
//       numEdges++;
//   }
//
//   cudaMallocManaged(*toRet, sizeof((int) * ((numEdges+5)*numVertices)* (4*numVertices+4));
//   // toRet = (int*) malloc(sizeof(int) * ((numEdges+5)*numVertices)* (4*numVertices+4));//Y dimension is 6*numVertices and X dimension is 2*numVertices;
//
//   int row;
//   int col;
//   //Constraint 1
//   for (int i = 0; i < numEdges * numVertices; i++)
//   {
//
//   }
//
//   //Constraint 2
//   int startConstraint2 = i % (4*numVertices + 4); //The row at which const 2 matches.
//   for (; i < numVertices; i++)
//   {
//     row = i % (4*numVertices + 4);
//     col = i / (4*numVertices + 4);
//     toRet[]
//     for (int j = 0; j < 4; j++)
//     {
//       toRet[]
//     }
//   }
// }

// int* const1(int * matrix, int numVertices)
// {
//   //Find # of edges.
//   int * toRet;
//   int numEdges = 0;
//   for (int i = 0; i < numVertices * numVertices; i++)
//   {
//     if(matrix[i])
//       numEdges++;
//   }
//   toRet = (int*) malloc(sizeof(int) * (numVertices+numEdges));
//
//   for (int i = 0; i < numVertices+numEdges; i++)
//   {
//     toRet[i] = 0;
//   }
//   //Populate the matrix to return.
//   numEdges = 0;
//   int row;
//   int col;
//   for (int i = 0; i < numVertices * numVertices; i++)
//   {
//
//     if (matrix[i])
//     {
//       row = i % numVertices;
//       col = i / numVertices;
//       toRet[numEdges* numVertices + row] = 1;
//       toRet[numEdges* numVertices + col] = 1;
//       numEdges++;
//
//     }
//
//   }
//   return toRet;
// }
void readEdgesPosition(const char filename[], int* edgeList)
{
   string line;
   ifstream infile(filename);
   if (infile.fail()) {
      printf("Failed to open %s\n", filename);
      return;
   }

   int i=0;
   while (getline(infile, line)) {
      istringstream iss(line);
      string s;
      int node1, node2;
      iss >> s;
      if (s != "e")
         continue;

      iss >> node1 >> node2;

      edgeList[2*i] = node1;
      edgeList[2*i+1] = node2;
      i++;

   }
   infile.close();
}

void getInfo(const char filename[], int* numEdges, int* numVertices)
{
   string line;
   ifstream infile(filename);
   if (infile.fail()) {
      printf("Failed to open %s\n", filename);
      return;
   }

   while (getline(infile, line)) {
      istringstream iss(line);
      string s,node1;
      iss >> s;
      if (s != "p")
         continue;

      iss >> node1 >> *numVertices >> *numEdges;

      // Assume node numbering starts at 1
   }
   infile.close();
}



__global__ void constraint1 (double *simplexTable, int* edgeList, int numColors, int numEdges, int numVertices)
{
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  int stride = blockDim.x * gridDim.x;
  for(int i = index; i < numEdges; i+=stride)
  {
    for(int j = 0; j < numColors; j++)
    {
      simplexTable[(numColors*i+j)*(numVertices*numColors+numColors)+(edgeList[2*i]-1) * numColors+j] = 1.0;
      simplexTable[(numColors*i+j)*(numVertices*numColors+numColors)+(edgeList[2*i+1]-1) * numColors+j] = 1.0;
    }
  }
}

__global__ void constraint2 (double *simplexTable, int numVertices, int numColors)
{
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  int stride = blockDim.x * gridDim.x;

  for (int i = index; i < numVertices; i += stride)
  {
    for (int j = 0; j < numColors; j++)
    {
      simplexTable[i*(numVertices*numColors+numColors)+numColors*i+j] = 1.0;
    }
  }
}

__global__ void constraint3 (double *simplexTable, int numVertices, int numColors)
{
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  int stride = blockDim.x * gridDim.x;

  for (int i = index; i < numVertices*numColors; i += stride)
  {
    simplexTable[i*(numVertices*numColors+numColors) + i] = 1.0;
    simplexTable[i*(numVertices*numColors+numColors) + numColors * numVertices + (i%numColors)] = -1.0;
  }
}

//===================================Main=======================================
void constraintGenerator(const char filename[])
{
  //cudaMallocManaged(&simplexTable, sizeof((int) * ((numEdges+5)*numVertices)* (4*numVertices+4))
  int numColors = 2;
  int numEdges;
  int numVertices;
  if (string(filename).find(".col") != string::npos)
  {
    getInfo(filename, &numEdges, &numVertices);
  }

  int* edgeList;
  cudaMallocManaged(&edgeList,sizeof(int)*2*numEdges);
  readEdgesPosition(filename,edgeList);

  int * simp;
  cudaMallocManaged(&simp, sizeof(int) *numEdges*2*10);
  for (int i = 0; i < 16*20; i++)
  {
    simp[i] = 0;
  }
   double* simplexTable;
   cudaMallocManaged(&simplexTable, (numColors*numVertices + numVertices) * (numEdges*numColors + numVertices + numColors*numVertices));
  constraint1<<<1,1>>>(simplexTable, edgeList,numColors,numEdges,numVertices);
  constraint2<<<1,1>>>(simplexTable + ((numColors*numVertices + numColors)*(numEdges*numColors)), numVertices, numColors);
  constraint3<<<1,1>>>(simplexTable + ((numColors*numVertices + numColors)*(numEdges*numColors + numVertices)) , numVertices, numColors);
  cudaDeviceSynchronize();
  for (int i = 0; i < (numColors*numVertices + numColors)*(numEdges*numColors + numVertices + numColors*numVertices); i++)
  {
    cout << simplexTable[i] << " ";
    if (i%10==9) cout << endl;
  }
  cout<<endl;
  cudaFree(simp);
  cudaFree(simplexTable);

}

int main(int argc, char const *argv[]) {


  constraintGenerator(argv[1]);
  return 0;
}
