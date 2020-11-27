/*
* FILE: isingV2.cu
* THMMY, 7th semester, Parallel and Distributed Systems: 3rd assignment
* Parallel Implementation  of the Ising Model
* Authors:
*   Moustaklis Apostolos, 9127, amoustakl@ece.auth.gr
*   Papadakis Charis , 9128, papadakic@ece.auth.gr
* Compile command with :
*   make all
* Run command example:
*   ./src/isingV2
* It will calculate the evolution of the ising Model
* for a given number n  of points and k steps
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

// Defines for the block and grid calculation
#define BLOCK_SIZE_X 16
#define BLOCK_SIZE_Y 16
//The dimensions of the lattice
#define N_X 517
#define N_Y 517

// The size of the weights
#define WSIZE 5

//Helper Defines to access easier the arrays
#define old(i,j,n) *(old+i*n+j)
#define current(i,j,n) *(current+i*n+j)
#define w(i,j) *(w+i*5+j)
#define d_w(i,j) *(d_w+i*5+j)
#define G(i,j,n) *(G+i*n+j)
#define d_current(i,j,n) *(d_current+i*n+j)
#define d_old(i,j,n) *(d_old+i*n+j)

//Functions Declaration
void swapElement(int  ** one, int  ** two);
__global__
   void kernel2D(int *d_current, int *d_old, double *d_w, int n , int * d_flag);
void ising( int *G, double *w, int k, int n);


//! Ising model evolution
/*!

  \param G      Spins on the square lattice             [n-by-n]
  \param w      Weight matrix                           [5-by-5]
  \param k      Number of iterations                    [scalar]
  \param n      Number of lattice points per dim        [scalar]

  NOTE: Both matrices G and w are stored in row-major format.
*/

void ising( int *G, double *w, int k, int n){

  //Grid and block construction
  dim3 block(BLOCK_SIZE_X,BLOCK_SIZE_Y);
  int grid_size_x  = (N_X + BLOCK_SIZE_X - 1) / BLOCK_SIZE_X;
  int grid_size_y  = (N_Y + BLOCK_SIZE_Y - 1) / BLOCK_SIZE_Y;
  dim3 grid(grid_size_x,grid_size_y);
  //Device memory allocation
  int * old = (int*) malloc(n*n*(size_t)sizeof(int)); // old spin lattice
  int * current = (int*) malloc(n*n*(size_t)sizeof(int)); // current spin lattice
  //Leak check
  if(old==NULL || current == NULL){
      printf("Problem at memory allocation at host \n");
        exit(0);
      }

  int * d_old;
  int * d_current;
  double * d_w;// na valw void ** skatakia
  int *d_flag ;
  int flag ;
  //Host memory allocation and leak check
  if( cudaMalloc((void **)&d_old ,n*n*(size_t)sizeof(int)) != cudaSuccess  || cudaMalloc((void **)&d_current,n*n*(size_t)sizeof(int))   != cudaSuccess   || cudaMalloc((void **)&d_w, WSIZE*WSIZE*(size_t)sizeof(double))   != cudaSuccess || cudaMalloc(&d_flag,(size_t)sizeof(int)) !=cudaSuccess){
    printf("Problem at memory allocation");
    exit(0);
  }
  //Copy memory from host to device
  cudaMemcpy(d_w, w, WSIZE*WSIZE*sizeof(double), cudaMemcpyHostToDevice );
  cudaMemcpy(d_old, G, n*n*sizeof(int), cudaMemcpyHostToDevice );


  // run for k steps
  for(int l=0; l<k; l++){

    flag = 0;

    kernel2D<<<grid,block>>>(d_current, d_old, d_w, n  , d_flag );
    cudaDeviceSynchronize();

  //cudaMemcpy(old, d_old, n*n*sizeof(int), cudaMemcpyDeviceToHost );
    cudaMemcpy(current, d_current, n*n*sizeof(int), cudaMemcpyDeviceToHost );
    // save result in G
    memcpy(G , current , n*n*sizeof(int));

    // swap the pointers for the next iteration
    swapElement(&d_old,&d_current);

    cudaMemcpy(&flag , d_flag , (size_t)sizeof(int), cudaMemcpyDeviceToHost);
    // terminate if no changes are made
    if(flag){
        printf("terminated: spin values stay same (step %d)\n" , l );
      exit(0);
    }
  }
  //Memory deallocation
  free(old);
  free(current);
  cudaFree(d_old);
  cudaFree(d_current);
  cudaFree(d_w);
}

//Helper function to swap the pointers of the arrays
void swapElement(int  ** one, int  ** two) {
  int  * temp = * one;
  * one = * two;
  * two = temp;
}

 //The kernel function that updates the values of the ising model
__global__
void kernel2D(int *d_current, int *d_old, double *d_w, int n , int * d_flag)
{

  double influence = 0;
  // Compute global  column and row indices.
  int r = blockIdx.x * blockDim.x + threadIdx.x;
  int c = blockIdx.y * blockDim.y + threadIdx.y;

  // Check if within bounds.
  if ((c >= n) || (r >= n))
  return;

  for(int i = r; i<n; i+=blockDim.x*gridDim.x){
    for(int j = c; j<n; j+=blockDim.y*gridDim.y){

      for(int ii=0; ii<5; ii++){
        for(int jj=0; jj<5; jj++){
          influence +=  d_w(ii,jj) * d_old((i-2+n+ii)%n, (j-2+n+jj)%n, n);
        }
      }
      // magnetic moment gets the value of the SIGN of the weighted influence of its neighbors
      if(fabs(influence) < 10e-7){
        d_current(i,j,n) = d_old(i,j,n); // remains the same in the case that the weighted influence is zero
      }
      else if(influence > 10e-7){
        d_current(i,j,n) = 1;
        *d_flag = 0;
      }
      else if(influence < 0){
        d_current(i,j,n) = -1;
        *d_flag=0;
      }
      influence = 0;
    }
  }
}
