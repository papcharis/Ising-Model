/*
* FILE: isingV0.c
* THMMY, 7th semester, Parallel and Distributed Systems: 3rd assignment
* Sequential Implementation with shared memory of the Ising Model
* Authors:
*   Moustaklis Apostolos, 9127, amoustakl@ece.auth.gr
*   Papadakis Charis , 9128, papadakic@ece.auth.gr
* Compile command with :
*   make all
* Run command example:
*   ./src/isingV0
* It will calculate the evolution of the ising Model
* for a given number n  of points and k steps
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

//Helper Defines to access easier the arrays
#define old(i,j,n) *(old+(i)*n+j)
#define current(i,j,n) *(current+(i)*n+j)
#define w(i,j) *(w+(i)*5+j)
#define G(i,j,n) *(G+(i)*n+j)


//Functions Declaration
void swapElement(int  ** one, int  ** two);
double influenceCalc(int *old , double *w , int n , int i , int j  );
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

  int * old = (int*) malloc(n*n*sizeof(int)); // old spin lattice
  int * current= (int*) malloc(n*n*sizeof(int)); // current spin lattice
//  int * tmp;
  double influence; // weighted influence of the neighbors
  double flag = 1;

//Elearning tester checks the values of the G so by swaping
// The "head" pointer it can not pass the validation
// So we manual copy
  memcpy(old,G,n*n*sizeof(int));

  // run for k steps
  for(int l=0; l<k; l++){

    // for each current[i][j] point
    for(int i=0; i<n; i++){
      for(int j=0; j<n; j++){

        // calculation of weighted influence
        influence = influenceCalc(old , w ,n ,i , j);

        // magnetic moment gets the value of the SIGN of the weighted influence of its neighbors
        if(fabs(influence) < 10e-7){
          current(i,j,n) = old(i,j,n); // remains the same in the case that the weighted influence is zero
        }
        else if(influence > 10e-7){
          current(i,j,n) = 1;
          flag = 0;
        }
        else if(influence < 0){
          current(i,j,n) = -1;
          flag = 0;
        }
        // save result in G
        G(i,j,n) = current(i,j,n);
      }
    }

    // swap the pointers for the next iteration
    swapElement(&old,&current);
    // terminate if no changes are made
    if(flag ==1){
       printf("terminated: spin values stay same (step %d)\n" , l);
    }
  }
}

void swapElement(int  ** one, int  ** two) {
  int  * temp = * one;
  * one = * two;
  * two = temp;
}

double influenceCalc(int *old , double *w , int n , int i , int j  ){
  double influence = 0;
  for(int ii=0; ii<5; ii++){
    for(int jj=0; jj<5; jj++){
      influence +=  w(ii,jj) * old((i-2+n+ii)%n, (j-2+n+jj)%n, n);
    }
  }
  return influence;
}

int main(int argc, const char* argv[]){


  int n = 517;
  int fr;
  double weights[] = {0.004, 0.016, 0.026, 0.016, 0.004, 0.016, 0.071, 0.117, 0.071, 0.016, 0.026, 0.117, 0.0, 0.117, 0.026, 0.016, 0.071, 0.117, 0.071, 0.016, 0.004, 0.016, 0.026, 0.016, 0.004};
  int G1[n*n]; // G that changes k times
  int G2[n*n]; // G that changes k times
  int G3[n*n]; // G that changes k times
  int Gk1[n*n];
  int Gk2[n*n];
  int Gk3[n*n];
  FILE *ptr;

  // read initial G
  ptr = fopen("./conf-init.bin","rb");
  fr = fread(G1,sizeof(G1),1,ptr);
  fclose(ptr);
  // read initial G
  ptr = fopen("./conf-init.bin","rb");
  fr = fread(G2,sizeof(G2),1,ptr);
  fclose(ptr);
  // read initial G
  ptr = fopen("./conf-init.bin","rb");
  fr = fread(G3,sizeof(G3),1,ptr);
  fclose(ptr);
  // read k-th Gk
  ptr = fopen("./conf-1.bin","rb"); // allazo onoma arxeiou gia allagi k
  fr = fread(Gk1,sizeof(Gk1),1,ptr);
  fclose(ptr);

  ptr = fopen("./conf-4.bin","rb"); // allazo onoma arxeiou gia allagi k
  fr = fread(Gk2,sizeof(Gk2),1,ptr);
  fclose(ptr);

  ptr = fopen("./conf-11.bin","rb"); // allazo onoma arxeiou gia allagi k
  fr = fread(Gk3,sizeof(Gk3),1,ptr);
  fclose(ptr);

  // execution
  ising(G1, weights, 1, n);

  // check correctness
  int c = 0;
  for(int i=0; i<n; i++){
    for(int j=0; j<n; j++){
      if( *(G1+i*n+j) != *(Gk1+i*n+j) ){
        c++;
      }
    }
  }

  if(c!=0){
    printf("k=1 Wrong\n");
  }
  else{
    printf("k=1 Correct\n");
  }

  // execution
  ising(G2, weights, 4, n);

  // check correctness
   c = 0;
  for(int i=0; i<n; i++){
    for(int j=0; j<n; j++){
      if( *(G2+i*n+j) != *(Gk2+i*n+j) ){
        c++;
      }
    }
  }
  if(c!=0){
    printf("k=4 Wrong\n");
  }
  else{
    printf("k=4 Correct\n");
  }


  // // execution
  // ising(G, weights, 100, 1000);
  // execution
  ising(G3, weights, 11 , n);



  // check correctness
   c = 0;
  for(int i=0; i<n; i++){
    for(int j=0; j<n; j++){
      if( *(G3+i*n+j) != *(Gk3+i*n+j) ){
        c++;
      }
    }
  }
  if(c!=0){
    printf("k=11 Wrong\n");
  }
  else{
    printf("k=11 Correct\n");
  }


  }
