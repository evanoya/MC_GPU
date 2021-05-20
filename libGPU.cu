#include "errhand.h"
#include <stdio.h>
#include <math.h>
#define NBLOCKMAX  65535



extern "C" void gpucount_(int *nogpus) {
  int ndev;
  // get the number of GPUs
  cudaError_t cudaRescode = cudaGetDeviceCount( &ndev ) ;
if (cudaRescode != 0){
  *nogpus = 0;
 } else {
  *nogpus = ndev;
 }
}

