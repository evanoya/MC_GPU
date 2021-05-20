#include <stdio.h>
#include <curand.h>

/*
  Add in this handle any code you consider appropriate to manage CUDA error
  returns
*/
static void ManejaError( cudaError_t error, const char *file, int line ) {
    if (error != cudaSuccess) {
        printf( "Cuda Error exception !! :: %s in %s , line %d\n", cudaGetErrorString( error ), file, line );
        exit( EXIT_FAILURE );
    }
}
#define MANEJA_ERROR( error ) (ManejaError( error, __FILE__, __LINE__ ))

    const char* curandGetErrorString(curandStatus_t status)
    {
        switch(status)
        {
        case CURAND_STATUS_SUCCESS: return "CURAND_STATUS_SUCCESS";
    	case CURAND_STATUS_VERSION_MISMATCH: return "CURAND_STATUS_VERSION_MISMATCH";
    	case CURAND_STATUS_NOT_INITIALIZED: return "CURAND_STATUS_NOT_INITIALIZED";
    	case CURAND_STATUS_ALLOCATION_FAILED: return "CURAND_STATUS_ALLOCATION_FAILED";
    	case CURAND_STATUS_TYPE_ERROR: return "CURAND_STATUS_TYPE_ERROR";
    	case CURAND_STATUS_OUT_OF_RANGE: return "CURAND_STATUS_OUT_OF_RANGE";
    	case CURAND_STATUS_LENGTH_NOT_MULTIPLE: return "CURAND_STATUS_LENGTH_NOT_MULTIPLE";
    	case CURAND_STATUS_DOUBLE_PRECISION_REQUIRED: return "CURAND_STATUS_DOUBLE_PRECISION_REQUIRED";
    	case CURAND_STATUS_LAUNCH_FAILURE: return "CURAND_STATUS_LAUNCH_FAILURE";
    	case CURAND_STATUS_PREEXISTING_FAILURE: return "CURAND_STATUS_PREEXISTING_FAILURE";
    	case CURAND_STATUS_INITIALIZATION_FAILED: return "CURAND_STATUS_INITIALIZATION_FAILED";
    	case CURAND_STATUS_ARCH_MISMATCH: return "CURAND_STATUS_ARCH_MISMATCH";
    	case CURAND_STATUS_INTERNAL_ERROR: return "CURAND_STATUS_INTERNAL_ERROR";
        }
        return "Unknown cuRAND error";
    }

static void ManejaRndError( curandStatus_t error, const char *file, int line ) {
    if (error != CURAND_STATUS_SUCCESS) {
        printf( "Curand Error exception !! :: %s in %s , line %d\n", curandGetErrorString(error) , file, line );
        exit( EXIT_FAILURE );
    }
}
#define MANEJA_RND_ERROR( error ) (ManejaRndError( error, __FILE__, __LINE__ ))


    
