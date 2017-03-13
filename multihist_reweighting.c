#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

#include "io.h"
#include "single_run.h"

int main( int argc, char** argv ) {
  if( argc != 4 ) {
    printf( "ERROR: Need 3 input parameters: lambdas.txt sf_paths.txt action_paths.txt\n" );
    exit(1);
  }
  
  if( sizeof(double) >= sizeof(long double) ){
     printf("WARNING: long double seems no longer than double: %lu, long double: %lu", sizeof(double), sizeof(long double));
  }

  double* lambdas = NULL;
  size_t nlambda = readLambdasFromFile( argv[1], &lambdas );
  char* sfNames[nlambda];
  char* actionNames[nlambda];
  
  for( size_t n = 0; n < nlambda; ++n ) {
    printf("lamb=%.3f\n", lambdas[n]);
  }
  
  readPathsFromFile( argv[2], nlambda, sfNames );
  readPathsFromFile( argv[3], nlambda, actionNames );
  
  double* sfVals = NULL;
  double* actionVals = NULL;
  int lengths[nlambda];

  size_t numThermal = 100;
  size_t len_total = readData( numThermal, nlambda, sfNames, &sfVals, actionNames, &actionVals, lengths );
  printf("Read a total of %zu data points.\n", len_total);

  // Set parameters and calculate solution
  struct rparams p = {
    lambdas,
    actionVals,
    lengths,
    nlambda,
    len_total
  };
  
  single_run( &p, sfVals );
  
  // binning and bootstrapping for error estimates
//   int seed = 12;
//   srand(seed);
//   size_t bin_size = 100;
//   double* actionSelect = malloc( len_total * sizeof *actionVals );
//   double* sfSelect = malloc( len_total * sizeof *sfVals );
//   random_select( actionVals, sfVals, lengths, nlambda, bin_size, actionSelect, sfSelect );
  
//   for( size_t i = 0; i < len_total; ++i ) {
//     printf("%.6f %.6f\n",actionSelect[i], sfSelect[i]);
//   }
  
  // Cleanup
  free( sfVals );
  free( actionVals );

  return EXIT_SUCCESS;
}