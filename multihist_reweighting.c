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
  size_t numInterpol = 101;
  
  struct rparams p = {
    lambdas,
    actionVals,
    lengths,
    nlambda,
    len_total
  };
  
  double ip_lam   [numInterpol];
  double ip_sfabs [numInterpol];
  double ip_sus   [numInterpol];
  double ip_bc    [numInterpol];
  double ip_dlog  [numInterpol];
  
  single_run( &p, sfVals, numInterpol, ip_lam, ip_sfabs, ip_sus, ip_bc, ip_dlog );
  FILE * file = fopen("ip.dat", "w");
  for( size_t ip = 0; ip < numInterpol; ++ip ) {
    fprintf( file, "%.16f %.16f %.16f %.16f %.16f\n", ip_lam[ip], ip_sfabs[ip], ip_sus[ip], ip_bc[ip], ip_dlog[ip]);
  }
  
  // binning and bootstrapping for error estimates
  int seed = 12;
  srand(seed);
  size_t bin_size = 100;
  size_t Nboot = 2;
  
  double* actionSelect = malloc( len_total * sizeof *actionVals );
  double* sfSelect = malloc( len_total * sizeof *sfVals );
  
  p.actions = actionSelect;
  
  for( size_t boot = 0; boot < Nboot; ++boot ) {
//     random_select( actionVals, sfVals, lengths, nlambda, bin_size, actionSelect, sfSelect );
//     single_run( &p, sfSelect );
  }
  
  // Cleanup
  free( sfVals );
  free( actionVals );

  return EXIT_SUCCESS;
}