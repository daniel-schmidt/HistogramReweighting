#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

#include "io.h"
#include "solver.h"

int main( void ) {
  
  // Load data
  int nlambda = 4;
  double lambdas[nlambda];
  char* sfNames[nlambda];
  char* actionNames[nlambda];
  
  lambdas[0] = 0.46;
  lambdas[1] = 0.47;
  lambdas[2] = 0.48;
  lambdas[3] = 0.49;
  
  readPathsFromFile( "sf_paths.txt", nlambda, sfNames );
  readPathsFromFile( "action_paths.txt", nlambda, actionNames );
  
  double* sfVals[nlambda];
  double* actionVals[nlambda];
  int length[nlambda];

  readData( nlambda, sfNames, sfVals, actionNames, actionVals, length );
  // Set parameters and calculate solution
  struct rparams p = {
    lambdas,
    actionVals,
    length,
    nlambda,
  };
  
  double fasSolution[nlambda];
  calcSolution( &p, fasSolution );
  
  // Use solution to calculate interpolations
  double result = calcObservable( 0.44, sfVals, &p, fasSolution );
  printf("value for 0.44: %.6f\n", result);
  
  // Cleanup
  for( int numLambda = 0; numLambda < nlambda-1; ++numLambda ){
    free( sfVals[numLambda] );
    free( actionVals[numLambda] );
  }

  return EXIT_SUCCESS;
}