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
  size_t numInterpol = 100;
  double lam_min = lambdas[0] - 0.02;
  double lam_max = lambdas[nlambda-1] + 0.02;
  double d_lam = (lam_max - lam_min) / numInterpol;
  printf( "Calculating interpolation from %.3f to %.3f in steps of %.3f\n", lam_min, lam_max, d_lam);
  
  double interpol_x[numInterpol];
  double interpol_y[numInterpol];
  for( size_t n = 0; n <= numInterpol; ++n )
  {
    interpol_x[n] = lam_min + n * d_lam;
    interpol_y[n] = calcObservable( interpol_x[n], sfVals, &p, fasSolution );
    printf("%.6f %.6f\n", interpol_x[n], interpol_y[n]);
  }
  
  // Cleanup
  for( int numLambda = 0; numLambda < nlambda-1; ++numLambda ){
    free( sfVals[numLambda] );
    free( actionVals[numLambda] );
  }

  return EXIT_SUCCESS;
}