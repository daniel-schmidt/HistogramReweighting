#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

#include "io.h"
#include "solver.h"

int main( void ) {
  
  // Load data
  int nlambda = 3;
  double lambdas[nlambda];
  char* sfNames[nlambda];
  char* actionNames[nlambda];
  
  lambdas[0] = 0.46;
  lambdas[1] = 0.47;
  lambdas[2] = 0.48;
    
  sfNames[0] = "/data2/Results/GN/red/24x23x23/results_1/Configs/ScalarField_ScalarOnConfig_.46.dat";
  sfNames[1] = "/data2/Results/GN/red/24x23x23/results_1/Configs/ScalarField_ScalarOnConfig_.47.dat";
  sfNames[2] = "/data2/Results/GN/red/24x23x23/results_1/Configs/ScalarField_ScalarOnConfig_.48.dat";
  actionNames[0] = "/data2/Results/GN/red/24x23x23/results_1/Configs/BosonicAction_ScalarOnConfig_.46.dat";
  actionNames[1] = "/data2/Results/GN/red/24x23x23/results_1/Configs/BosonicAction_ScalarOnConfig_.47.dat";
  actionNames[2] = "/data2/Results/GN/red/24x23x23/results_1/Configs/BosonicAction_ScalarOnConfig_.48.dat";
  
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