#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

#include "io.h"
#include "solver.h"
#include "observables.h"

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
  
  double* sfVals = NULL;
  double* actionVals = NULL;
  int lengths[nlambda];

  size_t len_total = readData( nlambda, sfNames, &sfVals, actionNames, &actionVals, lengths );
  printf("Read a total of %zu data points to %p and %p\n", len_total, sfVals, actionVals);
//   for( size_t i = 0; i < len_total; ++i ) {
//     printf("sf: %.2f, act: %.2f\n", sfVals[i], actionVals[i]);
//   }
  // Set parameters and calculate solution
  struct rparams p = {
    lambdas,
    actionVals,
    lengths,
    nlambda,
    len_total
  };
  
  double fasSolution[nlambda];
  calcSolution( &p, fasSolution );
  
  // Calculate input data for observables
  double* sfabs  = malloc( len_total * sizeof *sfVals );
  double* square = malloc( len_total * sizeof *sfVals );
  double* fourth = malloc( len_total * sizeof *sfVals );
  double* abs_Sb = malloc( len_total * sizeof *sfVals );
  calculateOnConfigData( sfVals, actionVals, len_total, sfabs, square, fourth, abs_Sb ); 
  
  // Use solution to calculate interpolations
  size_t numInterpol = 100;
  double lam_min = lambdas[0] - 0.02;
  double lam_max = lambdas[nlambda-1] + 0.02;
  double d_lam = (lam_max - lam_min) / numInterpol;
  printf( "Calculating interpolation from %.3f to %.3f in steps of %.3f\n", lam_min, lam_max, d_lam);
  
  double interpol_lam    [numInterpol];
  double interpol_sfabs  [numInterpol];
  double interpol_square [numInterpol];
  double interpol_fourth [numInterpol];
  double interpol_Sb     [numInterpol];
  double interpol_absSb  [numInterpol];
  for( size_t n = 0; n <= numInterpol; ++n )
  {
    interpol_lam[n]    = lam_min + n * d_lam;
    interpol_sfabs[n]  = calcObservable( interpol_lam[n], sfabs,      &p, fasSolution );
    interpol_square[n] = calcObservable( interpol_lam[n], square,     &p, fasSolution );
    interpol_fourth[n] = calcObservable( interpol_lam[n], fourth,     &p, fasSolution );
    interpol_Sb[n]     = calcObservable( interpol_lam[n], actionVals, &p, fasSolution );
    interpol_absSb[n]  = calcObservable( interpol_lam[n], abs_Sb,     &p, fasSolution );
    printf("%.6f %.6f %.6f %.6f %.6f\n"
      , interpol_lam[n]
      , interpol_sfabs[n]
      , interpol_square[n] - interpol_sfabs[n] * interpol_sfabs[n]
      , 1.-interpol_fourth[n] / (3 * interpol_square[n]*interpol_square[n] )
      , interpol_absSb[n] / interpol_sfabs[n] - interpol_Sb[n]
    );
  }
  
  // Cleanup
  free( sfVals );
  free( actionVals );
  free( sfabs );
  free( square );
  free( fourth );
  return EXIT_SUCCESS;
}