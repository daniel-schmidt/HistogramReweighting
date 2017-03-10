#ifndef SINGLE_RUN_H
#define SINGLE_RUN_H

#include "solver.h"
#include "observables.h"

void single_run( struct rparams * p, double const * const sfVals ) {
  double* lambdas = p->lambdas;
  double* actionVals = p->actions;
  size_t nlambda = p->nlambda;
  size_t len_total = p->naction;
  
  double fasSolution[nlambda];
  calcSolution( p, fasSolution );
  
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
    interpol_sfabs[n]  = calcObservable( interpol_lam[n], sfabs,      p, fasSolution );
    interpol_square[n] = calcObservable( interpol_lam[n], square,     p, fasSolution );
    interpol_fourth[n] = calcObservable( interpol_lam[n], fourth,     p, fasSolution );
    interpol_Sb[n]     = calcObservable( interpol_lam[n], actionVals, p, fasSolution );
    interpol_absSb[n]  = calcObservable( interpol_lam[n], abs_Sb,     p, fasSolution );
    printf("%.16f %.16f %.16f %.16f %.16f\n"
      , interpol_lam[n]
      , interpol_sfabs[n]
      , interpol_square[n] - interpol_sfabs[n] * interpol_sfabs[n]
      , 1.-interpol_fourth[n] / (3 * interpol_square[n]*interpol_square[n] )
      , interpol_absSb[n] / interpol_sfabs[n] - interpol_Sb[n]
    );
  }
  free( sfabs );
  free( square );
  free( fourth );
  free( abs_Sb );
}
#endif