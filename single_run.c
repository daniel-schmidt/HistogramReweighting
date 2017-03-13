#include "single_run.h"

/* Returns an integer in the range [0, n).
 * Uses rand(), and so is affected-by/affects the same seed.
 * see: http://stackoverflow.com/questions/822323/how-to-generate-a-random-number-in-c
 */
size_t randint(size_t n) {
  if ((n - 1) == RAND_MAX) {
    return rand();
  } else {
    // Chop off all of the values that would cause skew...
    long end = RAND_MAX / n; // truncate skew
    end *= n;

    // ... and ignore results from rand() that fall above that limit.
    size_t r;
    while ((r = rand()) >= end);

    return r % n;
  }
}

void random_select( double const * const actionVals, double const * const sfVals, int* lengths, size_t nlambda,        size_t bin_size, double* actionSelect, double* sfSelect) {
  
  size_t offset = 0;
  
  for( size_t a = 0; a < nlambda; ++a ) {
    size_t num_bins = lengths[a]/bin_size;
    if( lengths[a] % bin_size != 0 ) {
      printf("WARNING: in random_select: bin size %zu is not a divider of data length %d.\n", bin_size, lengths[a]);
      num_bins++;
    }
    
    for( size_t b = 0; b < num_bins; ++b ) {
      size_t bin_idx = randint( num_bins );
      printf( "Chose block %zu.\n", bin_idx );
      
      memcpy( actionSelect + offset + b * bin_size
            , actionVals + offset + bin_idx * bin_size
            , bin_size * sizeof *actionVals
            );
      
      memcpy( sfSelect + offset + b * bin_size
            , sfVals + offset + bin_idx * bin_size
            , bin_size * sizeof *sfVals
            );
      
    }
    offset += lengths[a];
  }
}


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