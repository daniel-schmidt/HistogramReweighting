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


void single_run( struct rparams * p, double const * const sfVals, size_t const numInterpol, double* const ip_lam, double* const ip_sfabs, double* const ip_sus, double* const ip_bc, double* const ip_dlog ) {
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
  double lam_min = lambdas[0] - 0.02;
  double lam_max = lambdas[nlambda-1] + 0.02;
  double d_lam = (lam_max - lam_min) / (numInterpol-1);
  printf( "Calculating interpolation from %.3f to %.3f in steps of %.3f\n", lam_min, lam_max, d_lam);
  
  for( size_t n = 0; n < numInterpol; ++n )
  {
    ip_lam[n]    = lam_min + n * d_lam;
    ip_sfabs[n]  = calcObservable( ip_lam[n], sfabs,      p, fasSolution );
    double interpol_square = calcObservable( ip_lam[n], square,     p, fasSolution );
    double interpol_fourth = calcObservable( ip_lam[n], fourth,     p, fasSolution );
    double interpol_Sb     = calcObservable( ip_lam[n], actionVals, p, fasSolution );
    double interpol_absSb  = calcObservable( ip_lam[n], abs_Sb,     p, fasSolution );
    
    ip_sus[n] = interpol_square - ip_sfabs[n] * ip_sfabs[n];
    ip_bc[n] = 1.-interpol_fourth / (3 * interpol_square * interpol_square );
    ip_dlog[n] = interpol_absSb / ip_sfabs[n] - interpol_Sb;
  }
  free( sfabs );
  free( square );
  free( fourth );
  free( abs_Sb );
}