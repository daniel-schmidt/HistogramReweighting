#include "observables.h"

void calculateOnConfigData( double const* const sf, double const* const action, const size_t len_total, double* sfabs, double* square, double* fourth_power, double* abs_times_action ) {
  for( size_t ai = 0; ai < len_total; ++ai ) {
    double valOnConf = sf[ai];
    double absOnConf = fabs( valOnConf );
    sfabs[ai] = absOnConf;
    double valSquare = valOnConf*valOnConf;
    square[ai] = valSquare;
    fourth_power[ai] = valSquare * valSquare;
    abs_times_action[ai] = absOnConf * action[ai];
  }
}

double calcObservable( double lambda, double* observableData, void* params, double* fasSolution ) {
  size_t naction = ( ( struct rparams* ) params )->naction;
  
  long double denom = 0.L;
  for( size_t bi = 0; bi < naction; ++bi ) {
    denom += P( lambda, bi, params, fasSolution );
    
  }
  
  long double numerator = 0.L;
  for( size_t bi = 0; bi < naction; ++bi ) {
    numerator += observableData[bi] * P( lambda, bi, params, fasSolution );
  }
//   printf( "numerator: %.10Le, denom: %.10Le, quotient: %.10Le\n", numerator, denom, numerator/denom);
  return (double) (numerator / denom);
}