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

long double calcPTable( double lambda, void* params, double const * const fasSolution, double* const PTable ) {
  size_t naction = ( ( struct rparams* ) params )->naction;
  long double denom = 0.L;
  for( size_t bi = 0; bi < naction; ++bi ) {
     PTable[bi] = P( lambda, bi, params, fasSolution );
     denom += PTable[bi];
  }
  return denom;
}

double calcObservable( double const * const observableData, const long double denom, double const * const PTable, const size_t naction ) {
  long double numerator = 0.L;
  for( size_t bi = 0; bi < naction; ++bi ) {
    numerator += observableData[bi] * PTable[bi];
  }
  //   printf( "numerator: %.10Le, denom: %.10Le, quotient: %.10Le\n", numerator, denom, numerator/denom);
  return (double) (numerator / denom);
}