#include "observables.h"

void calculateOnConfigData( double const* const* sf, const size_t nlambda, const size_t* lengths, double** sfabs, double** square, double** fourth_power ) {
  for( size_t a = 0; a < nlambda; ++a ) {
    for( size_t i = 0; i < lengths[a]; ++i ) {
      double valOnConf = sf[a][i];
      sfabs[a][i] = fabs( valOnConf );
      double valSquare = valOnConf*valOnConf;
      square[a][i] = valSquare;
      fourth_power[a][i] = valSquare * valSquare;
    }
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