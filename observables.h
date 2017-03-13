#ifndef OBSERVABLES_H
#define OBSERVABLES_H

#include <math.h>
#include "solver.h"

void calculateOnConfigData( double const* const sf, double const* const action, const size_t len_total, double* sfabs, double* square, double* fourth_power, double* abs_times_action );

long double calcPTable( const double lambda
                      , void* params
                      , double const * const fasSolution
                      , double * const PTable
                      );

double calcObservable( double const * const observableData
                     , const long double denom
                     , double const * const PTable
                     , const size_t naction
                     );
#endif