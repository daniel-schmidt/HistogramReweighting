#ifndef OBSERVABLES_H
#define OBSERVABLES_H

#include <math.h>
#include "solver.h"

void calculateOnConfigData( double const* const sf, double const* const action, const size_t len_total, double* sfabs, double* square, double* fourth_power, double* abs_times_action );

double calcObservable( double lambda, double* observableData, void* params, double* fasSolution );
#endif