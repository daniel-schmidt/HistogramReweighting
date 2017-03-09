#ifndef OBSERVABLES_H
#define OBSERVABLES_H

#include <math.h>
#include "solver.h"

void calculateOnConfigData( double const* const* sf, const size_t nlambda, const size_t* length, double** sfabs, double** square, double** fourth_power );
double calcObservable( double lambda, double* observableData, void* params, double* fasSolution );
#endif