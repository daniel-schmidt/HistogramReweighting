#ifndef SINGLE_RUN_H
#define SINGLE_RUN_H

#include <string.h>
#include "solver.h"
#include "observables.h"

size_t randint( size_t max );

void random_select( double const * const actionVals
                  , double const * const sfVals
                  , int* lengths
                  , size_t nlambda
                  , size_t bin_size
                  , double* actionSelect
                  , double* sfSelect
                  );

void single_run( struct rparams * p, double const * const sfVals );

#endif