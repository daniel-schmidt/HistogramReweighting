#ifndef SOLVER_H
#define SOLVER_H

#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

struct rparams {
  double* lambdas;
  double** actions;
  int* lengths;
  int nlambda;
};

void print_state (size_t iter, gsl_multiroot_fsolver * s);

long double P( double lambda, int b, int i, void * params, double* fas );

int equation( const gsl_vector * x, void * params, gsl_vector *eqn );

void calcSolution( struct rparams* params, double* sol );

double calcObservable( double lambda, double** observableData, void* params, double* fasSolution );

#endif