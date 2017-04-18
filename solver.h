#ifndef SOLVER_H
#define SOLVER_H

#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

struct rparams {
  double* lambdas;
  double* autocorr;
  double* actions;
  int* lengths;
  int nlambda;
  size_t naction;
  double f0;
};

void print_state (size_t iter, gsl_multiroot_fsolver * s);

long double P( double lambda, int bi, void * params, double const * const fas );

int equation( const gsl_vector * x, void * params, gsl_vector *eqn );

void calcSolution( struct rparams * params, double* sol );

void freeSolver();

#endif