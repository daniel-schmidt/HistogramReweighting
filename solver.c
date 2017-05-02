#include "solver.h"
gsl_vector* fa = NULL;

void print_state (size_t iter, gsl_multiroot_fsolver * s)
{
  printf ("iter = %2zu x = % .6f, % .6f "
  "f(x) = % .6e, % .6e \n",
          iter,
          gsl_vector_get (s->x, 0), 
          gsl_vector_get (s->x, 1),
          gsl_vector_get (s->f, 0),
          gsl_vector_get (s->f, 1)
  );
}

long double P( double lambda, int bi, void * params, double const * const fas ) {  
  double* lambdas = ( ( struct rparams* ) params )->lambdas;
  double* g = ( ( struct rparams* ) params )->autocorr;
  double* actions = ( ( struct rparams* ) params )->actions;
  int* lengths = ( ( struct rparams* ) params )->lengths;
  int nlambda = ( ( struct rparams * ) params )->nlambda;
  
  long double denom = 0.L;
  for( int a = 0; a < nlambda; ++a ) {
    denom += lengths[a] * g[a] * expl( (long double) (actions[bi] * (lambda - lambdas[a]) + fas[a]));
  }
  
  // determine b for g[b] from bi by subtracting previous lengths of data
  int b;
  for( b = 0; b < nlambda; ++b ) {
    bi -= lengths[b];
    if( bi < 0 ) {
      break;
    }
  }
  return g[b]/denom;
}

int equation( const gsl_vector * x, void * params, gsl_vector *eqn ) {
  double* lambdas = ( ( struct rparams* ) params )->lambdas;
  int nlambda = ( ( struct rparams * ) params )->nlambda;
  size_t naction = ( ( struct rparams* ) params )->naction;
  
  double fas[nlambda];
  double eqns[nlambda-1];
  fas[0] = ( ( struct rparams* ) params )->f0;
  for( int a = 1; a < nlambda; ++a ) {
    fas[a] = gsl_vector_get( x, a-1 );
  }
  
  for( int c = 1; c < nlambda; ++c ) {
    long double sum = 0.L;
    for( int bi = 0; bi < naction; ++bi ) {
      sum += P( lambdas[c], bi, params, fas );
    }
    eqns[c-1] = fas[c] + (double) logl(sum);
  }
  
  for( int a = 0; a < nlambda-1; ++a ) {
    gsl_vector_set( eqn, a, eqns[a] );
  }
  
  return GSL_SUCCESS;
}

void calcSolution( struct rparams * params, double* sol ) {
  int nlambda = params->nlambda;
  
  // setting initial values, differences of 10 seem to work quite general
  if( fa == NULL ) {
    const double del_fa = 1.;
    fa = gsl_vector_alloc( nlambda-1 );
    for( int numLambda = 0; numLambda < nlambda-1; ++numLambda ) {
      gsl_vector_set( fa, numLambda, del_fa*(numLambda+1) + (params->f0) );
    }
  }
  
  // testing equation evaluation on initial values
  gsl_vector *eqns = gsl_vector_alloc( nlambda-1 );
  equation( fa, params, eqns );
    
  const size_t numEqns = nlambda - 1;
  size_t iter = 0;
  int status;
  
  const gsl_multiroot_fsolver_type *T;
  gsl_multiroot_fsolver *s;
    
  gsl_multiroot_function f = {&equation, numEqns, params};
  T = gsl_multiroot_fsolver_hybrids;
  s = gsl_multiroot_fsolver_alloc(T, numEqns);
  gsl_multiroot_fsolver_set( s, &f, fa );
  
  print_state(iter, s);
  
  do
  {
    iter++;
    status = gsl_multiroot_fsolver_iterate (s);
    
    print_state (iter, s);
    
    if (status)   /* check if solver is stuck */
      break;
    
    status = gsl_multiroot_test_residual (s->f, 1e-7);
  }
  while (status == GSL_CONTINUE && iter < 1000);
  
  printf ("status = %s\n", gsl_strerror (status));
  
  sol[0] = params->f0;
  for( int a = 1; a < nlambda; ++a ) {
    sol[a] = gsl_vector_get( s->x, a-1 );
    gsl_vector_set( fa, a-1, sol[a] );
  }
  gsl_multiroot_fsolver_free (s);
  gsl_vector_free( eqns );
}

void freeSolver() {
  gsl_vector_free( fa );
}
  