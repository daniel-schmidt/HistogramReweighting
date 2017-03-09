#include "solver.h"

void print_state (size_t iter, gsl_multiroot_fsolver * s)
{
  printf ("iter = %3zu x = % .3f, % .3f "
  "f(x) = % .3e, % .3e \n",
          iter,
          gsl_vector_get (s->x, 0), 
          gsl_vector_get (s->x, 1),
          gsl_vector_get (s->f, 0),
          gsl_vector_get (s->f, 1)
  );
}

double P( double lambda, int b, int i, void * params, double* fas ) {  
  double* lambdas = ( ( struct rparams* ) params )->lambdas;
  double** actions = ( ( struct rparams* ) params )->actions;
  int* lengths = ( ( struct rparams* ) params )->lengths;
  int nlambda = ( ( struct rparams * ) params )->nlambda;
  
  
  double denom = 0.;
  for( int a = 0; a < nlambda; ++a ) {
    denom += lengths[a] * exp(actions[b][i] * (lambda - lambdas[a]) + fas[a]);
  }
  return 1./denom;
}

int equation( const gsl_vector * x, void * params, gsl_vector *eqn ) {
  double* lambdas = ( ( struct rparams* ) params )->lambdas;
  int nlambda = ( ( struct rparams * ) params )->nlambda;
  int* lengths = ( ( struct rparams* ) params )->lengths;
  double fas[nlambda];
  double eqns[nlambda-1];
  fas[0] = 0.;
  for( int a = 1; a < nlambda; ++a ) {
    fas[a] = gsl_vector_get( x, a-1 );
  }
  
  for( int c = 1; c < nlambda; ++c ) {
    double sum = 0.;
    for( int b = 0; b < nlambda; ++b ) {
      for( int i = 0; i < lengths[b]; ++i ) {        
        sum += P( lambdas[c], b, i, params, fas );
      }
    }
    eqns[c-1] = fas[c] + log(sum);
  }
  
  for( int a = 0; a < nlambda-1; ++a ) {
    gsl_vector_set( eqn, a, eqns[a] );
  }
  
  return GSL_SUCCESS;
}

void calcSolution( struct rparams* params, double* sol ) {
  int nlambda = params->nlambda;
  
  // setting initial values
  gsl_vector *fa = gsl_vector_alloc( nlambda-1 );
  for( int numLambda = 0; numLambda < nlambda-1; ++numLambda ) {
    gsl_vector_set( fa, numLambda, 10*(numLambda+1) );
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
  
  sol[0] = 0.;
  for( int a = 1; a < nlambda; ++a ) {
    sol[a] = gsl_vector_get( s->x, a-1 );
  }
  
  gsl_multiroot_fsolver_free (s);
  gsl_vector_free( fa );
  gsl_vector_free( eqns );
}

double calcObservable( double lambda, double** observableData, void* params, double* fasSolution ){
  int nlambda = ( ( struct rparams * ) params )->nlambda;
  int* lengths = ( ( struct rparams* ) params )->lengths;
  
  double denom = 0.;
  for( int b = 0; b < nlambda; ++b ) {
    for( int i = 0; i < lengths[b]; ++i ) {        
      denom += P( lambda, b, i, params, fasSolution );
    }
  }
  
  double numerator = 0.;
  for( int b = 0; b < nlambda; ++b ) {
    for( int i = 0; i < lengths[b]; ++i ) {        
      numerator += observableData[b][i] * P( lambda, b, i, params, fasSolution );
    }
  }
  
  return numerator / denom;
}