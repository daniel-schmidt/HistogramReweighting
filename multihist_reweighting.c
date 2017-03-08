#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

int readOnConfigFile( char* filename, int** firstCol, double** secondCol ) {
  FILE* file = fopen( filename, "r" );
  if( file == NULL ) {
    return EXIT_FAILURE;
  }
  
  int linesCount = 0;
  while( !feof( file ) ) {
    char ch = fgetc( file );
    if( ch=='\n' ) {
      linesCount++;
    }
  }
  
  printf( "File has %d lines.\n", linesCount );
  *firstCol = malloc( linesCount * sizeof **firstCol );
  *secondCol = malloc( linesCount * sizeof **secondCol );
  
  rewind(file);
  
  for( int i = 0; i < linesCount; ++i ) {
    fscanf( file, "%d%lf", &(*firstCol)[i], &(*secondCol)[i] );
    //     printf( "read %d and %.10lf\n", (*firstCol)[i], (*secondCol)[i] );
  }
  fclose(file);
  
  return linesCount;
}

struct rparams {
  double* lambdas;
  double** actions;
  int* lengths;
  int nlambda;
};

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
    printf("eqn %d: %f\n", a, eqns[a]);
    gsl_vector_set( eqn, a, eqns[a] );
  }
  
  return GSL_SUCCESS;
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

int main( void ) {
  //   char* filename = "/data2/Results/GN/red/24x23x23/results_1/Configs/ScalarField_ScalarOnConfig_.47.dat";
  
  int nlambda = 3;
  double lambdas[nlambda];
  lambdas[0] = 0.46;
  lambdas[1] = 0.47;
  lambdas[2] = 0.48;
  //   int V = 24 * 23 * 23;
  
  char* sfNames[nlambda];
  sfNames[0] = "/data2/Results/GN/red/24x23x23/results_1/Configs/ScalarField_ScalarOnConfig_.46.dat";
  sfNames[1] = "/data2/Results/GN/red/24x23x23/results_1/Configs/ScalarField_ScalarOnConfig_.47.dat";
  sfNames[2] = "/data2/Results/GN/red/24x23x23/results_1/Configs/ScalarField_ScalarOnConfig_.48.dat";
  char* actionNames[nlambda];
  actionNames[0] = "/data2/Results/GN/red/24x23x23/results_1/Configs/BosonicAction_ScalarOnConfig_.46.dat";
  actionNames[1] = "/data2/Results/GN/red/24x23x23/results_1/Configs/BosonicAction_ScalarOnConfig_.47.dat";
  actionNames[2] = "/data2/Results/GN/red/24x23x23/results_1/Configs/BosonicAction_ScalarOnConfig_.48.dat";
  
  double* sfVals[nlambda];
  double* actionVals[nlambda];
  int length[nlambda];
  int numData = 0;
  
  for( int numLambda = 0; numLambda < nlambda; ++numLambda ){
    
    int* count1 = NULL;         // TODO: count is not needed anywhere, can we skip reading this?
    int* count2 = NULL;
    //     double* sfOnConf = sfVals[numLambda];
    //     double* actionOnConf = actionVals[numLambda];
    
    length[numLambda] = readOnConfigFile( sfNames[numLambda], &count1, &sfVals[numLambda] );
    int linesCountAction = readOnConfigFile( actionNames[numLambda], &count2, &actionVals[numLambda] );
    
    if( linesCountAction != length[numLambda] ) {
      printf("Error: files have different lengths.");
      return EXIT_FAILURE;
    }
    
    numData += length[numLambda];
    
    //TODO: thermalization
    free( count1 );
    free( count2 );
  }
  
  printf("We have %d data points.\n", numData);
  
  struct rparams p = {
    lambdas,
    actionVals,
    length,
    nlambda,
  };
  
  gsl_vector *fa = gsl_vector_alloc( nlambda-1 );
  gsl_vector *eqns = gsl_vector_alloc( nlambda-1 );
  for( int numLambda = 0; numLambda < nlambda-1; ++numLambda ) {
    gsl_vector_set( fa, numLambda, 10*(numLambda+1) );
  }
  
  // testing equation evaluation on initial values
  equation( fa, &p, eqns );
  
  
  const gsl_multiroot_fsolver_type *T;
  gsl_multiroot_fsolver *s;
  
  const size_t n = nlambda - 1;
  size_t iter = 0;
  int status;
  
  gsl_multiroot_function f = {&equation, n, &p};
  T = gsl_multiroot_fsolver_hybrids;
  s = gsl_multiroot_fsolver_alloc(T, n);
  gsl_multiroot_fsolver_set( s, &f, fa );
  
  print_state(iter, s);
  
  do
  {
    iter++;
    status = gsl_multiroot_fsolver_iterate (s);
    
    print_state (iter, s);
    
    if (status)   /* check if solver is stuck */
      break;
    
    status = 
    gsl_multiroot_test_residual (s->f, 1e-7);
  }
  while (status == GSL_CONTINUE && iter < 1000);
  
  printf ("status = %s\n", gsl_strerror (status));
  
  double fasSolution[nlambda];
  fasSolution[0] = 0.;
  for( int a = 1; a < nlambda; ++a ) {
    fasSolution[a] = gsl_vector_get( s->x, a-1 );
  }
  double result = calcObservable( 0.44, sfVals, &p, fasSolution );
  printf("value for 0.44: %.6f\n", result);
  
  gsl_multiroot_fsolver_free (s);
  
  for( int numLambda = 0; numLambda < nlambda-1; ++numLambda ){
//     printf( "equation value: %lf\n", gsl_vector_get( eqns, numLambda ) );
    free( sfVals[numLambda] );
    free( actionVals[numLambda] );
  }
  gsl_vector_free( fa );
  gsl_vector_free( eqns );
  return EXIT_SUCCESS;
}