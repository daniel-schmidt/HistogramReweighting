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
  printf ("iter = %3u x = % .3f % .3f "
  "f(x) = % .3e % .3e\n",
          iter,
          gsl_vector_get (s->x, 0), 
          gsl_vector_get (s->x, 1),
          gsl_vector_get (s->f, 0), 
          gsl_vector_get (s->f, 1));
}

int equation( const gsl_vector * x, void * params, gsl_vector *eqn ) {
  double* lambdas = ( ( struct rparams* ) params )->lambdas;
  double** actions = ( ( struct rparams* ) params )->actions;
  int* lengths = ( ( struct rparams* ) params )->lengths;
  int nlambda = ( ( struct rparams * ) params )->nlambda;
  
  double fas[nlambda];
  double eqns[nlambda];
  for( int a = 0; a < nlambda; ++a ) {
    fas[a] = gsl_vector_get( x, a );
    eqns[a] = -fas[a];
  }
  
  for( int c = 0; c < nlambda; ++c ) {
    double sum = 0.;
    for( int b = 0; b < nlambda; ++b ) {
      for( int i = 0; i < lengths[b]; ++i ) {
        double denom = 0.;
        for( int a = 0; a < nlambda; ++a ) {
          denom += lengths[a] * exp(actions[b][i] * (lambdas[c] - lambdas[a]) + fas[a]);
        }
        sum += 1./denom;
      }
    }
    eqns[c] = fas[c] + log(sum);
  }
  
  for( int a = 0; a < nlambda; ++a ) {
    printf("eqn %d: %f\n", a, eqns[a]);
    gsl_vector_set( eqn, a, eqns[a] );
  }
  
  return GSL_SUCCESS;
}


int main( void ) {
  //   char* filename = "/data2/Results/GN/red/24x23x23/results_1/Configs/ScalarField_ScalarOnConfig_.47.dat";
  
  int nlambda = 2;
  double lambdas[nlambda];
  lambdas[0] = 0.46;
  lambdas[1] = 0.47;
  //   int V = 24 * 23 * 23;
  
  char* sfNames[nlambda];
  sfNames[0] = "/data2/Results/GN/red/24x23x23/results_1/Configs/ScalarField_ScalarOnConfig_.46.dat";
  sfNames[1] = "/data2/Results/GN/red/24x23x23/results_1/Configs/ScalarField_ScalarOnConfig_.47.dat";
  char* actionNames[2];
  actionNames[0] = "/data2/Results/GN/red/24x23x23/results_1/Configs/BosonicAction_ScalarOnConfig_.46.dat";
  actionNames[1] = "/data2/Results/GN/red/24x23x23/results_1/Configs/BosonicAction_ScalarOnConfig_.47.dat";
  
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
  
  gsl_vector *fa = gsl_vector_alloc( nlambda );
  gsl_vector *eqns = gsl_vector_alloc( nlambda );
  for( int numLambda = 0; numLambda < nlambda; ++numLambda ) {
    gsl_vector_set( fa, numLambda, 10*numLambda );
  }
  
  // testing equation evaluation on initial values
  equation( fa, &p, eqns );
  
  for( int numLambda = 0; numLambda < nlambda; ++numLambda ){
    printf( "equation value: %lf\n", gsl_vector_get( eqns, numLambda ) );
    free( sfVals[numLambda] );
    free( actionVals[numLambda] );
  }
  
  const gsl_multiroot_fsolver_type *T;
  gsl_multiroot_fsolver *s;
  
  const size_t n = 2;
  size_t iter = 0;
  int status;
  
  gsl_multiroot_function f = {&equation, n, &p};
  T = gsl_multiroot_fsolver_hybrids;
  s = gsl_multiroot_fsolver_alloc(T, 2);
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
  
  gsl_multiroot_fsolver_free (s);
  
  gsl_vector_free( fa );
  gsl_vector_free( eqns );
  return EXIT_SUCCESS;
}