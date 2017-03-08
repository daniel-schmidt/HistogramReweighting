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
  double* denom;
  int V;
};

int equation( const gsl_vector * x, void * params, gsl_vector *eqn ) {
  double* lambdas = ( ( struct rparams* ) params )->lambdas;
  double** actions = ( ( struct rparams* ) params )->actions;
  int* lengths = ( ( struct rparams* ) params )->lengths;
  int nlambda = ( ( struct rparams * ) params )->nlambda;
  double* denom = ( ( struct rparams* ) params )->denom;
  int V = ( ( struct rparams* ) params )->V;
    
  double fas[nlambda];
  double eqns[nlambda];
  for( int a = 0; a < nlambda; ++a ) {
    fas[a] = gsl_vector_get( x, a );
    eqns[a] = -fas[a];
  }
  
  int offset = 0;
  for( int b = 0; b < nlambda; ++b ){
    int currLen = lengths[b];
    for( int i = 0; i < currLen; ++i ){
      for( int a = 0; a < nlambda; ++a ){ 
        denom[ offset + i ] += lengths[a] * exp(-lambdas[a] * actions[b][i] + fas[a]);
      }
    }
    offset += currLen;
  }
  
  for( int a = 0; a < nlambda; ++a ){
    double fullSum = 0;
    offset = 0;
    for( int b = 0; b < nlambda; ++b ){
      int currLen = lengths[b];
      for( int i = 0; i < currLen; ++i ){
        fullSum += exp( -lambdas[a] * actions[b][i] ) / denom[ offset + i ];
      }
      offset += currLen;
    }
    eqns[a] -= log( fullSum );
  }
  
  for( int a = 0; a < nlambda; ++a ) {
    printf("eqn %d: %f", a, eqns[a]);
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
  int V = 24 * 23 * 23;
  
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
  
  double* denom = malloc( numData * sizeof (double) );
  printf("We have %d data points.\n", numData);
  
  struct rparams p = {
    lambdas,
    actionVals,
    length,
    nlambda,
    denom,
    V
  };
  
  gsl_vector *fa = gsl_vector_alloc( nlambda );
  gsl_vector *eqns = gsl_vector_alloc( nlambda );
  for( int numLambda = 0; numLambda < nlambda; ++numLambda ) {
    gsl_vector_set( fa, numLambda, 50*numLambda );
  }
  
  equation( fa, &p, eqns );
  
  for( int numLambda = 0; numLambda < nlambda; ++numLambda ){
    printf( "equation value: %lf\n", gsl_vector_get( eqns, numLambda ) );
    free( sfVals[numLambda] );
    free( actionVals[numLambda] );
  }
  free( denom);
  gsl_vector_free( fa );
  gsl_vector_free( eqns );
  return EXIT_SUCCESS;
}