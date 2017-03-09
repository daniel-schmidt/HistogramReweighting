#include "io.h"

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

void readData( int nlambda, char** sfNames, double** sfVals, char** actionNames, double** actionVals, int* length ) {
  for( int numLambda = 0; numLambda < nlambda; ++numLambda ){
    
    int* count1 = NULL;         // TODO: count is not needed anywhere, can we skip reading this?
    int* count2 = NULL;
    
    length[numLambda] = readOnConfigFile( sfNames[numLambda], &count1, &sfVals[numLambda] );
    int linesCountAction = readOnConfigFile( actionNames[numLambda], &count2, &actionVals[numLambda] );
    
    if( linesCountAction != length[numLambda] ) {
      printf("Error: files have different lengths.");
      exit(1);
    }
    
    //TODO: thermalization
    free( count1 );
    free( count2 );
  }
}