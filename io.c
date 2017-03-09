#include "io.h"

void readPathsFromFile( const char* filename, const size_t nlambda, char** paths ) {
  FILE* file = fopen( filename, "r" );
  if( file == NULL ) {
    printf("ERROR: file %s not found for import.", filename);
    exit(1);
  }
  size_t len = 0;
  ssize_t read;
  
  for( size_t i = 0; i < nlambda; ++i ) {
    char* line = NULL;
    read = getline( &line, &len, file );   // allocates memory for line
    line[ strcspn( line, "\n" ) ] = 0;     // remove trailing newline
    printf("read line: %s\n", line);
    paths[i] = line;
  }
  fclose(file);
}

int readOnConfigFile( char* filename, int** firstCol, double** secondCol ) {
  if( filename == NULL ) {
    printf("ERROR: readOnConfigFile got an empty file name.");
    exit(1);
  }
  
  FILE* file = fopen( filename, "r" );
  if( file == NULL ) {
    printf("ERROR: file not found for import: %s\n", filename);
    exit(1);
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
      printf("ERROR: files have different lengths.");
      exit(1);
    }
    
    //TODO: thermalization
    free( count1 );
    free( count2 );
  }
}