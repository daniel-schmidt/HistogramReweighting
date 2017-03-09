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

int readOnConfigFile( char* filename, int* firstCol, double** secondCol, size_t offset ) {
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
  firstCol = malloc( linesCount * sizeof *firstCol );
  double* newSecondCol = realloc( *secondCol, (offset + linesCount) * sizeof(double));
  if( newSecondCol == NULL ) {
    printf("ERROR: memory allocation failed.");
    exit(1);
  }
  *secondCol = newSecondCol;
  rewind(file);
  
  for( int i = 0; i < linesCount; ++i ) {
    fscanf( file, "%d%lf", firstCol + i, (*secondCol)+i+offset );
//     printf( "read %d and %.10lf\n", firstCol[i], secondCol[i+offset] );
  }
  fclose(file);
  
  return linesCount;
}

size_t readData( int nlambda, char** sfNames, double** sfVals, char** actionNames, double** actionVals, int* lengths ) {
  size_t offset = 0;
  for( int numLambda = 0; numLambda < nlambda; ++numLambda ){
    
    int* count1 = NULL;         // TODO: count is not needed anywhere, can we skip reading this?
    int* count2 = NULL;
       
    lengths[numLambda] = readOnConfigFile( sfNames[numLambda], count1, sfVals, offset );
    int linesCountAction = readOnConfigFile( actionNames[numLambda], count2, actionVals, offset );
    
    if( linesCountAction != lengths[numLambda] ) {
      printf("ERROR: files have different lengths.");
      exit(1);
    }
    offset += lengths[numLambda];
    //TODO: thermalization
    free( count1 );
    free( count2 );
  }  
  return offset;
}