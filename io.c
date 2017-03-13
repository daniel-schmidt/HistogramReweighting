#include "io.h"

size_t countLines( FILE* file ) {
  size_t linesCount = 0;
  while( !feof( file ) ) {
    char ch = fgetc( file );
    if( ch=='\n' ) {
      linesCount++;
    }
  }
  rewind(file);
  
  printf( "File has %zu lines.\n", linesCount );
  return linesCount;
}

size_t readLambdasFromFile( char const * const filename, double** lambdas ) {
  if( filename == NULL ) {
    printf("ERROR: readOnConfigFile got an empty file name.");
    exit(1);
  }
  
  FILE* file = fopen( filename, "r" );
  if( file == NULL ) {
    printf("ERROR: file not found for import in readLambdasFromFile: %s\n", filename);
    exit(1);
  }
  
  size_t linesCount = countLines( file );

  *lambdas = malloc( linesCount * sizeof *lambdas );
  if( *lambdas == NULL ) {
    printf("ERROR: memory allocation failed.");
    exit(1);
  }
  
  for( int i = 0; i < linesCount; ++i ) {
    fscanf( file, "%lf", *lambdas + i );
//     printf( "read %d and %.10lf\n", firstCol[i], secondCol[i+offset] );
  }
  fclose(file);
  
  return linesCount;
}

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

int readOnConfigFile( const size_t numThermal, char* filename, double** secondCol, size_t offset ) {
  if( filename == NULL ) {
    printf("ERROR: readOnConfigFile got an empty file name.");
    exit(1);
  }
  
  FILE* file = fopen( filename, "r" );
  if( file == NULL ) {
    printf("ERROR: file not found for import in readOnConfigFile: %s\n", filename);
    exit(1);
  }
  
  size_t linesCount = countLines( file );
  
  // skipping the first numThermal lines
  for( size_t line = 0; line < numThermal; ++line ) {
    fscanf(file, "%*[^\n]\n", NULL);
  }
  
  if( numThermal >= linesCount ) {
    printf("ERROR in readOnConfigFile: numThermal is larger than number of data lines!");
    exit(1);
  }
  
  linesCount -= numThermal;
    
  // enlarge the data array and append data from file to it
  double* newSecondCol = realloc( *secondCol, (offset + linesCount) * sizeof(double));
  if( newSecondCol == NULL ) {
    printf("ERROR: memory allocation failed.");
    exit(1);
  }
  *secondCol = newSecondCol;
  
  for( int i = 0; i < linesCount; ++i ) {
    fscanf( file, "%*d%lf", (*secondCol)+i+offset );    // read only second col and omit first
  }
  fclose(file);
  
  return linesCount;
}

size_t readData( const size_t numThermal, int nlambda, char** sfNames, double** sfVals, char** actionNames, double** actionVals, int* lengths ) {
  size_t offset = 0;
  for( int numLambda = 0; numLambda < nlambda; ++numLambda ){      
    lengths[numLambda] = readOnConfigFile( numThermal, sfNames[numLambda], sfVals, offset );
    int linesCountAction = readOnConfigFile( numThermal, actionNames[numLambda], actionVals, offset );
    
    if( linesCountAction != lengths[numLambda] ) {
      printf("ERROR: files have different lengths.");
      exit(1);
    }
    offset += lengths[numLambda];
  }  
  return offset;
}