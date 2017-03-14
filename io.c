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

size_t readAutocorrFile( char const * const filename, double** lambdas, double** autocorr ) {
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

  *lambdas  = malloc( linesCount * sizeof *lambdas );
  *autocorr = malloc( linesCount * sizeof *autocorr );
  if( *lambdas == 0 || *autocorr == 0 ) {
    printf("ERROR: memory allocation failed.");
    exit(1);
  }
   
  for( int i = 0; i < linesCount; ++i ) {
    int read = fscanf( file, "%lf %*f %*f %*f %lf %*f", *lambdas + i, *autocorr + i );
    if( read != 2 ) {
      puts( "ERROR in readAutocorrFile: read a wrong number of items from autocorrelation file." );
      exit(1);
    }
    (*autocorr)[i] = 1./(1. + 2. * (*autocorr)[i] );
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
  
  for( size_t i = 0; i < nlambda; ++i ) {
    char* line = NULL;
    int read = getline( &line, &len, file );   // allocates memory for line
    if( read == -1 ) {
      puts( "ERROR reading line from PathFile." );
      exit(1);
    }
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
    int read = fscanf(file, "%*[^\n]\n");
    if( read != 0 ) {
      puts( "ERROR in readOnConfigFile: could not skip first lines" );
      exit(1);
    }
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
    int read = fscanf( file, "%*d%lf", (*secondCol)+i+offset );    // read only second col and omit first
    if( read != 1 ) {
      puts( "ERROR in readOnConfigFile: read a wrong number of items from on-config-file." );
      exit(1);
    }
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