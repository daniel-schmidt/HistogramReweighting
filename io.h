#ifndef IO_H
#define IO_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

size_t countLines( FILE* file );

size_t readAutocorrFile( char const * const filename, double** lambdas, double** autocorr );

void readPathsFromFile( const char* filename, const size_t nlambda, char** paths );

int readOnConfigFile( const size_t numThermal, char* filename, double** secondCol, size_t offset );

size_t readData( const size_t numThermal, int nlambda, char** sfNames, double** sfVals, char** actionNames, double** actionVals, int* length );

#endif