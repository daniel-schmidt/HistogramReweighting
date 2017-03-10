#ifndef IO_H
#define IO_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

size_t countLines( FILE* file );

size_t readLambdasFromFile( char const * const filename, double** lambdas );

void readPathsFromFile( const char* filename, const size_t nlambda, char** paths );

int readOnConfigFile( char* filename, int* firstCol, double** secondCol, size_t offset );

size_t readData( int nlambda, char** sfNames, double** sfVals, char** actionNames, double** actionVals, int* length );

#endif