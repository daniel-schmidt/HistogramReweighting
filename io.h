#ifndef IO_H
#define IO_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

void readPathsFromFile( const char* filename, const size_t nlambda, char** paths );

int readOnConfigFile( char* filename, int** firstCol, double** secondCol );

void readData( int nlambda, char** sfNames, double** sfVals, char** actionNames, double** actionVals, int* length );

#endif