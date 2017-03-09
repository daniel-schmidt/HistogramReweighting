#ifndef IO_H
#define IO_H

#include <stdlib.h>
#include <stdio.h>

int readOnConfigFile( char* filename, int** firstCol, double** secondCol );

void readData( int nlambda, char** sfNames, double** sfVals, char** actionNames, double** actionVals, int* length );

#endif