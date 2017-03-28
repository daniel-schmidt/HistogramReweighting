#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include <time.h>

#include "io.h"
#include "single_run.h"

int main( int argc, char** argv ) {
  if( argc != 7 ) {
    printf( "ERROR: Need 6 input parameters: lambdas.txt sf_paths.txt action_paths.txt subfolder_name L N_boot\n" );
    exit(1);
  }
  
  if( sizeof(double) >= sizeof(long double) ){
    printf("WARNING: long double seems no longer than double: %lu, long double: %lu", sizeof(double), sizeof(long double));
  }
  const size_t L = atoi(argv[5]);
  const size_t V = L*(L-1)*(L-1);
  double* lambdas = 0;
  double* autocorr = 0;
  size_t nlambda = readAutocorrFile( argv[1], &lambdas, &autocorr );
  char* sfNames[nlambda];
  char* actionNames[nlambda];
  
  for( size_t n = 0; n < nlambda; ++n ) {
    printf("lamb=%.3f, 1/(1+2*autocorr)=%.3f\n", lambdas[n], autocorr[n]);
  }
  
  readPathsFromFile( argv[2], nlambda, sfNames );
  readPathsFromFile( argv[3], nlambda, actionNames );
  
  double* sfVals = NULL;
  double* actionVals = NULL;
  int lengths[nlambda];
  
  size_t numThermal = 200;
  size_t len_total = readData( numThermal, nlambda, sfNames, &sfVals, actionNames, &actionVals, lengths );
  printf("Skipped %zu thermalisation each, have a total of %zu data points.\n", numThermal, len_total);
  
  for( size_t a = 0; a < nlambda; ++a ) {    
    free( sfNames[a] );
    free( actionNames[a] );
  }
  
  // Set parameters and calculate solution
  size_t numInterpol = 101;
  
  struct rparams p = {
    lambdas,
    autocorr,
    actionVals,
    lengths,
    nlambda,
    len_total
  };
  
  double ip_lam   [numInterpol];
  double ip_sfabs [numInterpol];
  double ip_sus   [numInterpol];
  double ip_bc    [numInterpol];
  double ip_dlog  [numInterpol];
  
  single_run( V, &p, sfVals, numInterpol, ip_lam, ip_sfabs, ip_sus, ip_bc, ip_dlog );
  
  // binning and bootstrapping for error estimates
  srand(time(0));
  size_t bin_size = 100;
  size_t Nboot = atoi(argv[6]);
  
  double* actionSelect = malloc( len_total * sizeof *actionVals );
  double* sfSelect = malloc( len_total * sizeof *sfVals );
  
  p.actions = actionSelect;
  
  double* err_sfabs = calloc( numInterpol * sizeof(double), sizeof(double) );
  double* err_sus = calloc( numInterpol * sizeof(double), sizeof(double) );
  double* err_bc = calloc( numInterpol * sizeof(double), sizeof(double) );
  double* err_dlog = calloc( numInterpol * sizeof(double), sizeof(double) );
  
  double bin_ip_sfabs [numInterpol];
  double bin_ip_sus   [numInterpol];
  double bin_ip_bc    [numInterpol];
  double bin_ip_dlog  [numInterpol];
  
 
  const size_t numObservables = 4;
  char* filenames[numObservables];
  filenames[0] = "/BinnedScalarFieldAbs.dat";
  filenames[1] = "/BinnedSusceptibility.dat";
  filenames[2] = "/BinnedBinderCumulant.dat";
  filenames[3] = "/BinnedDLogScalarField.dat";
  FILE* files[numObservables];
  
  for( size_t k = 0; k < numObservables; ++k ) {
    char outpath[80] = { 0 };
    strcat( outpath, argv[4] );
    strcat( outpath, filenames[k] );
    char mkdircommand[80] = { 0 };
    strcat( mkdircommand, "mkdir -p ");
    strcat( mkdircommand, argv[4] );
    if( system( mkdircommand ) != 0 ) {
      puts( "ERROR: could not execute command to create file." );
    }
    files[k] = fopen( outpath, "w" );
    for( size_t ip = 0; ip < numInterpol; ++ip ) {
      fprintf( files[k],  "%.10f ", ip_lam[ip] );
    }
    fprintf( files[k], "\n");
  }
  
  for( size_t boot = 0; boot < Nboot; ++boot ) {
    printf( "Calculating bootstrap sample %zu...\n", boot );
    random_select( actionVals, sfVals, lengths, nlambda, bin_size, actionSelect, sfSelect );
    single_run( V, &p, sfSelect, numInterpol, ip_lam, bin_ip_sfabs, bin_ip_sus, bin_ip_bc, bin_ip_dlog );
    
    for( size_t ip = 0; ip < numInterpol; ++ip ) {
      err_sfabs[ip] += (bin_ip_sfabs[ip] - ip_sfabs[ip]) * (bin_ip_sfabs[ip] - ip_sfabs[ip]);
      err_sus[ip]   += (bin_ip_sus[ip] - ip_sus[ip])     * (bin_ip_sus[ip] - ip_sus[ip]);
      err_bc[ip]    += (bin_ip_bc[ip] - ip_bc[ip])       * (bin_ip_bc[ip] - ip_bc[ip]);
      err_dlog[ip]  += (bin_ip_dlog[ip] - ip_dlog[ip])   * (bin_ip_dlog[ip] - ip_dlog[ip]);
      
      fprintf( files[0],  "%.10f ", bin_ip_sfabs[ip] );
      fprintf( files[1],  "%.10f ", bin_ip_sus[ip] );
      fprintf( files[2],   "%.10f ", bin_ip_bc[ip] );
      fprintf( files[3], "%.10f ", bin_ip_dlog[ip] );
    }
    
    for( size_t k = 0; k < numObservables; ++k ) {
      fprintf( files[k], "\n");
    }
  }
  
  for( size_t k = 0; k < numObservables; ++k ) {
    fclose( files[k] );
  }
  
  // Writing full interpolations with error to files
  for( size_t ip = 0; ip < numInterpol; ++ip ) {
    err_sfabs[ip] = sqrt( err_sfabs[ip] / Nboot );
    err_sus[ip]   = sqrt( err_sus[ip]   / Nboot );
    err_bc[ip]    = sqrt( err_bc[ip]    / Nboot );
    err_dlog[ip]  = sqrt( err_dlog[ip]  / Nboot );
  }
  
  filenames[0] = "/InterpolScalarFieldAbs.dat";
  filenames[1] = "/InterpolSusceptibility.dat";
  filenames[2] = "/InterpolBinderCumulant.dat";
  filenames[3] = "/InterpolDLogScalarField.dat";
  
  for( size_t k = 0; k < numObservables; ++k ) {
    char outpath[80] = { 0 };
    strcat( outpath, argv[4] );
    strcat( outpath, filenames[k] );
    files[k] = fopen( outpath, "w" );
  }
  
  for( size_t ip = 0; ip < numInterpol; ++ip ) {
    fprintf( files[0], "%.10f %.10f %.10f \n", ip_lam[ip], ip_sfabs[ip], err_sfabs[ip] );
    fprintf( files[1], "%.10f %.10f %.10f \n", ip_lam[ip], ip_sus[ip],   err_sus[ip] );
    fprintf( files[2], "%.10f %.10f %.10f \n", ip_lam[ip], ip_bc[ip],    err_bc[ip] );
    fprintf( files[3], "%.10f %.10f %.10f \n", ip_lam[ip], ip_dlog[ip],  err_dlog[ip] );
  }
  
  for( size_t k = 0; k < numObservables; ++k ) {
    fclose( files[k] );
  }
  
  // Cleanup
  free( sfVals );
  free( actionVals );
  free( actionSelect );
  free( sfSelect );
  free( err_sfabs );
  free( err_sus );
  free( err_bc );
  free( err_dlog );
  
  free( lambdas );
  freeSolver();
  
  return EXIT_SUCCESS;
}