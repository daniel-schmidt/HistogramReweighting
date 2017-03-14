#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

#include "io.h"
#include "single_run.h"

int main( int argc, char** argv ) {
  if( argc != 4 ) {
    printf( "ERROR: Need 3 input parameters: lambdas.txt sf_paths.txt action_paths.txt\n" );
    exit(1);
  }
  
  if( sizeof(double) >= sizeof(long double) ){
    printf("WARNING: long double seems no longer than double: %lu, long double: %lu", sizeof(double), sizeof(long double));
  }
  
  double* lambdas = 0;
  double* autocorr = 0;
  size_t nlambda = readAutocorrFile( argv[1], &lambdas, &autocorr );
  char* sfNames[nlambda];
  char* actionNames[nlambda];
  
  for( size_t n = 0; n < nlambda; ++n ) {
    printf("lamb=%.3f, autocorr=%.3f\n", lambdas[n], autocorr[n]);
  }
  
  readPathsFromFile( argv[2], nlambda, sfNames );
  readPathsFromFile( argv[3], nlambda, actionNames );
  
  double* sfVals = NULL;
  double* actionVals = NULL;
  int lengths[nlambda];
  
  size_t numThermal = 9000;
  size_t len_total = readData( numThermal, nlambda, sfNames, &sfVals, actionNames, &actionVals, lengths );
  printf("Read a total of %zu data points.\n", len_total);
  
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
  
  single_run( &p, sfVals, numInterpol, ip_lam, ip_sfabs, ip_sus, ip_bc, ip_dlog );
  
  // binning and bootstrapping for error estimates
  int seed = 12;
  srand(seed);
  size_t bin_size = 100;
  size_t Nboot = 0;
  
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
  
  FILE* fileBinAbs = fopen( "BinnedScalarFieldAbs.dat", "w" );  
  FILE* fileBinSus = fopen( "BinnedSusceptibility.dat", "w" );
  FILE* fileBinBC = fopen( "BinnedBinderCumulant.dat", "w" );
  FILE* fileBinDlog = fopen( "BinnedDLogScalarField.dat", "w" );
  
  for( size_t ip = 0; ip < numInterpol; ++ip ) {
    fprintf( fileBinAbs,  "%.10f ", ip_lam[ip] );
    fprintf( fileBinSus,  "%.10f ", ip_lam[ip] );
    fprintf( fileBinBC,   "%.10f ", ip_lam[ip] );
    fprintf( fileBinDlog, "%.10f ", ip_lam[ip] );
  }
  fprintf( fileBinAbs, "\n");
  fprintf( fileBinSus, "\n");
  fprintf( fileBinBC, "\n");
  fprintf( fileBinDlog, "\n");
  
  for( size_t boot = 0; boot < Nboot; ++boot ) {
    printf( "Calculating bootstrap sample %zu...\n", boot );
    random_select( actionVals, sfVals, lengths, nlambda, bin_size, actionSelect, sfSelect );
    single_run( &p, sfSelect, numInterpol, ip_lam, bin_ip_sfabs, bin_ip_sus, bin_ip_bc, bin_ip_dlog );
    
    for( size_t ip = 0; ip < numInterpol; ++ip ) {
      err_sfabs[ip] += (bin_ip_sfabs[ip] - ip_sfabs[ip]) * (bin_ip_sfabs[ip] - ip_sfabs[ip]);
      err_sus[ip]   += (bin_ip_sus[ip] - ip_sus[ip])     * (bin_ip_sus[ip] - ip_sus[ip]);
      err_bc[ip]    += (bin_ip_bc[ip] - ip_bc[ip])       * (bin_ip_bc[ip] - ip_bc[ip]);
      err_dlog[ip]  += (bin_ip_dlog[ip] - ip_dlog[ip])   * (bin_ip_dlog[ip] - ip_dlog[ip]);
      
      fprintf( fileBinAbs,  "%.10f ", bin_ip_sfabs[ip] );
      fprintf( fileBinSus,  "%.10f ", bin_ip_sus[ip] );
      fprintf( fileBinBC,   "%.10f ", bin_ip_bc[ip] );
      fprintf( fileBinDlog, "%.10f ", bin_ip_dlog[ip] );
    }
    fprintf( fileBinAbs, "\n");
    fprintf( fileBinSus, "\n");
    fprintf( fileBinBC, "\n");
    fprintf( fileBinDlog, "\n");
  }
  
  fclose( fileBinAbs );
  fclose( fileBinSus );
  fclose( fileBinBC );
  fclose( fileBinDlog );
  
  for( size_t ip = 0; ip < numInterpol; ++ip ) {
    err_sfabs[ip] = sqrt( err_sfabs[ip] / Nboot );
    err_sus[ip]   = sqrt( err_sus[ip]   / Nboot );
    err_bc[ip]    = sqrt( err_bc[ip]    / Nboot );
    err_dlog[ip]  = sqrt( err_dlog[ip]  / Nboot );
  }
  
  FILE * fileAbs = fopen("InterpolScalarFieldAbs.dat", "w");
  FILE * fileSus = fopen("InterpolSusceptibility.dat", "w");
  FILE * fileBC = fopen("InterpolBinderCumulant.dat", "w");
  FILE * fileDlog = fopen("InterpolDLogScalarField.dat", "w");
  
  for( size_t ip = 0; ip < numInterpol; ++ip ) {
    fprintf( fileAbs,  "%.10f %.10f %.10f \n", ip_lam[ip], ip_sfabs[ip], err_sfabs[ip] );
    fprintf( fileSus,  "%.10f %.10f %.10f \n", ip_lam[ip], ip_sus[ip],   err_sus[ip] );
    fprintf( fileBC,   "%.10f %.10f %.10f \n", ip_lam[ip], ip_bc[ip],    err_bc[ip] );
    fprintf( fileDlog, "%.10f %.10f %.10f \n", ip_lam[ip], ip_dlog[ip],  err_dlog[ip] );
  }
  
  fclose( fileAbs );
  fclose( fileSus );
  fclose( fileBC );
  fclose( fileDlog );
  
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