import subprocess as subp
import time
import re
import glob
import numpy as np

def getLambdaFromFile(filename):
  """Use a regular expression to find the value of the coupling in the filename.
  
     The file is expected to have the coupling value at the end of the filename
     as a number preceeded by '.dat'.
  """
  findLambda=re.compile('\d*\.\d+(?=.dat)')
  findName=re.compile('^(.*)_.*')
  m = findLambda.search(filename)
  n = findName.search(filename)
  try:
    lam = float(m.group(0))
    name = n.group(1)
  except AttributeError:
    print("No matching pattern to extract coupling found in filename: %s" % filename)
    exit(1)    
  return  lam, name


def run_interpolation( Nf, L, lambda_min, lambda_max, N_boot, N_thermal, f0, N_interpol ):  
  """ Creates the necessary input files and runs the C code to obtain interpolations.
  """
  base_path = "/data2/Results/GN/red/%dx%dx%d/results_%d/Configs/" % (L, L-1, L-1, Nf)
  base_files = ["ScalarField_ScalarOnConfig_", "BosonicAction_ScalarOnConfig_"]
  out_names = []
  
  # load autocorrelation and filter for coupling range
  #try:
    #autocorr = np.loadtxt(autocorr_file)
  #except FileNotFoundError:
    # If the autocorrelation file does not exist, we call the octave script to calculate it.
  #print("Autocorrelation must be computed beforehand!")
  command = [ "./calcAutocorr", str(Nf), str(L), str(N_thermal)]
  subp.run(command)
  autocorr_file = base_path + "IntAutocorrTime.dat"
  autocorr = np.loadtxt(autocorr_file)
  
  # filter out the autocorrelations where the coupling is in the given range
  ac_filter = np.logical_and(autocorr[:, 0] >= lambda_min, autocorr[:, 0] <= lambda_max)
  autocorr = autocorr[ac_filter]
  out_names.append("IntAutocorrTime_%d_%d.dat" % (Nf, L))
  np.savetxt(out_names[-1], autocorr, "%.5f")

  # write paths to OnConfigFiles into files
  for base_file in base_files:
    files = glob.glob(base_path + base_file + "*")
    files.sort()
    print(getLambdaFromFile(files[0]))
    # filter list for couplings in the given range
    files = [filename for filename in files if (
        (getLambdaFromFile(filename)[0] >= lambda_min)
        and
        (getLambdaFromFile(filename)[0] <= lambda_max)
      ) ]
  
    out_name = base_file + "%d_%d.txt" % (Nf, L)
    if files:
      with open(out_name, 'w' ) as file_handle:
        for filename in files:
          print(filename, file=file_handle)
    else:
      print("ERROR: No file found matching given coupling range from %.3f to %.3f." % (lam_min, lam_max) )
      print("File name patter is: " + base_path + base_file)
      exit()
    out_names.append(out_name)

  # set up C command and execute it
  command = ["./multihist",
            *out_names,
            "%dx%dx%d/" % (L, L-1, L-1),
            str(L),
            str(N_boot),
            str(N_thermal),
            str(f0),
            str(N_interpol)]
  #time_start = time.time()
  #subp.run(command)
  process = subp.Popen(command)
  #time_end = time.time()
  #print("Main loop took %.2f seconds." % (time_end-time_start) )
  legend = ["Nf", "L", "lambda_min", "lambda_max", "N_boot", "N_thermal", "f0", "N_interpol"]
  parameters = [Nf, L, lambda_min, lambda_max, N_boot, N_thermal, f0, N_interpol]
  with open( "%dx%dx%d/parameters.dat" % (L, L-1, L-1), 'w' ) as file_handle:
    print( parameters, file = file_handle)
  
  return process
  
if __name__ == "__main__":
  #LList = [10, 12, 16, 20, 24]
  #f0List = [-40, -70, -170, -350, -600]
  #LList = [8]
  #f0List = [-20] 
  #Nf = 1
  #lam_min = 0.44
  #lam_max = 0.52
  
  #LList = [8, 10, 12, 16, 20, 24]
  #f0List = [-20, -40, -70, -170, -350, -475]
  #LList = [24]
  #f0List = [-475]
  #Nf = 2
  #lam_min = 1.136
  #lam_max = 1.312
  
  LList = [8, 12, 16, 20, 24]
  f0List = [-20, -70, -170, -350, -475]
  Nf = 8
  lam_min = 5.312
  lam_max = 6.208
  
  N_boot = 10
  N_thermal = 200
  N_interpol = 101
  
  time_start = time.time()
  proc_list = []
  for f0, L in zip(f0List, LList):
    proc_list.append(run_interpolation(Nf, L, lam_min, lam_max, N_boot, N_thermal, f0, N_interpol))
  
  for proc in proc_list:
    proc.wait()
  time_end = time.time()
  print("Total execution took %.2f seconds." % (time_end-time_start) )