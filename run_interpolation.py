import subprocess as subp
import time
import re
import glob
import numpy as np

def getLambdaFromFile(filename):
  findLambda=re.compile('\.\d+(?=.dat)')
  findName=re.compile('^(.*)_.*')
  m = findLambda.search(filename)
  n = findName.search(filename)
  return float(m.group(0)), n.group(1)


def run_interpolation( Nf, L, lambda_min, lambda_max, N_boot ):
  
  base_path = "/data2/Results/GN/red/%dx%dx%d/results_%d/Configs/" % (L, L-1, L-1, Nf)
  base_files = ["ScalarField_ScalarOnConfig_", "BosonicAction_ScalarOnConfig_"]
  out_names = []
  
  # load autocorrelation and filter for coupling range
  autocorr_file = base_path + "IntAutocorrTime.dat"
  autocorr = np.loadtxt(autocorr_file)
  ac_filter = np.logical_and(autocorr[:, 0] >= lambda_min, autocorr[:, 0] <= lambda_max)
  autocorr = autocorr[ac_filter]
  print(autocorr)
  out_names.append("IntAutocorrTime_%d_%d.dat" % (Nf, L))
  np.savetxt(out_names[-1], autocorr, "%.5f")

  # write paths to OnConfigFiles into files
  for base_file in base_files:
    files = glob.glob(base_path + base_file + "*")
    files.sort()
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
      print("ERROR: file not found")
      exit()
    out_names.append(out_name)
  
  
  #command = ["./multihist",
             #"IntAutocorrTime-%d-6.dat" % L,
             #"sf%d-6.txt" % L,
             #"action%d-6.txt" % L,
             #"%dx%dx%d/" % (L, L-1, L-1),
             #"%d" % L]
  command = ["./multihist",
            *out_names,
            "%dx%dx%d/" % (L, L-1, L-1),
            "%d" % L,
            "%d" % N_boot]
  time_start = time.time()
  #print(command)
  subp.run(command)
  time_end = time.time()
  print("Main loop took %.2f seconds." % (time_end-time_start) )
  
if __name__ == "__main__":
  LList = [10, 12, 16, 24]
  N_boot = 10
  for L in LList:
    run_interpolation(1, L, 0.46, 0.488, N_boot)