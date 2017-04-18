function means = autocorr_GN(Nf, L, numThermal)
  lat_size_string = [ num2str(L) "x" num2str(L-1) "x" num2str(L-1)];
  base_path=['/data2/Results/GN/red/' lat_size_string '/results_' num2str(Nf) '/Configs/'];
  observable_name = 'ScalarField_ScalarOnConfig_';
  paths=glob([base_path observable_name '*.dat']);

  numFiles = size(paths)(1);
  means = zeros(numFiles, 6);
  for i = 1:numFiles
    endIdx = rindex(paths{i},".")-1;
    frontIdx = rindex(paths{i},"_")+1;
    coupling = str2num(paths{i}( frontIdx : endIdx ) );
    data=load(paths{i});

    dataEqui=abs(data(numThermal:end, 2));
    [meanval, err, derr, tau, dtau] = UWerrTexp( dataEqui, [], [], 1 );
    means(i, :) = [coupling, meanval, err, derr, tau, dtau];
  endfor

%    errorbar(means(:,1),means(:,5),means(:,6))
%    errorbar(means(:,1),means(:,2),means(:,3))
  
%    out_name = ['Autocorrelations/tau_int_' num2str(Nf) '_' num2str(L) '.dat']
  save_header_format_string ("# coupling, meanval, error, error of error, int. autocorr. time, autocorr. error");
  out_name = [base_path 'IntAutocorrTime.dat'];
  save(out_name, "means");