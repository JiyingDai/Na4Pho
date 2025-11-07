function SNR = SNR_calc(sig_window,noisewindow,xaxis,spectra,delta_ppm)

disp('Calculating SNR')
noisewindow = find(xaxis>5+delta_ppm & xaxis<20+delta_ppm);
sig_window = find(xaxis>-5+delta_ppm & xaxis<5+delta_ppm); % for phantom signal

noisemap = squeeze(std(real(spectra(noisewindow,:,:,:))));

SNR = squeeze(max(real(spectra(sig_window,:,:,:)),[],1))./noisemap;

disp('Finished SNR');fprintf('\n');

end