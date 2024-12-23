function PhasedSPECTRA=PhaseCorrJD_SingleVox(SPECTRA)


PhasedSPECTRA=zeros(size(SPECTRA));

spectrum=SPECTRA;
phasevector=angle(spectrum);
% [~, ind]=max(abs(spectrum));
ind=129;
PhasedSPECTRA=complex(abs(spectrum).*cos(phasevector-phasevector(ind)) ,abs(spectrum).*sin(phasevector-phasevector(ind)) );


PhasedSPECTRA=reshape(PhasedSPECTRA,size(SPECTRA));

