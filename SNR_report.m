function SNR_report(SNR,dim1,dim2,dim3,mask)

% voxel_tot = dim1*dim2*dim3;
voxel_tot = nnz(mask);

SNR_threshold = 8;

SNR=SNR.*mask;
SNR_tot = sum(SNR(:),'all');
SNR8upNum = size(find(SNR(:)>=SNR_threshold),1);
SNR8upAVG = sum(SNR(find(SNR(:)>=SNR_threshold)))/SNR8upNum;


disp(['The total SNR of ',num2str(voxel_tot),' voxels is ',num2str(SNR_tot),'.']);
disp(['There are ',num2str(SNR8upNum),' voxels with SNR >= ',num2str(SNR_threshold),'.']);
disp(['The average SNR of the qualified voxels (where SNR>=8) is ',num2str(SNR8upAVG)])
fprintf('\n');
end