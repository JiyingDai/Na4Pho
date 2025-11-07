%% This code is to plot Figure S4 of the Na4P manuscript.
clear all; close all; clc;
%%
% General notations:
% P4P (phospurus for phosphorus): 
% corresponds to the self-weighted method, for which the combination 
% sensitity is the average of 2:5 FID points of the denoised 31P signal (i.e., denoised raw 31P data)

% dSP4SP (denoised shifted phosphorus for shifted phosphrus):
% corresponds to the x-axis of Figure S4, for which the combination
% sensitivity is the average of 2:5 FID points of the to-be-combined 31P signal (i.e., raw 31P data + added noise)

% Na4SP (sodium for shifted phosphorus):
% correspond to the 23Na for 31P method, for which the combination
% sensitivity is the average of 2:5 FID points of the same-resolution 23Na signal

% Load 4 reconstructed datasets from the big phantom.
% Each dataset contains the common information (mask,etc.) and the combined
% spectra using three types of sensitivities respectively.
% N20: 20 refers to the standard deviation of the added noise (to the raw 31P FID). Same for N40 to N160.

load("N20.mat")
load("N40.mat")
load("N80.mat")
load("N160.mat")



% load the common specs (same for all four datasets)
mask = N20.mask;
xaxis = N20.xaxis;
sig_window = N20.sig_window;
noisewindow = N20.noisewindow;
delta=-45;
delta_ppm = delta/256*(max(N20.xaxis)-min(N20.xaxis));
dim1 = size(N20.P4P_spectra,2);
dim2 = size(N20.P4P_spectra,3);
dim3 = size(N20.P4P_spectra,4);

% Concatenate spectra from all four datasets
dSP4SP_spectra = cat(5,N20.dSP4SP_spectra,N40.dSP4SP_spectra,N80.dSP4SP_spectra,N160.dSP4SP_spectra);
P4P_spectra = cat(5,N20.P4P_spectra,N40.P4P_spectra,N80.P4P_spectra,N160.P4P_spectra);
Na4SP_Corr_spectra = cat(5,N20.Na4SP_Corr_spectra,N40.Na4SP_Corr_spectra,N80.Na4SP_Corr_spectra,N160.Na4SP_Corr_spectra);

dSP4SP_SNR = zeros(dim1,dim2,dim3,4);
P4P_SNR = zeros(dim1,dim2,dim3,4);
Na4SP_Corr_SNR = zeros(dim1,dim2,dim3,4);
%%
% Calculate SNR according to the spectra
for nn=1:4
    dSP4SP_SNR(:,:,:,nn) = SNR_calc(sig_window,noisewindow,xaxis,dSP4SP_spectra(:,:,:,:,nn),delta_ppm);
    P4P_SNR(:,:,:,nn) = SNR_calc(sig_window,noisewindow,xaxis,P4P_spectra(:,:,:,:,nn),0);
    Na4SP_Corr_SNR(:,:,:,nn) = SNR_calc(sig_window,noisewindow,xaxis,Na4SP_Corr_spectra(:,:,:,:,nn),delta_ppm);
end
%% Plot Figure S4
cat_mask=repmat(mask,[1 1 1 4]);
figure(),set(gcf, 'Position',[00, 50, 500, 400],'name','Figure S4')
for k=1:size(P4P_SNR(:))
    plot(cat_mask(k)*dSP4SP_SNR(k),cat_mask(k)*P4P_SNR(k),'r.',cat_mask(k)*dSP4SP_SNR(k),cat_mask(k)*Na4SP_Corr_SNR(k),'b.'),xlim([0 40]),hold on
end
legend('SW,4-FID-points','^{23}Na-W')
%% CSI plot of a single z-slice with SNR background, P4P

% select dataset and combination method
spectra = N20.Na4SP_Corr_spectra; 

% select transversal slice
zslc = 5;

im = squeeze(max(real(spectra(:,:,:,zslc)),[],1));

figure(),set(gcf, 'Position',[550, 50, 700, 650],'color','r','name',strcat('Single transversal slice CSI'));

for nx = 1:dim1
    for ny = 1:dim2
        subplot('Position',[1/dim2*(ny-1) 1-nx*1/dim1 1/dim2  1/dim1]),plot(xaxis,real(squeeze(spectra(:,nx,ny,zslc))),'color',[0,im(nx,ny)/max(max(im)),1],'Linewidth',1.2),xlim([-20 20]),ylim([-0.2*max(max(im)) 1.2*max(max(im))]),set(gca,'XDir','reverse');
        axis square,set(subplot('Position',[1/dim2*(ny-1) 1-nx*1/dim1 1/dim2  1/dim1]),'xtick',[],'ytick',[],'xcolor','r','ycolor','r'),hold on; 
    end
end