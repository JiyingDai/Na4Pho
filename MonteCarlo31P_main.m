clear all; close all;


% update the address
load('...\fit_dummy.mat');
load('...\xaxis.mat');

M=max(real(fit_dummy(:)));
figure(1),set(gcf, 'Position',  [0, 50, 600, 400],'name','Original spectrum'),plot(xaxis,real(fit_dummy),'LineWidth',1.5),xlim([-20 20]),ylim([-0.2*M 1.2*M]),set(gca,'XDir','reverse');

signal=ifft(ifftshift(fit_dummy,1),[],1);

figure(2),plot(real(fftshift(fft(squeeze(signal).',[],1),1)),'LineWidth',1.5);


ks=size(fit_dummy,1);
Nch=16;

Namp=0.3*max(abs(signal(:)));
nsa=1;

fidPick=[1:5];
% Namp=0.1:0.1:1.5;

ratio_ture=zeros(size(Namp));
ratio_peak=zeros(size(Namp));
peak_fake=zeros(size(Namp));
E_chi=zeros(size(Namp));
E_peak_real_mag=zeros(size(Namp));
peak_real=zeros(size(Namp));


for nn=1:size(Namp,2)
    %define signal
    sens=zeros(ks,Nch);
    signal_true=zeros(ks,Nch);

    for ch=1:Nch
        sens(:,ch)=1*exp(1i*2*pi/Nch*ch);
        signal_true(:,ch)=signal.*sens(:,ch);
    end

    SNR_ch(nn)=Namp(nn)/max(abs(signal(:)));
    disp(strcat('The SNR per channel:',num2str(SNR_ch(nn))));


    sig_com_true=zeros(nsa,ks);
    sig_com=zeros(nsa,ks);
    Na_com=zeros(nsa,ks);
    signal_Rx=zeros(nsa,ks,Nch);
    Na_Rx=zeros(nsa,ks,Nch);

    tic

    for na=1:nsa
        %define noise
        noise=Namp(nn)*(randn(ks,Nch)+1i*randn(ks,Nch));
%         noise=Namp(nn)*(randn(ks,Nch));
        noise_Na=Namp(nn)*(randn(ks,Nch)+1i*randn(ks,Nch));
        

        %received signal
        signal_Rx(na,:,:)=signal_true+noise;
        Na_Rx(na,:,:)=signal_true+noise_Na; 

        %combine with the true signal
        S=(mean(signal_true(fidPick,:),1)).';
        U=pinv(sqrt(S'*S))*S';
        V=U*squeeze(signal_Rx(na,:,:)).';
        sig_com_true(na,:)=V;


        %combine with the received signal
        S=(mean(squeeze(signal_Rx(na,fidPick,:)),1)).';
        U=pinv(sqrt(S'*S))*S';
        V=U*squeeze(signal_Rx(na,:,:)).';
        sig_com(na,:)=V;

        %combine with the received signal
        S=(mean(squeeze(Na_Rx(na,fidPick,:)),1)).';
        U=pinv(sqrt(S'*S))*S';
        V=U*squeeze(signal_Rx(na,:,:)).';
        Na_com(na,:)=V;

    end

    toc
%% mean -> fft -> phasing
     % evaluate spectrum combined with perfect B1-
    Parameters.FirstOrdPhaseFunct=ones(256,1);
    sig_com_true_av=mean(sig_com_true,1);
    spec_com_true_av=PhaseCorrJD_SingleVox(fftshift(fft(sig_com_true_av.',[],1),1));
    % tested to be the same as above no matter ind fixed as 129 or not, so
    % phaseing first then FFT or vice versa does not matter


    M=max(real(spec_com_true_av(:)));
    figure(),set(gcf, 'Position',  [600, 50, 600, 400],'name','Opt combined spectrum'),plot(xaxis,real(spec_com_true_av),'LineWidth',1.5),xlim([-20 20]),ylim([-0.2*M 1.2*M]),set(gca,'XDir','reverse');

    % evaluate spectrum combined with averaged (received) FID
    sig_com_av=mean(sig_com,1);
    spec_com_sw_av=PhaseCorrJD_SingleVox(fftshift(fft(sig_com_av.',[],1),1));
%     % tested to be the same as above no matter ind fixed as 129 or not, so
    % phaseing first then FFT or vice versa does not matter
    
    M=max(real(spec_com_sw_av(:)));

    figure(),set(gcf, 'Position',  [1200, 50, 600, 400],'name','Sw combined spectrum'),plot(xaxis,real(spec_com_sw_av),'LineWidth',1.5),xlim([-20 20]),ylim([-0.2*M 1.2*M]),set(gca,'XDir','reverse');

    % evaluate spectrum combined with averaged (received) FID with
    % independent noise
    Na_com_av=mean(Na_com,1);
    spec_com_Na_av=PhaseCorrJD_SingleVox(fftshift(fft(Na_com_av.',[],1),1));
%     % tested to be the same as above no matter ind fixed as 129 or not, so
    % phaseing first then FFT or vice versa does not matter

    M=max(real(spec_com_Na_av(:)));
    figure(),set(gcf, 'Position',  [1800, 50, 600, 400],'name','Na combined spectrum'),plot(xaxis,real(spec_com_Na_av),'LineWidth',1.5),xlim([-20 20]),ylim([-0.2*M 1.2*M]),set(gca,'XDir','reverse');

 
    
    
end

