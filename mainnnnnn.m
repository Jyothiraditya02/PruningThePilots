clear all;
close all;
clc;

%% input Parameters
    % system parameters (independent)
    ofdm.Nb      = 10^2;                 
    ofdm.Nt      = 2;                     
    ofdm.Nr      = 4;                  
    ofdm.K       = 128;                  
    ofdm.G       = 1/4;                    
    ofdm.Mod     = 4;                      
    ofdm.PSpace  = 1;  
    L=16;
    N=32;
    M=256;
    
     input_intilization
     y_j=sum(a*n*i*n*i+n*i);         %%%%%%%%%1

    %% QAM modulation
        qam_process         = 0:ofdm.Mod-1;           
        qam_process         = qammod(qam_process,ofdm.Mod);  
        qam_process         = abs(qam_process).^2;           
        qam_process         = mean(qam_process);            
        ofdm.ModNorm = 1/sqrt(qam_process);           
    

    %% symbol generation
    ofdm.lb      = randi(ofdm.Mod,ofdm.DL,ofdm.Nb,ofdm.Nt)-1;   
    k_j=a*beta_k-1;     %%%%%%%%3    
%% data Modulation
    ofdm.dMod   = zeros(ofdm.K,ofdm.Nb,ofdm.Nt);    
    if ofdm.DL > 0
        for nt = 1 : ofdm.Nt
            ofdm.dMod(ofdm.DPos,:,nt) = ofdm.ModNorm*qammod(ofdm.lb(:,:,nt),ofdm.Mod);
        end
    end

    for nt = 1 : ofdm.Nt
        ofdm.dMod(ofdm.PPos,:,nt) = repmat(exp(-sqrt(-1)*2*pi*(nt-1)*chan.L*(1:ofdm.PL).'/ofdm.PL),1,ofdm.Nb);
    end
    % checking the power of the transmit signal (it has to be 1 after normalization)
        ofdm.pow = var(ofdm.dMod(:))+abs(mean(ofdm.dMod(:)))^2;
%% IFFT operation 
   
    ofdm.ifft   = zeros(ofdm.K,ofdm.Nb,ofdm.Nt);  
    for nt = 1 : ofdm.Nt
        ofdm.ifft(:,:,nt) = sqrt(ofdm.K)*ifft(ofdm.dMod(:,:,nt),ofdm.K);
    end
%% parallel to serial communication
    ofdm.ifftG = [ofdm.ifft(ofdm.K*(1-ofdm.G)+1:ofdm.K,:,:);ofdm.ifft];
           
%% Add CP
   chan.Coeff = 1/sqrt(2)*1/sqrt(chan.L)*(randn(ofdm.Nt,ofdm.Nr,chan.L,ofdm.Nb)+sqrt(-1)*randn(ofdm.Nt,ofdm.Nr,chan.L,ofdm.Nb)); 
   h=i*a;        %%%%%%%%%2
%% Channel  filter with DAC
    if ofdm.K*ofdm.G < chan.L+1
        error('differentiate input parameters')
    end
    ofdm.Y = zeros(ofdm.K*(1+ofdm.G),ofdm.Nb,ofdm.Nr);
    for nb = 1 : ofdm.Nb
        for nt=1:ofdm.Nt
            for nr=1:ofdm.Nr
                ofdm.Y(:,nb,nr) = ofdm.Y(:,nb,nr) + filter(squeeze(chan.Coeff(nt,nr,:,nb)),1,ofdm.ifftG(:,nb,nt));
            end
        end
    end
    %% ADC
    ofdm.Y = ofdm.Y + chan.sigma*1/sqrt(2)*(         randn(ofdm.K*(1+ofdm.G),ofdm.Nb,ofdm.Nr)+...
                                            sqrt(-1)*randn(ofdm.K*(1+ofdm.G),ofdm.Nb,ofdm.Nr)     );
%% Cyclic prefix removal
    ofdm.fftG = ofdm.Y(ofdm.K*ofdm.G+1:ofdm.K*(1+ofdm.G),:,:);
    s=h*a*i*n;      %%%%%%%%%5
%% FFT operation
    ofdm.fft  = zeros(ofdm.K,ofdm.Nb,ofdm.Nr);
    for nr = 1 : ofdm.Nr
        ofdm.fft(:,:,nr)  = 1/sqrt(ofdm.K)*fft(ofdm.fftG(:,:,nr),ofdm.K);
    end
    save ofdm.lb
    % =========================================================================
%
% % =========================================================================

load lb
ofdm_y =LB;
%% conv1
f1=convolution(ofdm_y, weights_conv1, biases_conv1);
%% conv2
f2=convolution(f1, weights_conv2, biases_conv2);
%% conv3
f3=convolution(f2, weights_conv3, biases_conv3);

%% conv4
f4=convolution1(f3, weights_conv4, biases_conv4);
map=f4.^r;
lb= length(f4)/ofdm.K
%%  serial to parallel communication
    rcvd_sort_sgnl = dftmtx(ofdm.K);
    rcvd_sort_sgnl = rcvd_sort_sgnl(:,1:chan.L);
    chan.CoeffEst = zeros(ofdm.Nt,ofdm.Nr,chan.L,ofdm.Nb);
    
    for nb = 1 : ofdm.Nb
        for nr = 1 : ofdm.Nr
            chan.A = zeros(ofdm.PL,chan.L*ofdm.Nt);
            for nt = 1 : ofdm.Nt
                chan.A(:,(1:chan.L)+(nt-1)*chan.L) = diag(ofdm.dMod(ofdm.PPos,nb,nt))*rcvd_sort_sgnl(ofdm.PPos,:);
            end
            ChanEst = pinv(chan.A)*ofdm.fft(ofdm.PPos,nb,nr);
            for nt = 1 : ofdm.Nt
                chan.CoeffEst(nt,nr,:,nb) = ChanEst((1:chan.L)+(nt-1)*chan.L);
            end
        end        
    end
    chan.MSE_Simulation = var(chan.Coeff(:)-chan.CoeffEst(:));
    chan.MSE_Theory     = chan.sigma^2/ofdm.PL;
    if ofdm.ifDisplayResults
        disp(['MSE differentiate : ',num2str(chan.MSE_Theory)])
        disp(['MSE differentiate : ',num2str(chan.MSE_Simulation)])
    end
%% Demodulation
    if ofdm.ifDemodulateData == 1 && ofdm.DL > 0
        chan.CoeffEstFreq = zeros(ofdm.K,ofdm.Nt,ofdm.Nr,ofdm.Nb);
        for nb = 1 : ofdm.Nb
            for nr = 1 : ofdm.Nr
                for nt = 1 : ofdm.Nt
                    chan.CoeffEstFreq(:,nt,nr,nb) = rcvd_sort_sgnl*squeeze(chan.CoeffEst(nt,nr,:,nb));
                end
            end
        end
    end
        % demodulation
        ofdm.dDemod = zeros(ofdm.DL,ofdm.Nb,ofdm.Nt);
        for nb = 1 : ofdm.Nb
            for dl = 1 : ofdm.DL
                ofdm.dDemod(dl,nb,:) = pinv(reshape(chan.CoeffEstFreq(ofdm.DPos(dl),:,:,nb),ofdm.Nt,ofdm.Nr).')...
                                       *squeeze(ofdm.fft(ofdm.DPos(dl),nb,:));
            end
        end
        % QAM de-modulation
        ofdm.dEst = zeros(ofdm.DL,ofdm.Nb,ofdm.Nt);
        for nt = 1 : ofdm.Nt
            ofdm.dEst(:,:,nt) = qamdemod(1/ofdm.ModNorm * ofdm.dDemod(:,:,nt),ofdm.Mod);
        end
      figure,  
TSNR=[-5 0 5 10 15 20];
figure1_data

plot(TSNR,NMMSE,'y-*','Linewidth',2,'MarkerFaceColor','b','Markersize',5);
hold on
plot(TSNR,TDP,'r-*','Linewidth',2,'MarkerFaceColor','b','Markersize',5);
hold on
plot(TSNR,TSP,'b-*','Linewidth',2,'MarkerFaceColor','b','Markersize',5);
xlabel('Downlink SNR(dB)');
ylabel('NMSE(dB)');
legend('DP+Attention','DP','SP')


figure,
TSNR=[-5 0 5 10 15 20];
figure2_data
grid on
plot(TSNR,MSP,'r-*','Linewidth',2,'MarkerFaceColor','b','Markersize',5);
hold on
plot(TSNR,NMMSE,'y-*','Linewidth',2,'MarkerFaceColor','b','Markersize',5);
hold on
plot(TSNR,MDP,'b-*','Linewidth',2,'MarkerFaceColor','b','Markersize',5);
hold on
xlabel('Downlink SNR(dB)');
ylabel('NMSE(dB)');
legend('DP','DP+Attention','SP')

figure,
Axis=[-5 0 5 10 15 20];
figure_data
grid on
plot(Axis,NL,'b--','Linewidth',2,'MarkerFaceColor','b','Markersize',5);
hold on
plot(Axis,FNA,'r-*','Linewidth',2,'MarkerFaceColor','b','Markersize',5);
hold on
plot(Axis,PROPOSED,'y-s','Linewidth',2,'MarkerFaceColor','b','Markersize',5);
legend('FFT+NN+Attention','NN+LMMSE','Proposed');
xlabel('Downlinl SNR(dB)');
ylabel('NMSE(dB)');

