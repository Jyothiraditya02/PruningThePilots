% channel parameters
    chan.SNR_dB  = 15;                  
    chan.L       = 6;                   
    
    % control parameters
    ofdm.ifDemodulateData = 1;          
    ofdm.ifDisplayResults = 1;       
    M=0.5;i=1;
    n=0.01;beta_k=0.2;
    a=0.02;r=10;
    
    
    if nargin > 2
%         error('Only two set as inputs thresholding')
    elseif nargin == 2
        % updating the set parameters
        freq_coeff = fieldnames(ofdmIn);
        for nS = 1:length(freq_coeff)
             ofdm.(freq_coeff{nS}) = ofdmIn.(freq_coeff{nS});
        end
        freq_coeff = fieldnames(chanIn);
        for nS = 1:length(freq_coeff)
             chan.(freq_coeff{nS}) = chanIn.(freq_coeff{nS});
        end
    elseif nargin == 1
        freq_coeff = fieldnames(ofdmIn);
        for nS = 1:length(freq_coeff)
             ofdm.(freq_coeff{nS}) = ofdmIn.(freq_coeff{nS});
        end
    end
    % dependent parameters
    ofdm.PPos    = 1:(ofdm.PSpace+1):ofdm.K;    
    ofdm.PL      = length(ofdm.PPos);           
    ofdm.DPos    = setxor(1:ofdm.K,ofdm.PPos);  
    ofdm.DL      = length(ofdm.DPos);           
    ofdm.BER     = 0;                          
    chan.sigma   = sqrt(10^(-0.1*chan.SNR_dB)); 