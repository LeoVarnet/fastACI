function [metric,description,extra] = Get_envelope_metric(insig,fs,type)
% function [metric,description,extra] = Get_envelope_metric(insig,fs,type)
%
% 1. Description:
%       Assesses one of the measures of envelope fluctuations.
%
%       'V': ratio between the standard deviation and the mean of the envelope
%            for an AM sinusoid V is always 3 dB below the modulation depth,
%            and it is -inf for a flat Gaussian noise.
%
%       'W': fourth-moment of the waveform 'insig'. Kohlrausch et al. (1997,
%            Eq. 1) attribute this metric to Hartmann and Pumplin and ranges
%            between 1.5 and 3.
%
%       'varnet2021_env': Gammatone filtered envelope, followed by a HWR and
%            a low pass filter, using a 30-Hz first-order butterworth filter.
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses
% Created on    : 13/09/2021
% Last update on: 13/09/2021 
% Last use on   : 13/09/2021 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

extra = [];

switch type
    case {'kohlrausch2021_env'}
        fs_env = 1000;
        dBFS = 100; % just a reference
        
        % windowtype = 'rectangular';
        windowtype = 'hanning';
        yenv_tmp    = abs(hilbert(insig));
    
        yenv = resample(yenv_tmp,fs_env,fs);    
        K = size(yenv,1);
        
        [~,yenv_dB,f_env] = freqfft2(yenv,K,fs_env,windowtype,dBFS);
        
        metric = yenv_dB;
        extra.fs_env = fs_env;
        extra.f_env = f_env;
        description = 'Envelope spectrum (hilbert envelope + FFT) as used by Kohlrausch et al 2021';
        
    case {'kohlrausch2021_env_noDC'}
        fs_env = 1000;
        dBFS = 100; % just a reference
        
        % windowtype = 'rectangular';
        windowtype = 'hanning';
        yenv_tmp    = abs(hilbert(insig));
    
        yenv = resample(yenv_tmp,fs_env,fs);    
        K = size(yenv,1);
        
        [~,yenv_dB,f_env] = freqfft2(yenv-mean(yenv),K,fs_env,windowtype,dBFS);
        
        metric = yenv_dB;
        extra.fs_env = fs_env;
        extra.f_env = f_env;
        description = 'Envelope spectrum (hilbert envelope + FFT) as used by Kohlrausch et al 2021 but excluding the DC';
        
    case {'varnet2021_env'}
        basef = 1000;
        outsig = auditoryfilterbank(insig,fs,'basef',basef,'flow',basef,'fhigh',basef);
        outsig = ihcenvelope(outsig,fs,'ihc_king2019');
        
        % Extra LPF filter at 30 Hz:
        cutofffreq=30;
        ihc_filter_order = 2;
        [b, a] = butter(1, cutofffreq*2/fs);
        for ii=1:ihc_filter_order
            outsig = filter(b,a, outsig);
        end
        metric = outsig;
        description = 'Gammatone + HWR + LPF at 30 Hz, as used by Varnet and Lorenzi (2021)';
        
        extra.basef = basef;
        
    case {'varnet2017_AMi'}
        metric = varnet2017_AMi(insig,fs);
        description = 'AMi as described by Varnet et al. (2017)';
        
    case {'V','v'}
        yenv = abs(hilbert(insig));
        extra.yenv = yenv;
        metric = 20*log10(std(yenv)/mean(yenv));
        description = 'V as defined by Kohlrausch et al. (1997)';
        
    case {'W','w'}
        metric = mean(insig.^4)/( mean(insig.^2).^2 );
        description = 'W as defined by Kohlrausch et al. (1997)';
        
    otherwise
        error('Metric not found...')
end
