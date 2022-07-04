function [outsig, mfc, step] = osses2022_modfilterbank(insig,fs,fc,varargin)
%osses2022_modfilterbank modulation filterbank used by King et al. 2019
%
%   Usage: [outsig, mfc, step] = king2019_modfilterbank(insig,fs,varargin)
% 
%   Input parameters:
%     insig: Input signal(s)
%     fs     : Sampling rate in Hz
%     fc     : Centre frequencies of the input signals
%     mflow  : minimum modulation centre frequency in Hz
%     mfhigh : maximum modulation centre frequency in Hz
%     N      : Number of logarithmically-spaced (between fmin and fmax) 
%              modulation filters (default N = 10)
%     Qfactor: Quality factor of the filters (default Qfactor = 1).
%
%   Output parameters:
%     outsig: Modulation filtered signals 
%     mfc   : Centre frequencies of the modulation filters
%     step  : Contains some intermediate outputs.
%
%   #Author: Leo Varnet and Andrew King (2020)
%   #Author: Alejandro Osses (2020) Original implementation for the AMT
%   #Author: Clara Hollomey (2021) Adapted for AMT
%   #Author: Piotr Majdak (2021) Further adaptations to AMT 1.0

definput.keyvals.mfc=[];
definput.import = {'modfilterbank_local'};
definput.importdefaults = {'mfb_osses2022a'};
[flags,kv]=ltfatarghelper({},definput,varargin);

nfreqchannels=length(fc);
if nfreqchannels == 0
    % This can be the case if the filter bank stage was by-passed, then:
    nfreqchannels = size(insig,2);
    fc = fs/2; % frequency that is high enough to ensure that all modfilters are assessed
end
outsig=cell(nfreqchannels,1);

% first order modulation Butterworth lowpass filter with a cut-off
% frequency of 150 Hz. This is to remove all modulation frequencies
% above 150 Hz. The motivation behind this filter can be found in kohlrausch2000
if flags.do_LP_150_Hz || flags.do_LP_150_Hz_att % modbank_LPfilter
    [b_lp_150_Hz,a_lp_150_Hz] = butter(1,150/(fs/2));
    if flags.do_LP_150_Hz
        % In this case, the filtering is directly applied
        insig = filter(b_lp_150_Hz,a_lp_150_Hz,insig);
    end
    %%%
    if flags.do_LP_150_Hz_att
        % In this case the frequency response is assessed
        K =  8192; % arbitrary number
        f = (0:K-1)/K * fs/2;
        hLPF    = freqz(b_lp_150_Hz,a_lp_150_Hz, K); % frequency response of the 150 Hz filter
        hLPF_dB = 20*log10(abs(hLPF));
    end
end
mfc_upper_limit_max = kv.mfc_upper_limit_max; % Hz
if flags.do_mfc_upper_limit
    umf = min(fc.*0.25, mfc_upper_limit_max);
end

if flags.do_no_mfc_upper_limit
    umf = mfc_upper_limit_max*ones(1,nfreqchannels);
end

% second order modulation Butterworth lowpass filter with a cut-off frequency 
% of 2.5 Hz.
[b_lowpass,a_lowpass] = butter(2,2.5/(fs/2));

% Parameters modulation filter:
mflow           = kv.mflow; 
mfhigh          = kv.mfhigh; 
modbank_Nmod    = kv.modbank_Nmod;
Q_mfb           = kv.Q_mfb;
    
% -------------------------------------------------------------------------
% -- 1. modulation_filterbank.m
if isempty(modbank_Nmod)
    step_mfc = (sqrt(4*Q_mfb^2+1)+1)/(sqrt(4*Q_mfb^2+1)-1);
    logfmc = log(mflow):log(step_mfc):log(mfhigh);
    
    modbank_Nmod = length(logfmc);
else
    logfmc = linspace(log(mflow), log(mfhigh), modbank_Nmod); %log(fmin):log((sqrt(4*Qfactor^2+1)+1)/(sqrt(4*Qfactor^2+1)-1)):log(fmax);
end
mfc = exp(logfmc);

for ichan = 1:modbank_Nmod
    flim(ichan,:) = mfc(ichan)*sqrt(4+1/Q_mfb^2)/2 +  [-1 +1]*mfc(ichan)/Q_mfb/2;
    [b(ichan,:),a(ichan,:)] = butter(2,2*[flim(ichan,:)]/fs); 
end
mfc = [2.5 mfc];
b = [b_lowpass 0 0; b];
a = [a_lowpass 0 0; a];
Nchannels = modbank_Nmod+1;

flim = [0 2.5; flim];

if flags.do_LP_150_Hz_att
    mfc_gains = interp1(f(:)',hLPF_dB(:)',mfc);
end
% -------------------------------------------------------------------------
% 2. apply_filterbank.m
if size(a,1)~= Nchannels
    error('Number of lines in ''a'' and ''b'' must be equal');
end

for j = 1:nfreqchannels
    for i=1:Nchannels
        if mfc(i)<=umf(j)
            outsig{j}(:,i) = filter(b(i,:),a(i,:),insig(:,j));
            
            if flags.do_LP_150_Hz_att
                outsig{j}(:,i) = gaindb( outsig{j}(:,i), mfc_gains(i));
            end
        end
    end
end

if nargout >= 3 % Intermediate outputs
    step.mfb_a = a;
    step.mfb_b = b;
    step.a_lowpass = a_lowpass;
    step.b_lowpass = b_lowpass;
    step.flim = flim;
    step.flim_description = 'Limits of the modulation filters (Hz)';
    step.fs_design = fs;
    
    step.E_mod = outsig;
    step.fmc = mfc;
end
    
% -------------------------------------------------------------------------
% 3. Phase insensitivity
% ---
if flags.do_phase_insens_hilbert
    % amt_disp('  Phase insensitivity (hilbert)',flags.disp);
    
    if flags.do_att_factor
        Factor = 1/sqrt(2); % This is according to Jepsen et al. (2008), page 426
    end
    if flags.do_no_att_factor
        Factor = 1;
    end
    
    phase_insens_cut = kv.phase_insens_cut; % Hz
    
    for j=1:nfreqchannels
        Nmod_here = size(outsig{j},2);
        for i=1:Nmod_here
            if mfc(i)>phase_insens_cut
                outsig_env = il_hilbert_extraction(squeeze(outsig{j}(:,i)));
                outsig{j}(:,i) = outsig_env * Factor;
            end
       end
    end
    
    if nargout>=3
        step.E_phase_ins = outsig;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function outsig_env = il_hilbert_extraction(D)
% Author: Leo Varnet 2016 - last modified 11/10/2018

if isreal(D)
    hilbert_responses = hilbert(D);
else
    hilbert_responses = D;
end

outsig_env = abs(hilbert_responses);