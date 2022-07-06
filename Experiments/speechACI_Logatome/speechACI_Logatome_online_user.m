function [str_stim,data_passation] = speechACI_Logatome_online_user(cfg,data_passation)
% function [str_stim,data_passation] = speechACI_Logatome_online_user(cfg,data_passation)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

str_stim = [];

i_current = data_passation.i_current;
SNR       = data_passation.expvar(i_current);

n_stim    = data_passation.n_stim(i_current);
n_signal  = cfg.n_targets_sorted(n_stim);

if isfield(cfg,'Rove_level')
    % Rove_range has to be adjusted in the _set.m script.
    data_passation.Rove_level(i_current) = cfg.Rove_level(n_stim);
    presentation_gain = 10^(data_passation.Rove_level(i_current)/20);
else
    presentation_gain = 1;
end

[signal,fs] = audioread([cfg.dir_target cfg.filename_target{n_signal}]); % will load one of the four utterances

dBFS = cfg.dBFS;
lvl_target = cfg.SPL;
dur_ramp   = cfg.dur_ramp; 
%%% Fixing the seed
if isfield(cfg,'seeds_order')
    bStore_seeds = 1;
    s_current = rng; % gets current seed
    seed_number = cfg.seeds_order(i_current); % needs to be a positive integer...
    rng(seed_number);
else
    bStore_seeds = 0;
end

N_samples = length(signal);
N_ramp = round(dur_ramp*fs); % ramp duration in samples

rp = ones([N_samples 1]); 
rp(1:N_ramp)         = rampup(N_ramp);
rp(end-N_ramp+1:end) = rampdown(N_ramp);

switch cfg.noise_type
    case 'SSN' % apply the SSN
        error('obtain the b_fir coefficients (see init file)')
        noise = Generate_noise(N_samples,'white');
        noise = filter(b_fir,1,noise);
    case {'bumpv1p2_10dB','sMPSv1p3'}
        noise = Generate_noise(N_samples,cfg.noise_type,fs);
    case {'bumpv1p1','bumpv1p1_30dB','bumpv1p2','bumpv1p2_5dB', ...
          'sMPSv1p1','sMPSv1p2'}
        % These are previous versions of the noise algorithms (that 
        %     we decided to exclude
        noise = Generate_noise_debug(N_samples,cfg.noise_type,fs);    
    otherwise
        noise = Generate_noise(N_samples,cfg.noise_type,fs);
end
        
% lvls_before_noise(i) = rmsdb(noise)+dBFS;
noise = scaletodbspl(noise,lvl_target,dBFS);
% lvls_after_noise(i) = rmsdb(noise)+dBFS;
noise = rp.*noise;
        
if bStore_seeds
    %%% Seed set back
    rng(s_current);
    %%% 
end
 
bSpeech_level_variable = 1;
bNoise_level_variable = ~bSpeech_level_variable;
 
if bSpeech_level_variable
    gain_snr = 10^(SNR/20);
    signal = gain_snr * signal;
end
if bNoise_level_variable
    error('Not validated yet...')
    % gain_snr = 10^(-SNR/20);
    % noise = gain_snr * noise;
end

tuser_cal = noise+signal;
 
str_stim.tuser            = presentation_gain*tuser_cal;
str_stim.stim_tone_alone  = signal;
