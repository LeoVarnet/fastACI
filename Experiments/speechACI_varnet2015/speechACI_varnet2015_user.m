function str_stim = speechACI_varnet2015_user(str_inout,cfg,data_passation)
% function str_stim = speechACI_varnet2015_user(str_inout,cfg,data_passation)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% In the paper: 0 for ‘da’ and 1 for ‘ga’
% Here: 1 = 'da', 2 = 'ga'

istarget = str_inout.istarget;
 
% fc   = cfg.fc;
% fmod = cfg.fm;
% dur  = cfg.stim_dur;
SNR = str_inout.expvar;
% switch bLevel_norm_version
%     case 1
%         % Unintended denominator in the conversion from modulation depth to
%         %    modulation index:
%         m = 10^(m_dB/10); 
%     case 2
%         m = 10^(m_dB/20); % modulation index
% end

switch istarget
    case 0
        [signal,fs] = audioread([cfg.dir_speech 'Arda.wav']);
    case 1
        [signal,fs] = audioread([cfg.dir_speech 'Arga.wav']);
end

if cfg.warmup == 1
    % % a random noise is picked up:
    % n_noise = round( (cfg.N-1)*random('unif',0,1) )+1; 
    noise = audioread([cfg.dir_noise str_inout.filename]);
else
    i_current = data_passation.i_current;
    n_noise = cfg.randorder_idxs(i_current);
    noise = audioread([cfg.dir_noise cfg.ListStim(n_noise).name]);
end

bSpeech_level_variable = 1;
bNoise_level_variable = ~bSpeech_level_variable;

if bSpeech_level_variable
    gain_snr = 10^(SNR/20);
    signal = gain_snr * signal;
end
if bNoise_level_variable
    error('Not validated yet...')
    gain_snr = 10^(-SNR/20);
    noise = gain_snr * noise;
end
tuser_cal = noise+signal;

str_stim.tuser = tuser_cal;

str_stim.stim_tone_alone  = signal;
% str_stim.stim_noise_alone = signal2;