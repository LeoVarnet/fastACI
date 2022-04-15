function [str_stim,data_passation] = speechACI_Logatome_user(cfg,data_passation)
% function [str_stim,data_passation] = speechACI_Logatome_user(cfg,data_passation)
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
fname_noise = [cfg.dir_noise cfg.ListStim(n_stim).name];
switch cfg.noise_type
    case 'bump'
        %%% Fixing the seed
        if isfield(cfg,'seeds_order')
            bStore_seeds = 1;
            s_current = rng; % gets current seed
            seed_number = cfg.seeds_order(i_current); % needs to be a positive integer...
            rng(seed_number);
        else
            bStore_seeds = 0;
        end
        
        %%%         
        % if ~exist(fname_noise,'file') || cfg.warmup == 0
        % Needs creation
        sigma_t = 0.02; % temporal width of the bumps (in s)
        sigma_ERB = 0.5; % spectral width of the bumps (in ERB)
        A = 30; % amplitude of the bumps (in dB)

        lvl_bump_noise = 60;
        % bump noise generation
        noise = bumpnoisegen(length(signal), fs, Nb, sigma_t, sigma_ERB, A, lvl_bump_noise, cfg.dBFS);
        if exist(fname_noise,'file')
            % disp('The sound exists on disk and will be overwritten, press ctrl+C to stop this action (you have 5 seconds to react)');
            % pause(5);
        end
        audiowrite(fname_noise,noise,fs);

        data_passation.bump_numbers(i_current) = Nb;
        
        if bStore_seeds
            %%% Seed set back
            rng(s_current);
            %%% 
        end
        
    otherwise
        noise = audioread(fname_noise);
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
