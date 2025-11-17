function [str_stim,data_passation] = localisationILD_user(cfg,data_passation)
% function [str_stim,data_passation] = localisationILD_user(cfg,data_passation)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i_current = data_passation.i_current;

n_stim = data_passation.n_stim(i_current);
isleft  = (cfg.n_targets_sorted(n_stim) == 1);
isright = ~isleft;

if ~isfield(cfg,'bDebug')
    bDebug = 0;
else
    bDebug = cfg.bDebug;
end

fs   = cfg.fs;
ILD = data_passation.expvar(i_current); % Frequency deviation delta f
% fdev = 10^(m_dB/20); % modulation index

if ~isfield(cfg,'ListStim')
    bLoad = 0; % probably the init file is being called
else
    bLoad = 1;
    filename = cfg.ListStim(n_stim).name;
end

% ---
% create stim
if isfield(cfg,'noise_type')
    noise_type = cfg.noise_type;
else
    noise_type = 'white';
end
SPL = cfg.SPL;
dur_ramp_samples = cfg.fs*cfg.fadein_s;

if bLoad
    file2load = [cfg.dir_noise filename];
    noise = audioread(file2load);
else
    fprintf('%s: Generating noise...\n',upper(mfilename));
    N_samples = round(cfg.stim_dur * fs);
    switch cfg.noise_type
        case 'white'
            noise=randn(N_samples,1);
        case 'pink'
            error('Not validated yet...')
            noise=pinknoise(N_samples)';
        otherwise
            error('%s: Unknown type of noise. Possible options are ''pink'' or ''white''',upper(mfilename))
    end
    noise = scaletodbspl(noise,SPL,cfg.dBFS);
end

gain_ild = 10^(ILD/20);
%%% bLevel_norm_version, option '2' in the old versions of this script:
if isleft
    stim_normalised = [gain_ild*noise noise];
end
if isright
    stim_normalised = [noise gain_ild*noise];
end

% lvl_offset = (lvl_after-lvl_before);

% Applying the ramp:

str_stim.tuser = stim_normalised;
str_stim.stim_noise_alone = noise;
% str_stim.stim_tone_alone  = gaindb(extra.stim_S,lvl_offset);