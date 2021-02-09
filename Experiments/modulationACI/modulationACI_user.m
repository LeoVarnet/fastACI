function str_stim = modulationACI_user(str_inout,cfg)
% function str_stim = modulationACI_user(str_inout,cfg)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

istarget = str_inout.istarget;
if ~isfield(str_inout,'bLevel_norm_version')
    str_inout.bLevel_norm_version = 2;
end
bLevel_norm_version = str_inout.bLevel_norm_version;

if ~isfield(str_inout,'filename')
    str_inout.filename = '';
    bLoad = 0;
else
    bLoad = 1;
end

fc   = cfg.fc;
fmod = cfg.fm;
dur  = cfg.stim_dur;
fs   = cfg.fs;
m_dB = str_inout.m;
switch bLevel_norm_version
    case 1
        % Unintended denominator in the conversion from modulation depth to
        %    modulation index:
        m = 10^(m_dB/10); 
    case 2
        m = 10^(m_dB/20); % modulation index
end

if bLoad
    file2load = [cfg.dir_stim cfg.folder_name filesep str_inout.filename];
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
end

% ---
signal = create_AM(fc, fmod, m*istarget, dur, fs)';

% ADD SILENCE FOR THE MODEL:
signal = [zeros(length(noise)-length(signal),1); signal];

% create stim
SNR = cfg.SNR;
noise_type = cfg.noise_type;
SPL = cfg.SPL;
dur_ramp_samples = cfg.fs*cfg.fadein_s;

switch bLevel_norm_version
    case 1 % default in ENS 'as received'
        [stim_normalised,extra] = generate_stim(signal,noise,SNR,dur_ramp_samples,noise_type);
        [tuser_cal,dBFS,gain_factor] = dBlvl(stim_normalised,SPL);
    case 2 % I suggest this calibration method
        [stim_normalised,extra] = generate_stim(signal,noise,SNR,0,noise_type);

        rp    = ones(size(noise)); 
        rp(1:dur_ramp_samples)         = rampup(dur_ramp_samples);
        rp(end-dur_ramp_samples+1:end) = rampdown(dur_ramp_samples);
        
        dBFS       = cfg.dBFS;
        lvl_before = rmsdb(stim_normalised);
        tuser_cal  = setdbspl(stim_normalised,SPL,dBFS);
        lvl_after  = rmsdb(tuser_cal);
        
        lvl_offset = (lvl_after-lvl_before);
        lvl_S = extra.lvl_S_dBFS+lvl_offset;
        lvl_N = extra.lvl_N_dBFS+lvl_offset;
        
        % Applying the ramp:
        tuser_cal = rp.*tuser_cal;
        extra.stim_N = rp.*extra.stim_N;
        extra.stim_S = rp.*extra.stim_S;
        
        if cfg.bDebug == 1
            fprintf('The exact level of the noise is %.1f dB (0 dB FS=%.1f)\n',lvl_N+dBFS,dBFS);
            fprintf('The exact level of the pure tone (modulated or not) is %.1f dB\n',lvl_S+dBFS);
        end
end

str_stim.tuser = tuser_cal;
str_stim.stim_noise_alone = gaindb(extra.stim_N,lvl_offset);
str_stim.stim_tone_alone  = gaindb(extra.stim_S,lvl_offset);