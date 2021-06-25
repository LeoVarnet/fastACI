function [str_stim,data_passation] = modulationACI_seeds_user(cfg,data_passation)
% function [str_stim,data_passation] = modulationACI_seeds_user(cfg,data_passation)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bLevel_norm_version = 3; % set to 1 for version 'as received'

i_current = data_passation.i_current;

n_stim = data_passation.n_stim(i_current);
istarget = (cfg.n_targets_sorted(n_stim) == 2);

if ~isfield(cfg,'bDebug')
    bDebug = 0;
else
    bDebug = cfg.bDebug;
end
fc   = cfg.fc;
fmod = cfg.fm;
dur  = cfg.stim_dur;
fs   = cfg.fs;
m_dB = data_passation.expvar(i_current);
m    = 10^(m_dB/20); % modulation index

% fprintf('%s: Generating noise...\n',upper(mfilename));
N_samples = round(cfg.stim_dur * fs);

%%% Fixing the seed
s_current = rng; % gets current seed
try
    seed_number = cfg.seeds_order(i_current); % needs to be a positive integer...
catch
    disp('')
end
rng(seed_number);
%%% 

noise = Generate_noise(N_samples,cfg.noise_type);

%%% Seed set back
rng(s_current);
%%% 

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
    case {2,3} % I suggest this calibration method
        
        rp    = ones(size(noise)); 
        rp(1:dur_ramp_samples)         = rampup(dur_ramp_samples);
        rp(end-dur_ramp_samples+1:end) = rampdown(dur_ramp_samples);
        
        dBFS       = cfg.dBFS;
        switch bLevel_norm_version
            case 2
                [stim_normalised,extra] = generate_stim(signal,noise,SNR,0,noise_type);
        
                lvl_before = rmsdb(stim_normalised);
                tuser_cal  = scaletodbspl(stim_normalised,SPL,dBFS);
                lvl_after  = rmsdb(tuser_cal);
        
                lvl_offset = (lvl_after-lvl_before);
                lvl_S = extra.lvl_S_dBFS+lvl_offset;
                lvl_N = extra.lvl_N_dBFS+lvl_offset;
                if bDebug == 1
                    fprintf('The exact level of the noise is %.1f dB (0 dB FS=%.1f)\n',lvl_N+dBFS,dBFS);
                    fprintf('The exact level of the pure tone (modulated or not) is %.1f dB\n',lvl_S+dBFS);
                end
            case 3
                noise  = scaletodbspl(noise ,SPL    ,dBFS); 
                signal = scaletodbspl(signal,SPL+SNR,dBFS); 
                extra.stim_N = noise;
                
                extra.stim_S = signal;
                tuser_cal = signal + noise;
        end
        
        % Applying the ramp:
        tuser_cal = rp.*tuser_cal;
        extra.stim_N = rp.*extra.stim_N;
        extra.stim_S = rp.*extra.stim_S;
        
        idx = find(cfg.sessionN_trials_validation>=i_current,1,'first');
        next_sessionN_trial_validation = cfg.sessionN_trials_validation(idx);
        
        if ~isempty(next_sessionN_trial_validation)
            if i_current == next_sessionN_trial_validation
                data_passation.waveforms_noise_alone(:,idx) = extra.stim_N;
                data_passation.waveforms_target_alone(:,idx) = extra.stim_S;
                data_passation.waveforms_trial(:,idx) = tuser_cal;
                data_passation.waveforms_i_current(:,idx) = i_current;
                data_passation.waveforms_n_stim(:,idx) = n_stim;
            end
        end       
        
end

try % Checking if the seed has changed since it was 'set back' and this line of the code
    bSame = Still_same_seed_status(s_current);
    if bSame == 0
        warning('The seed has changed since you generated the stimuli and the end of the script %s. PLEASE CHECK THIS!',mfilename)
    end
end
   
str_stim.tuser = tuser_cal;

switch bLevel_norm_version
    case 2
        str_stim.stim_noise_alone = gaindb(extra.stim_N,lvl_offset);
        str_stim.stim_tone_alone  = gaindb(extra.stim_S,lvl_offset);
    case 3
        str_stim.stim_noise_alone = extra.stim_N; 
        str_stim.stim_tone_alone  = extra.stim_S;
end

bCompare_with_stored = 0;

if bCompare_with_stored
    dir_where = [data_passation.dir_results cfg.Subject_ID filesep 'NoiseStims' filesep];
    if ~isfield(data_passation,'files_here')
        data_passation.files_here = Get_filenames(dir_where,'Noise*.wav');
        
        f2load = [dir_where 'Noises_all.mat'];
        if exist(f2load,'file')
            var = load(f2load);
            data_passation.noises_here = var.noise_all;
        end
    end
    figure;
    subplot(1,2,1)
    noise_here = str_stim.stim_noise_alone;
    plot(noise_here,'b-'); hold on;
    
    fname2load = [dir_where data_passation.files_here{i_current}];
    noise_stored_wav = audioread(fname2load);
    noise_stored_mat = data_passation.noises_here(:,i_current);
    % plot(noise_stored,'r--');
    
    ha = gca;
    
    subplot(1,2,2)
    plot(noise_here-noise_stored_wav,'k-'); hold on
    plot(noise_here-noise_stored_mat,'m--','LineWidth',2);
    legend('diff re. wav file','diff re. mat file')
    
    ha(end+1) = gca;
    linkaxes(ha,'x');
    
    warning('Comparing noises generated on the fly with stored versions... (set bCompare_with_stored to 0 for running the experiment normally)');
    
    pause();
end
disp('')