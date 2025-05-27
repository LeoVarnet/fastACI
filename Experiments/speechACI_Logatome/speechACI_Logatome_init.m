function cfg_inout = speechACI_Logatome_init(cfg_inout)
% function cfg_inout = speechACI_Logatome_init(cfg_inout)
%
% 1. Description:
%     This is the initialisation script for the experiment speechACI_Logatome.m.
%     This script should only be run once, to create the wave files.
%     
%     There are three types of validated noises (in cfg_inout.Condition) for
%     this experiment:
%           - Condition = 'white';
%           - Condition = 'bumpv1p2_10dB';
%           - Condition = 'sMPSv1p3';
%     Other possible conditions have not been validated for this experiment.
%
% The following processing is introduced here:
%     1. The speech samples are zero padded at the beginning and end by 
%        according to the duration of 'dur_ramp'
%     2. Up and down ramps are added to the noise samples. This is done by
%        adding a new section of noise at the beginning and end (defined by
%        dur_ramp), which is taken from a random section throughout the same noise,
%        i.e., (dur_ramp*fs) samples are randomly taken from the noise to be put
%        at the beginning, and a new selection of (dur_ramp*fs) samples is added
%        to the end (random selection done in Randomise_insig.m).
%     3. Finally a cosine ramp is added ONLY to the noise, and the noise is stored.
%     4. Note that all signals are scaled in level here already, depending
%        on the dBFS value (use dBFS = 93.6 dB for assuming HD 650 headphones)
%        and on lvl_target. lvl_target = 65 dB is a typical presentation level
%        in speech perception. The noises are, therefore, stored at an SNR=0 dB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
if nargin == 0
    help speechACI_Logatome_init;
    return;
end

dir_speech_orig = [fastACI_basepath 'Stimuli' filesep 'Logatome' filesep];

noise_type = cfg_inout.Condition;
switch noise_type
	case {'white','bumpv1p2_10dB','sMPSv1p3'}
        % Nothing to do...
    otherwise
        error('noise_type=%s not recognised. Run speechACI_LogatomeDebug to get more noise options',noise_type); % [cfg_inout.dir_logatome_src 'Noises' filesep];
end

dBFS       = cfg_inout.dBFS;
lvl_target = cfg_inout.SPL;
dur_ramp   = cfg_inout.dur_ramp; 

cfg_inout = Check_cfg_crea_dirs(cfg_inout); % updates the folders to the local computer if
                                            % crea comes from another computer
if ~exist(cfg_inout.dir_data_experiment,'dir')
    error('Please define an existing cfg_inout.dir_data_experiment (Not found: %s)',cfg_inout.dir_data_experiment);
end
dir_subject = [cfg_inout.dir_data_experiment cfg_inout.Subject_ID filesep];
if ~exist(dir_subject,'dir')
    mkdir(dir_subject);
end

%%% 1. Speech sounds:
dir_target = cfg_inout.dir_target;
 
if exist(dir_target,'dir')
    bGenerate_stimuli = 0;
else
    % The original speech sounds are zero padded:
    bGenerate_stimuli = 1;

    % If you are in this part of the code, 'dir_speech' does not exist
    mkdir(dir_target);
end
 
if bGenerate_stimuli
    
    files = il_get_waveforms_from_Cond_extra(dir_speech_orig,cfg_inout);
    
    for i = 1:length(files)
        [insig,fs]  = audioread([dir_speech_orig files{i}]);
                
        if i == 1
            if fs ~= cfg_inout.fs
                error('Sounds do not have the same sampling frequency as specified in the *._set.m file');
            end
        end
        lvls_before(i) = rmsdb(insig)+dBFS;
 
        insig = scaletodbspl(insig,lvl_target,dBFS); % level adjustment before the silence is added
                 
        lvls_after(i) = rmsdb(insig)+dBFS;
        
        sil = zeros(round(dur_ramp*fs),1);
        insig  = [sil; insig; sil];
        
        idx = strfind(files{i},'_eq');
        if ~isempty(idx)
            fprintf('\t%s: suffix ''eq'' found and being removed from the speech sample...\n',mfilename)
            files{i} = [files{i}(1:idx-1) '.wav'];
        end
        
        audiowrite([dir_target files{i}],insig,fs);
    end
else
    files = Get_filenames(dir_target,'*.wav');
    [insig,fs] = audioread([dir_target files{1}]); % reading one speech file
    if fs ~= cfg_inout.fs
        error('Sounds do not have the same sampling frequency as specified in the %s_set.m file');
    end
end
N_samples = length(insig);
 
%%% 2. Noise sounds:
dir_noise = cfg_inout.dir_noise;
if exist(dir_noise,'dir')
    % if dir is empty then it will generate the noises:
    bGenerate_stimuli = Check_if_dir_is_empty(dir_noise,'*.wav'); 
else
    bGenerate_stimuli = 1;
end

if bGenerate_stimuli == 1
    if ~strcmp(dir_noise(end),filesep) % if last character is \ or / (it should be the case always)
        dir_noise = [fileparts(dir_noise(1:end-1)) filesep];
    end
end
 
N_ramp = round(dur_ramp*fs); % ramp duration in samples
 
if bGenerate_stimuli
    [cfg_inout,s_current, bGenerate_stimuli] = Check_seeds_and_initialise(cfg_inout);
end

if bGenerate_stimuli == 0
    dir_noise = cfg_inout.dir_noise;
    
    files = Get_filenames(dir_noise,'*.wav');
else
    % If you are in this part of the code, 'dir_noise' does not exist:
    mkdir(dir_noise);
    
    s_current = rng; % gets current seed
end

ListStim = [];
for i = 1:cfg_inout.N
    if bGenerate_stimuli
        
        if i == 1
            % Generates the cosine ramp to be applied to the sounds:
            rp = ones([N_samples 1]); 
            rp(1:N_ramp)         = rampup(N_ramp);
            rp(end-N_ramp+1:end) = rampdown(N_ramp);
        end
        
        seed_number = cfg_inout.seeds_order(i);
        rng(seed_number); % insig = randn(N_samples,1); rng(seed_number)
        
        switch noise_type
            case {'bumpv1p2_10dB','sMPSv1p3'}
                insig = Generate_noise(N_samples,noise_type,fs);
            otherwise
                insig = Generate_noise(N_samples,noise_type,fs);
        end
        
        lvls_before_noise(i) = rmsdb(insig)+dBFS;
        insig = scaletodbspl(insig,lvl_target,dBFS);
        lvls_after_noise(i) = rmsdb(insig)+dBFS;
        insig = rp.*insig;
        
        fname_part1 = 'Noise'; % Bruit
        number = i;
 
        % Ensuring that the noise iteration number has four characters:
        stimnumber=num2str(number);
        if number<10000
            stimnumber=['0' stimnumber];
        end
        if number<1000
            stimnumber=['0' stimnumber];
        end
        if number<100
            stimnumber=['0' stimnumber];
        end
        if number<10
            stimnumber=['0' stimnumber];
        end
        % Ensuring 'End'
        fprintf(['Generating noise #' stimnumber '\n']);

        fname = [fname_part1 '_' stimnumber '.wav'];
    else
        if nargout >= 1
            fname = files{i};
        end
    end
    
    if nargout >= 1
        ListStim(i).name = fname;
    end
     
    if bGenerate_stimuli
        audiowrite([dir_noise fname],insig,fs);
    end
end
if bGenerate_stimuli
    try
        fsound = [fastACI_basepath 'Stimuli' filesep 'ready.wav'];
        [fsound,fs_sound] = audioread(fsound);
        sound(fsound,fs_sound);
    end
end

if isfield(cfg_inout,'bRove_level')
    
    if ~isfield(cfg_inout,'Rove_level')
        max_rove = cfg_inout.Rove_range; % dB
        cfg_inout.Rove_level = (rand(1,cfg_inout.N)-0.5)*2*max_rove;
    end

end

if nargout >= 1
    cfg_inout.ListStim = ListStim;
end

%%% Seed set back
if bGenerate_stimuli
    rng(s_current);
end
%%% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function files = il_get_waveforms_from_Cond_extra(dir_speech_orig,cfg_inout)

Speaker_ID = cfg_inout.Cond_extra_2;
Cond       = cfg_inout.Cond_extra_1;

files = Get_filenames(dir_speech_orig,[Speaker_ID '_*.wav']);

if length(files) < 2
    error('%s.m: Waveforms for condition %s (speaker %s) not found on disk',mfilename,Cond,Speaker_ID);
elseif length(files) > 2
    if length(Cond) == 4
        Cond1 = Cond(1:2);
        Cond2 = Cond(3:4);
    elseif length(Cond) == 5
        Cond1 = Cond(1:2);
        Cond2 = Cond(3:5);
    end
    files = Get_filenames(dir_speech_orig,[Speaker_ID '_*' Cond1 '*.wav']);
    
    extra_opts = [];
    extra_opts.bShow_header = 0;
    files(end+1) = Get_filenames(dir_speech_orig,[Speaker_ID '_*' Cond2 '*.wav'],extra_opts);
    
elseif length(files) == 2
    % Nothing to do: two files found
end
