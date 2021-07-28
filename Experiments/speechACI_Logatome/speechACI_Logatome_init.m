function cfg_inout = speechACI_Logatome_init(cfg_inout)
% function cfg_inout = speechACI_Logatome_init(cfg_inout)
%
% This should only be run once, to create the wave files.
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
 
dir_speech_orig = [fastACI_paths('dir_fastACI') 'Stimuli' filesep 'Logatome' filesep];
dir_noise_spectrum = []; warning('Temporal'); % [cfg_inout.dir_logatome_src 'Noises' filesep];

dBFS       = cfg_inout.dBFS;
lvl_target = cfg_inout.SPL;
dur_ramp   = cfg_inout.dur_ramp; 

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
    
    Speaker_ID = cfg_inout.Cond_extra_2;
    files = Get_filenames(dir_speech_orig,[Speaker_ID '*.wav']);
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
    bGenerate_stimuli = 0;
else
    bGenerate_stimuli = 1;
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
            noise_type = cfg_inout.Condition;
            switch noise_type
                case 'SSN' % apply the SSN
                    % See g20210401_generate_LTASS.m
                    
                    error('We are working for you: the option SSN will be enabled back soon...')
                    
                    if ~exist(dir_noise_spectrum,'dir')
                        dir_noise_spectrum = uigetdir(pwd,'Select the directory where ''Logatome-average-power-speaker-S46M_FR.mat'' is located');
                        dir_noise_spectrum = [dir_noise_spectrum filesep]; % file separator is the last character
                    end

                    % Make sure that S46 is the speaker used in the chosen logatomes
                    file = [dir_noise_spectrum 'Logatome-average-power-speaker-S46M_FR.mat']; % Logatome-average-power-speaker-S46M_FR.mat
                    var = load(file,'f','Pxx','fs');
                    if fs ~= var.fs
                        error('Wrong sampling frequency to generate the SSN')
                    end
                    
                    [~,~,~,b_fir] = Create_LTASS_noise(var.f,var.Pxx,fs);
            end
            
            % Generates the cosine ramp to be applied to the sounds:
            rp = ones([N_samples 1]); 
            rp(1:N_ramp)         = rampup(N_ramp);
            rp(end-N_ramp+1:end) = rampdown(N_ramp);
        end
        
        seed_number = cfg_inout.seeds_order(i);
        rng(seed_number); % insig = randn(N_samples,1); rng(seed_number)
        
        switch noise_type
            case 'SSN' % apply the SSN
                insig = Generate_noise(N_samples,'white');
                insig = filter(b_fir,1,insig);
            case 'bump'
                % Nothing to do
            otherwise
                insig = Generate_noise(N_samples,noise_type,fs);
        end
        
        switch noise_type
            case 'bump'
                % Nothing to do: Noises will be generated during the participant's session
            otherwise
                lvls_before_noise(i) = rmsdb(insig)+dBFS;
                insig = scaletodbspl(insig,lvl_target,dBFS);
                lvls_after_noise(i) = rmsdb(insig)+dBFS;
                insig = rp.*insig;
        end
        
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
        fname = files{i};
    end
    
    ListStim(i).name = fname;
     
    if bGenerate_stimuli
        switch noise_type
            case 'bump'
                % Nothing to do: Noises will be generated during the participant's session
            otherwise
                audiowrite([dir_noise fname],insig,fs);
        end
    end
end

if isfield(cfg_inout,'bRove_level')
    
    if ~isfield(cfg_inout,'Rove_level')
        max_rove = cfg_inout.Rove_range; % dB
        cfg_inout.Rove_level = (rand(1,cfg_inout.N)-0.5)*2*max_rove;
    end

end

cfg_inout.ListStim = ListStim;

%%% Seed set back
if bGenerate_stimuli
    rng(s_current);
end
%%% 