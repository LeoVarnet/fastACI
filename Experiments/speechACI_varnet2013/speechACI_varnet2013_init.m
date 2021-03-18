function cfg_inout = speechACI_varnet2013_init(cfg_inout)
% function cfg_inout = speechACI_varnet2013_init(cfg_inout)
%
% This script assumes that the sounds 'Aba.wav', 'Ada.wav'
%     are located in the folder: '/home/alejandro/Documents/Databases/data/fastACI/speechACI/speech-samples_orig/Alda.wav'
%
% This script assumes that the noises are located in:
%     /home/alejandro/Documents/Databases/data/fastACI/speechACI/NoiseStim_S11_orig/
%     i.e., temporarily I pick up the noises of participant S11 only, it will be
%     solved in the future to pick up the wavefiles (or the seeds) of the correct
%     and customised participant.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%%% Values that are defined in speechACI_varnet2013_set.m:
%     These parameters are used only if the sounds need to be re-stored
dBFS       = cfg_inout.dBFS;
lvl_target = cfg_inout.SPL;
dur_ramp   = cfg_inout.dur_ramp; 
fs         = cfg_inout.fs;    
%%% 1. Speech sounds:
dir_speech = cfg_inout.dir_speech;
if exist(dir_speech,'dir')
    bGenerate_stimuli = 0;
else
    % The original speech sounds are zero padded:
    bGenerate_stimuli = 1;
    
    if strcmp(dir_speech(end),filesep) % if last character is \ or / (it should be the case always)
        dir_main = fileparts(dir_speech(1:end-1));
        dir_main = [fileparts(dir_main) filesep];
    end
    
    % dir_speech_orig should contain the original speech samples:
    dir_speech_orig = [dir_main 'speech-samples_orig' filesep];
    
    if ~exist(dir_speech_orig,'dir')
        error('The original speech samples (aba and ada) were not found in your computer...');
    end
        
    % If you are in this part of the code, 'dir_speech' does not exist
    mkdir(dir_speech);
end

if bGenerate_stimuli
    files = Get_filenames(dir_speech_orig,'*.wav');
    for i = 1:length(files)
        [insig,fs] = audioread([dir_speech_orig files{i}]);

        if i == 1
            if fs ~= cfg_inout.fs;
                error('Sounds do not have the same sampling frequency as specified in the %s_set.m file',experiment);
            end
        end
        insig = setdbspl(insig,lvl_target,dBFS);
        
        if i == 1
            sil = zeros(round(dur_ramp*fs),1);
        end
        insig = [sil; insig; sil];
        audiowrite([dir_speech files{i}],insig,fs);
    end
else
    files = Get_filenames(dir_speech,'*.wav');
    [insig,fs] = audioread([dir_speech files{1}]); % reading one speech file
    if fs ~= cfg_inout.fs
        error('Sounds do not have the same sampling frequency as specified in the %s_set.m file',experiment);
    end
end
N_samples = length(insig);

%%% 2. Noise sounds:
dir_noise = cfg_inout.dir_noise;
if exist(dir_noise,'dir')
    bGenerate_stimuli = 0;
else
    bGenerate_stimuli = 1;
    if strcmp(dir_noise(end),filesep) % if last character is \ or / (it should be the case always)
        dir_main = [fileparts(dir_noise(1:end-1)) filesep];
    end
    
    % If you are in this part of the code, 'dir_noise' does not exist:
    mkdir(dir_noise);
end

N_ramp = round(dur_ramp*fs); % ramp duration in samples

if bGenerate_stimuli
    % Nothing to do...
else
    files = Get_filenames(dir_noise,'*.wav');
end

ListStim = [];
for i = 1:cfg_inout.N
    if bGenerate_stimuli
        
        insig = randn(N_samples,1); % This N_samples already includes the ramp times
        
        insig = setdbspl(insig,lvl_target,dBFS);
        
        if i == 1
            % Generates the cosine ramp:
            rp = ones(size(insig)); 
            rp(1:N_ramp)         = rampup(N_ramp);
            rp(end-N_ramp+1:end) = rampdown(N_ramp);
        end
        
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

        fname = [fname_part1 '_' stimnumber '.wav'];
    else
        fname = files{i};
    end
        
    ListStim(i).name = fname;
    
    if bGenerate_stimuli
        audiowrite([dir_noise fname],insig,fs);
    end
end
cfg_inout.ListStim = ListStim;

if nargout == 0
%     fprintf('List of directories where the sounds should be:\n');
%     fprintf('\ndir_main=%s\n',dir_main);
%     fprintf('\ndir_speech=[dir_main ''speech_samples'' filesep];\n');
%     fprintf('\ndir_speech=[dir_main ''NoiseStim_S11''  filesep];\n');
end