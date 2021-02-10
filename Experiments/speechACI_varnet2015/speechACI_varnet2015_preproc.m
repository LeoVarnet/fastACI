function speechACI_varnet2015_preproc
% function speechACI_varnet2015_preproc
%
% This should only be run once, to create the wave files.
%
% This script assumes that the sounds 'Alda.wav', 'Alga.wav', Arda.wav','Arga.wav'
%     are located in the folder: '/home/alejandro/Documents/Databases/data/fastACI/speechACI/speech-samples_orig/Alda.wav'
%
% This script assumes that the noises are located in:
%     /home/alejandro/Documents/Databases/data/fastACI/speechACI/NoiseStim_S11_orig/
%     i.e., temporarily I pick up the noises of participant S11 only, it will be
%     solved in the future to pick up the wavefiles (or the seeds) of the correct
%     and customised participant.
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
 
% cfg = [];
% cfg_out = cfg_in; % copying input to output struct

dBFS = 100;
lvl_target = 65;

dir_main = '/home/alejandro/Documents/Databases/data/fastACI/speechACI/';
dir_speech     = [dir_main 'speech-samples_orig' filesep];
dir_speech_out = [dir_main 'speech-samples'      filesep];

dur_ramp = 75e-3;
%%%
if isdir(dir_speech_out)
    bGenerate_sounds = 0;
else
    bGenerate_sounds = 1;
end
    
if bGenerate_sounds
    files = Get_filenames(dir_speech,'*.wav');
    for i = 1:length(files)
        [insig,fs] = audioread([dir_speech files{i}]);

        lvls_before(i) = rmsdb(insig)+dBFS;

        insig = setdbspl(insig,lvl_target,dBFS);
        
        lvls_after(i) = rmsdb(insig)+dBFS;
        
        if i == 1
            sil = zeros(round(dur_ramp*fs),1);
            mkdir(dir_speech_out);
        end
        audiowrite([dir_speech_out files{i}],[sil; insig; sil],fs);
    end
end
%%%

dir_noise     = [dir_main 'NoiseStim_S11_orig' filesep];
dir_noise_out = [dir_main 'NoiseStim_S11'      filesep];

files = Get_filenames(dir_noise,'*.wav');

if isdir(dir_noise_out)
    bGenerate_sounds = 0;
else
    bGenerate_sounds = 1;
end
if bGenerate_sounds
    for i = 1:length(files)
        [insig,fs] = audioread([dir_noise files{i}]);

        if i == 1
            N = round(dur_ramp*fs); % ramp duration in samples
        end
        noise_rampup = Randomise_insig(insig,N);
        noise_rampdn = Randomise_insig(insig,N);
        
        insig = [noise_rampup; insig; noise_rampdn];
        lvls_before_noise(i) = rmsdb(insig)+dBFS;

        insig = setdbspl(insig,lvl_target,dBFS);

        lvls_after_noise(i) = rmsdb(insig)+dBFS;
        
        if i == 1
            rp    = ones(size(insig)); 
            rp(1:N)         = rampup(N);
            rp(end-N+1:end) = rampdown(N);
            
            mkdir(dir_noise_out);
        end
        insig = rp.*insig;
        
        fname_part1 = strsplit(files{i}(1:end-4),'_');
        number = str2double(fname_part1{end});
        
        fname_part1 = fname_part1{1};
        
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
                
        audiowrite([dir_noise_out fname_part1 '_' stimnumber '.wav'],insig,fs);
    end
end

fprintf('Make sure that you set the correct paths in speech_ACI_varnet2015_set.m, if you don''t relocate the files, use:\n');
fprintf('\ndir_main=%s\n',dir_main);
fprintf('\ndir_speech=[dir_main ''speech_samples'' filesep];\n');
fprintf('\ndir_speech=[dir_main ''NoiseStim_S11''  filesep];\n');