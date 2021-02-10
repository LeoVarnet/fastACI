function speechACI_varnet2015_preproc
% function speechACI_varnet2015_preproc
%
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