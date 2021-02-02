function Script1_Initialisation_EN
% function Script1_Initialisation_EN
%
% 
% Changes by AO:
%   - cfg_crea.color_noise changed by cfg_crea.noise_type
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc, close all

%% Parameters to create the targets:
cfg_crea.fs              = 48000; % Hz, sampling frequency
cfg_crea.stim_dur        = 1.5;   % s,  stimulus suration (s)
cfg_crea.N_noise         = 1500;  % number of stimuli / condition
cfg_crea.N_signal        = 2;     % Number of conditions
cfg_crea.noise_type      = 'white';

[path,name,ext]=fileparts(which(mfilename)); % path will be the folder where this file is located...
cfg_crea.path            = [path filesep];    %'C:\Users\Varnet L�o\Dropbox\Professionnel\Matlab\MyScripts\modulationACI\AM';
cfg_crea.folder_name     = 'NoiseStims'; % nom du dossier � creer contenant les noises

dir_main = cfg_crea.path;
dir_where = [dir_main cfg_crea.folder_name filesep];

fs  = cfg_crea.fs;
dur = cfg_crea.stim_dur;
noise_type = cfg_crea.noise_type;

N_noise  = cfg_crea.N_noise;
N_signal = cfg_crea.N_signal;
N_total  = N_noise*N_signal; 

disp('The new stimuli will be generated using the following parameters:')
fprintf('\tNoise type=%s\n',noise_type);
fprintf('\tDuration=%.2f s\n',dur);
fprintf('\tSampling frequency=%.2f Hz\n',fs);
fprintf('\tA total of %.0f stimuli will be generated: %.0f ''noises'', %.0f ''conditions''\n',N_total,N_noise,N_signal);
fprintf('\tTarget folder: %s\n',dir_where);
disp('Press any button to continue, press ctrl+C to abort')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
pause()

if ~isdir(dir_where)
    mkdir(dir_where);
elseif ~isempty(dir(dir_where))
    error('%s: Folder %s is not empty. Choose a new ''folder_name'' or remove manually the existing content in the target folder',upper(mfilename),dir_where);
end
 
%% Noise samples generation
N_samples = fs*dur;
for i=1:N_total
    noise = []; % refreshing the noise
    clc
    fprintf('Creating noise stimulus # %.0f of %.0f\n',i,N_total);
    
    switch noise_type
        case 'white'
            noise=randn(N_samples,1);
        case 'pink'
            error('Not validated yet...')
            noise=pinknoise(N_samples)';
        otherwise
            error('%s: Unknown type of noise. Possible options are ''pink'' or ''white''',upper(mfilename))
    end
     
    % Ensuring that the noise iteration number has four characters:
    stimnumber=num2str(i);
    if i<1000
        stimnumber=['0' stimnumber];
    end
    if i<100
        stimnumber=['0' stimnumber];
    end
    if i<10
        stimnumber=['0' stimnumber];
    end
    % Ensuring 'End'
    
    audiowrite([dir_where 'Noise_' stimnumber '.wav'], 0.99*noise/max(abs(noise)), fs);
end

%%% Save parameters
 
clock_now=fix(clock);
cfg_crea.date = clock_now;
savename = ['cfgcrea_' num2str(clock_now(1)) '_' num2str(clock_now(2)) '_' num2str(clock_now(3)) '_' num2str(clock_now(4)) '_' num2str(clock_now(5))];
save([dir_main savename], 'cfg_crea');
fprintf(['cfg file saved: ' savename '.mat\n\n']);
