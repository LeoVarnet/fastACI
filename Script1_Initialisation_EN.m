function Script1_Initialisation_EN(experiment)
% function Script1_Initialisation_EN(experiment)
%
% Description:
%       It creates a new participant file.
%       It first runs:
%           experiment_set.m
%
% Changes by AO:
%   - cfg_crea.color_noise changed by cfg_crea.noise_type
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc, close all
if nargin == 0
    experiment = 'speechACI_varnet2015';
    % experiment = 'modulationACI';
end

Subject_ID = input('Enter a Subject ID (e.g., ''S01''): ');

cfg_crea = [];
exp2eval = sprintf('cfg_crea=%s_set;',experiment);
eval(exp2eval);

% case 'speechACI_varnet2015'
switch experiment
    case 'modulationACI'
        dir_stim = cfg_crea.dir_stim;
        cfg_crea.path        = dir_stim;    %'C:\Users\Varnet Leo\Dropbox\Professionnel\Matlab\MyScripts\modulationACI\AM';
end
cfg_crea.experiment  = experiment;
cfg_crea.Subject_ID  = Subject_ID;

[path,name,ext]=fileparts(which(mfilename)); % path will be the folder where this file is located...
dir_main = [path filesep];    %'C:\Users\Varnet L�o\Dropbox\Professionnel\Matlab\MyScripts\modulationACI\AM';
dir_results = [dir_main 'Interim_results' filesep];
switch experiment
    case 'modulationACI'        
        dir_where  = [cfg_crea.path cfg_crea.folder_name filesep]; 
        
        fs         = cfg_crea.fs;
        dur        = cfg_crea.stim_dur;
        noise_type = cfg_crea.noise_type;
        N_noise    = cfg_crea.N_noise;
        N_signal   = cfg_crea.N_signal;
        N_total    = N_noise*N_signal; 

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
            bGenerate_stimuli = 1;
        elseif ~isempty(dir(dir_where))
            warning('%s: Folder %s is not empty.',upper(mfilename),dir_where);
            bInit_participant_only = input('Do you to initialise a participant without generating new sounds? (1=yes; 0=no): ');
            bGenerate_stimuli = ~bInit_participant_only;

            if bGenerate_stimuli
                error('%s: Folder %s is not empty. Choose a new ''folder_name'' or remove manually the existing content in the target folder',upper(mfilename),dir_where);
            end
        end
        if bGenerate_stimuli
            %% Noise samples generation
            for i=1:N_total
                clc
                fprintf('Creating noise stimulus # %.0f of %.0f\n',i,N_total);

                str_stim = [];
                str_inout = []; 
                str_inout.istarget = 0;
                str_inout.expvar   = 0; % idle
                str2eval = sprintf('str_stim=%s_user(str_inout,cfg_crea);',cfg_crea.experiment);
                eval(str2eval);
                noise_cal = str_stim.stim_noise_alone;

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

                audiowrite([dir_where 'Noise_' stimnumber '.wav'], noise_cal, fs);
            end
        end
end

cfg_crea.bDebug = 1;

%%% Save parameters
[clock_str, clock_now] = Get_date_and_time_str;
cfg_crea.date = clock_now;

savename = ['cfgcrea_' clock_str '_' Subject_ID '_' cfg_crea.experiment];
save([dir_results savename], 'cfg_crea');
fprintf(['cfg file saved: ' savename '.mat\n\n']);
