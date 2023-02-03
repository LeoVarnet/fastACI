function outs = publ_osses2022b_utils(Subject_ID, noise_type, type_action)
% function outs = publ_osses2022b_utils(Subject_ID, noise_type, type_action)
%
% 1. Description:
%   Possible Subject_ID values: from 'S01' up to 'S14'
%   Possible noise_type values: 'white','bumpv1p2_10dB','sMPSv1p3', 
%                               or [] (to check all)
%   Possible type_actions:
%       - 'Check_sounds_for_a_participant': Independent of the variable 
%             noise_type, this action will check whether all background noises
%             for participant Subject_ID are stored on disk or not. If the 
%             sounds are not found they will be reconstructed from the 
%             cfg_crea MAT file.
%       - 'Run_experiment_for_a_participant'
%       - 'Get_participant_ACI'
%       - 'Play_intellitest'
%       - 'Copy_results_to_Dropbox'
%       - 'Copy_results_from_Dropbox'
%
% 2. Example:
%   To check the waveforms of participant 'S13', all noises:
%       Subject_ID = 'S13';
%       type_action = 'Check_sounds_for_a_participant';
%       noise_type = [];
%       publ_osses2022b_utils(Subject_ID, noise_type, type_action);
%
%   To run the experimental sessions of participant 'S13':
%       Subject_ID = 'S13';
%       type_action = 'Run_experiment_for_a_participant';
%       noise_type = [];
%       publ_osses2022b_utils(Subject_ID, noise_type, type_action);
%
% Author: Alejandro Osses
%
% Original name: g20210924_all_sessions_latin_square
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

outs = [];
% outs = publ_osses2022c_ARO_poster_utils(Subject_ID, noise_type, type_action)

if nargin < 2
    Conditions = {'white','bumpv1p2_10dB','sMPSv1p3'};
else
    if ischar(noise_type)
        Conditions = {noise_type}; % converting to a cell array if noise_type is a char
    end
    if isempty(noise_type)
        disp('All noise conditions will be used')
        Conditions = {'white','bumpv1p2_10dB','sMPSv1p3'};
    end
end
experiment = 'speechACI_Logatome-abda-S43M';

dir_dropbox = [];
if ismac
    dir_dropbox = '~/Dropbox/ENS_lab_shared/';
elseif isunix
    dir_dropbox = '/home/alejandro/Dropbox/ENS_lab_shared/';
elseif iswindows
    % Leo's folder should go here
end

if nargin == 0
    Subject_ID = input('Enter the participant ID (e.g. ''S01''): ');
end
clc

bDropbox_found = exist(dir_dropbox,'dir');
if nargin < 3
    type_actions = {'Play_intellitest', ...
                    'Check_sounds_for_a_participant'};
    if bDropbox_found
        type_actions{end+1} = 'Copy_results_to_Dropbox';
        type_actions{end+1} = 'Copy_results_from_Dropbox';
    end
    Show_cell(type_actions);
    bInput = input('Choose your action: ');

    type_action = type_actions{bInput};
end

switch type_action
    case 'Play_intellitest'
        dir_where = [fastACI_basepath 'Stimuli' filesep 'Intellitest' filesep];
        files = Get_filenames(dir_where,'*1a.wav');
        
        Show_cell(files);
        
        bInput = input('Input 0 to listen all sounds once, or the number of the sound to hear: ');
        
        if bInput == 0
        else
            files = files(bInput);
        end
        
        for i = 1:length(files)
            [insig,fs] = audioread([dir_where files{i}]);
            
            fprintf('Playing %s\n',files{i}(1:end-5));
            sound(insig,fs);
            pause(3);
        end
        
    case 'Check_sounds_for_a_participant'
        % No sounds should be generated, because we have put all the 
        %    waveforms in the test computers upfront.
        bOnly_init = 1;
        f20220119_all_sessions_latin_square(Subject_ID,bOnly_init,bDropbox_found);
    
    case 'Run_experiment_for_a_participant'
        bOnly_init = 0;
        f20220119_all_sessions_latin_square(Subject_ID,bOnly_init,bDropbox_found);
        
    case 'Get_participant_ACI'
        flags = {'no_plot'};
        publ_osses2022b_preregistration_2_participants_ACIs('subjectname',Subject_ID,flags{:});
        
    case 'Copy_results_to_Dropbox'
        N_copied = 0;
         
        dir_subj_dropbox = [dir_dropbox experiment filesep Subject_ID filesep];
        if ~exist(dir_subj_dropbox,'dir')
            mkdir(dir_subj_dropbox);
        end
        dir_results = [fastACI_dir_data experiment filesep Subject_ID filesep 'Results' filesep];
        
        for i = 1:length(Conditions)
            exp2filter = ['save*' Conditions{i} '.mat'];
            [load_name_full, load_name] = Get_savenames(dir_results, exp2filter);
            
            name_dst = [dir_subj_dropbox load_name];
            if ~exist(name_dst,'file')
                copyfile(load_name_full,name_dst);
                fprintf('%s successfully copied to %s\n',load_name,dir_subj_dropbox);
                N_copied = N_copied + 1;
            else
                fprintf('%s is already located in %s\n',load_name,dir_subj_dropbox);
                fprintf('(file not being copied)\n');
            end
        end
        
        % Then we need to copy the Conditions_file too...
        load_name = ['Conditions-for-' Subject_ID '+' experiment '.mat'];
        Cond_name2store = [fastACI_paths('dir_output') load_name];
        if exist(Cond_name2store,'file')
            name_dst = [dir_subj_dropbox load_name];
            if exist(name_dst,'file')
                info1 = dir(Cond_name2store);
                info2 = dir(name_dst);

                fprintf('%s:\n',load_name);
                fprintf('\tTrying to overwrite a new file dated: %s\n',info1.date);
                fprintf('\tby another dated : %s\n',info2.date);
                bOverwrite = input('Press 1 to proceed with the overwritting or 0 to not to: ');
                if bOverwrite
                    copyfile(Cond_name2store,name_dst);
                    fprintf('%s was successfully copied to %s\n',load_name,dir_subj_dropbox);
                else
                    fprintf('The file %s was NOT overwritten\n',load_name);
                end
            else
                % No file in destination, so no problem to directly
                %   copy the Cond_name2store file:
                copyfile(Cond_name2store,name_dst);
                fprintf('%s was successfully copied to %s\n',load_name,dir_subj_dropbox);
            end
        end
        
        disp('')
        
    case 'Get_filenames' % very similar to type_action in publ_osses2022c_ARO_poster_utils.m
        if length(Conditions) ~= 1
            error('Re-run this script only requesting one noise type');
        end
        bGenerate_stimuli = 0;
        
        dir_res = [fastACI_dir_data filesep experiment filesep Subject_ID filesep 'Results' filesep];
        if ~exist(dir_res,'dir')
            bGenerate_stimuli = 1;    
        else
            file = Get_filenames(dir_res,['cfgcrea*' noise_type '.mat']);
            if isempty(file)
                bGenerate_stimuli = 1;
            end
        end
        
        if bGenerate_stimuli == 0 % i.e., if dir_res exists
            file = Get_filenames(dir_res,['cfgcrea*' noise_type '.mat']);
            if length(file) ~= 1
                error('More than one creafile was found on disk');
            end
            load([dir_res file{1}],'cfg_crea');

            cfg_crea = Check_cfg_crea_dirs(cfg_crea);

            dir_target = cfg_crea.dir_target;
            dir_noise  = cfg_crea.dir_noise;

            if ~exist(dir_target,'dir')
                error('Directory %s not found on disk',dir_target);
            end
            if ~exist(dir_noise,'dir')
                bGenerate_stimuli = 1;
                warning('Directory %s not found on disk',dir_noise);
            end

            outs.dir_target = dir_target;
            outs.dir_noise  = dir_noise;
        end
        outs.bGenerate_stimuli = bGenerate_stimuli;
        % outs.fname_results = fname_results;
        
    case 'Copy_results_from_Dropbox'
        N_copied = 0;
         
        dir_subj_dropbox = [dir_dropbox experiment filesep Subject_ID filesep];
        
        dir_results = [fastACI_dir_data experiment filesep Subject_ID filesep 'Results' filesep];
        
        for i = 1:length(Conditions)
            exp2filter = ['save*' Conditions{i} '.mat'];
            [load_name_full, load_name] = Get_savenames(dir_subj_dropbox, exp2filter);
        
            name_dst = [dir_results load_name];
            if ~exist(name_dst,'file')
                copyfile(load_name_full,name_dst);
                fprintf('%s successfully copied to %s\n',load_name,dir_subj_dropbox);
                N_copied = N_copied + 1;
            else
                fprintf('%s is already located in %s\n',load_name,dir_subj_dropbox);
                fprintf('(file not being copied)\n');
            end
        end
         
        % if N_copied ~= 0
            % Then we need to copy the Conditions_file too...
            load_name = ['Conditions-for-' Subject_ID '+' experiment '.mat'];
            Cond_name2store = [dir_subj_dropbox load_name]; % [fastACI_paths('dir_output') load_name];
            if exist(Cond_name2store,'file')
                name_dst = [fastACI_paths('dir_output') load_name];
                if exist(name_dst,'file')
                    info1 = dir(Cond_name2store);
                    info2 = dir(name_dst);
         
                    fprintf('Trying to overwrite a new file dated: %s\n',info2.date);
                    fprintf('  by another dated : %s\n',info1.date);
                    if strcmp(info1.date,info2.date)
                        fprintf('The file %s is already up to date...\n',load_name);
                    else
                        try
                           bOverwrite = datenum(info1.date)>datenum(info2.date);
                        catch me
                            % In my computer, datenum has troubles to detect 'mrt' as March ('mar')
                            bOverwrite = input('Press 1 to proceed with the overwritting or 0 to not to: ');
                        end
                        if bOverwrite
                            bOverwrite = 1;
                        else
                            bOverwrite = input('Press 1 to proceed with the overwritting or 0 to not to: ');
                        end
                        if bOverwrite
                            copyfile(Cond_name2store,name_dst);
                            fprintf('%s was successfully copied to %s\n',load_name,dir_subj_dropbox);
                        else
                            fprintf('The file %s was NOT overwritten\n',load_name);
                        end
                    end
                else
                    % No file in destination, so no problem to directly
                    %   copy the Cond_name2store file:
                    copyfile(Cond_name2store,name_dst);
                    fprintf('%s was successfully copied to %s\n',load_name,dir_subj_dropbox);
                end
            end
        % end
        % disp('')        
end
