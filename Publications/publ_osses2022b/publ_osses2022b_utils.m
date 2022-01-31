function publ_osses2022b_utils(Subject_ID, type_action)
% function publ_osses2022b_utils(Subject_ID, type_action)
%
% Author: Alejandro Osses
%
% Original name: g20210924_all_sessions_latin_square
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Conditions = {'white','bumpv1p2_10dB','sMPSv1p3'};
experiment = 'speechACI_Logatome-abda-S43M';

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
if nargin < 2
    type_actions = { 'Check_sounds_for_a_participant', ...
                    'Copy_results_to_Dropbox', ...
                    'Copy_results_from_Dropbox'};
    Show_cell(type_actions);
    bInput = input('Choose your action: ');

    type_action = type_actions{bInput};
end

switch type_action
    case 'Check_sounds_for_a_participant'
        % No sounds should be generated, because we have put all the 
        %    waveforms in the test computers upfront.
        bOnly_init = 1;
        f20220119_all_sessions_latin_square(Subject_ID,bOnly_init);
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
        
        if N_copied ~= 0
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
        end
        disp('')
        
    case 'Copy_results_from_Dropbox'
        N_copied = 0;
         
        dir_subj_dropbox = [dir_dropbox experiment filesep Subject_ID filesep];
        % if ~exist(dir_subj_dropbox,'dir')
        %     mkdir(dir_subj_dropbox);
        % end
        dir_results = [fastACI_dir_data experiment filesep Subject_ID filesep 'Results' filesep];
        % 
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
         
        if N_copied ~= 0
            % Then we need to copy the Conditions_file too...
            load_name = ['Conditions-for-' Subject_ID '+' experiment '.mat'];
            Cond_name2store = [dir_subj_dropbox load_name]; % [fastACI_paths('dir_output') load_name];
            if exist(Cond_name2store,'file')
                name_dst = [fastACI_paths('dir_output') load_name];
                if exist(name_dst,'file')
                    info1 = dir(Cond_name2store);
                    info2 = dir(name_dst);
         
                    fprintf('Trying to overwrite a new file dated: %s\n',info1.date);
                    fprintf('  by another dated : %s\n',info2.date);
                    if strcmp(info1.date,info1.date)
                        fprintf('The file %s is already up to date...\n',load_name);
                    else
                        bOverwrite = input('Press 1 to proceed with the overwritting or 0 to not to: ');
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
        end
        % disp('')        
end
