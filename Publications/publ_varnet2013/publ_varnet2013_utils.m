function publ_varnet2013_utils% (Subject_ID, type_action)
% function publ_varnet2013_utils% (Subject_ID, type_action)
%
% Author: Alejandro Osses
%
% Original name: g20210924_all_sessions_latin_square
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dir_data = fastACI_dir_data;
% Conditions = {'white','bumpv1p2_10dB','sMPSv1p3'};
experiment = 'speechACI_varnet2013';

% dir_dropbox = '/home/alejandro/Dropbox/ENS_lab_shared/';
 
% if nargin == 0
%     Subject_ID = input('Enter the participant ID (e.g. ''S01''): ');
% end

% clc
% if nargin < 2
    type_actions = {'Check_varnet2013_S02_data'}; % , ...
%                     'Copy_results_to_Dropbox', ...
%                     'Copy_results_from_Dropbox'};
%     Show_cell(type_actions);
    bInput = 1; % input('Choose your action: ');

    type_action = type_actions{bInput};
% end
 
switch type_action
    case 'Check_varnet2013_S02_data'
        % switch Subject_ID
        %     case {'S01','varnet2013_S01'}
        %         % Nothing to do
        %     otherwise
        %         error('This function has only been tested for ''varnet2013_S01''...')
        % end
        
        dir_src = [dir_data 'data_varnet2013' filesep 'Sujet_Leo_S1' filesep];
        
        subj_dst = 'varnet2013_S02';
        
        dir_subj_dst = [dir_data experiment filesep subj_dst filesep];
        
        if ~exist(dir_subj_dst,'dir')
            bContinue = 1;
            % mkdir(dir_subj_dst);
        else
            bContinue = 0;
        end
            
        if bContinue == 1
            dirs2check_src = {'ListeBruit','ListeSignal'};

            dir_where = [dir_src dirs2check_src{1} filesep];
            files = Get_filenames(dir_where,'*.wav');

            dur_ramp = 75e-3; % 
            lvl = 65;
            dBFS = 100;
            for i = 1:length(files)
                [noise,fs] = audioread([dir_where files{i}]);
                lvl_bef(i) = rmsdb(noise)+dBFS;
            end
            gain2apply = From_dB(lvl-median(lvl_bef));

            dur_ramp_samples = round(dur_ramp*fs);

            if ~exist(dir_subj_dst,'dir')
                mkdir(dir_subj_dst);
            end
            dir_dst = [dir_subj_dst dirs2check_src{1} filesep];
            if ~exist(dir_dst,'dir')
                mkdir(dir_dst);
            end

            for i = 1:length(files)

                [noise,fs] = audioread([dir_where files{i}]);

                noi_bef = Randomise_insig(noise,dur_ramp_samples);
                noi_aft = Randomise_insig(noise,dur_ramp_samples);

                noise = [noi_bef; noise; noi_aft];
                noise   = gain2apply*noise;

                if i == 1
                    % Generates the cosine ramp:
                    rp = ones(size(noise)); 
                    rp(1:dur_ramp_samples)           = rampup(dur_ramp_samples);
                    rp(end-dur_ramp_samples+1:end) = rampdown(dur_ramp_samples);
                end
                noise = rp.*noise;

                audiowrite([dir_dst files{i}],noise,fs); 
            end

            dir_where = [dir_src dirs2check_src{2} filesep];
            files = Get_filenames(dir_where,'*.wav');
            dir_dst = [dir_subj_dst dirs2check_src{2} filesep];
            if ~exist(dir_dst,'dir')
                mkdir(dir_dst);
            end

            sil = zeros(dur_ramp_samples,1);
            for i = 1:length(files)
                [insig,fs] = audioread([dir_where files{i}]);
                insig = gain2apply*[sil; insig; sil];

                audiowrite([dir_dst files{i}],insig,fs);
            end
        else
            fprintf('The folder structure seems to be the correct one...\n');
        end
        
%     case 'Copy_results_to_Dropbox'
%         N_copied = 0;
%          
%         dir_subj_dropbox = [dir_dropbox experiment filesep Subject_ID filesep];
%         if ~exist(dir_subj_dropbox,'dir')
%             mkdir(dir_subj_dropbox);
%         end
%         dir_results = [fastACI_dir_data experiment filesep Subject_ID filesep 'Results' filesep];
%         
%         for i = 1:length(Conditions)
%             exp2filter = ['save*' Conditions{i} '.mat'];
%             [load_name_full, load_name] = Get_savenames(dir_results, exp2filter);
%             
%             name_dst = [dir_subj_dropbox load_name];
%             if ~exist(name_dst,'file')
%                 copyfile(load_name_full,name_dst);
%                 fprintf('%s successfully copied to %s\n',load_name,dir_subj_dropbox);
%                 N_copied = N_copied + 1;
%             else
%                 fprintf('%s is already located in %s\n',load_name,dir_subj_dropbox);
%                 fprintf('(file not being copied)\n');
%             end
%         end
%         
%         if N_copied ~= 0
%             % Then we need to copy the Conditions_file too...
%             load_name = ['Conditions-for-' Subject_ID '+' experiment '.mat'];
%             Cond_name2store = [fastACI_paths('dir_output') load_name];
%             if exist(Cond_name2store,'file')
%                 name_dst = [dir_subj_dropbox load_name];
%                 if exist(name_dst,'file')
%                     info1 = dir(Cond_name2store);
%                     info2 = dir(name_dst);
% 
%                     fprintf('%s:\n',load_name);
%                     fprintf('\tTrying to overwrite a new file dated: %s\n',info1.date);
%                     fprintf('\tby another dated : %s\n',info2.date);
%                     bOverwrite = input('Press 1 to proceed with the overwritting or 0 to not to: ');
%                     if bOverwrite
%                         copyfile(Cond_name2store,name_dst);
%                         fprintf('%s was successfully copied to %s\n',load_name,dir_subj_dropbox);
%                     else
%                         fprintf('The file %s was NOT overwritten\n',load_name);
%                     end
%                 else
%                     % No file in destination, so no problem to directly
%                     %   copy the Cond_name2store file:
%                     copyfile(Cond_name2store,name_dst);
%                     fprintf('%s was successfully copied to %s\n',load_name,dir_subj_dropbox);
%                 end
%             end
%         end
%         disp('')
%         
%     case 'Copy_results_from_Dropbox'
%         N_copied = 0;
%          
%         dir_subj_dropbox = [dir_dropbox experiment filesep Subject_ID filesep];
%         % if ~exist(dir_subj_dropbox,'dir')
%         %     mkdir(dir_subj_dropbox);
%         % end
%         dir_results = [fastACI_dir_data experiment filesep Subject_ID filesep 'Results' filesep];
%         % 
%         for i = 1:length(Conditions)
%             exp2filter = ['save*' Conditions{i} '.mat'];
%             [load_name_full, load_name] = Get_savenames(dir_subj_dropbox, exp2filter);
%         
%             name_dst = [dir_results load_name];
%             if ~exist(name_dst,'file')
%                 copyfile(load_name_full,name_dst);
%                 fprintf('%s successfully copied to %s\n',load_name,dir_subj_dropbox);
%                 N_copied = N_copied + 1;
%             else
%                 fprintf('%s is already located in %s\n',load_name,dir_subj_dropbox);
%                 fprintf('(file not being copied)\n');
%             end
%         end
%          
%         if N_copied ~= 0
%             % Then we need to copy the Conditions_file too...
%             load_name = ['Conditions-for-' Subject_ID '+' experiment '.mat'];
%             Cond_name2store = [dir_subj_dropbox load_name]; % [fastACI_paths('dir_output') load_name];
%             if exist(Cond_name2store,'file')
%                 name_dst = [fastACI_paths('dir_output') load_name];
%                 if exist(name_dst,'file')
%                     info1 = dir(Cond_name2store);
%                     info2 = dir(name_dst);
%          
%                     fprintf('Trying to overwrite a new file dated: %s\n',info1.date);
%                     fprintf('  by another dated : %s\n',info2.date);
%                     if strcmp(info1.date,info1.date)
%                         fprintf('The file %s is already up to date...\n',load_name);
%                     else
%                         bOverwrite = input('Press 1 to proceed with the overwritting or 0 to not to: ');
%                         if bOverwrite
%                             copyfile(Cond_name2store,name_dst);
%                             fprintf('%s was successfully copied to %s\n',load_name,dir_subj_dropbox);
%                         else
%                             fprintf('The file %s was NOT overwritten\n',load_name);
%                         end
%                     end
%                 else
%                     % No file in destination, so no problem to directly
%                     %   copy the Cond_name2store file:
%                     copyfile(Cond_name2store,name_dst);
%                     fprintf('%s was successfully copied to %s\n',load_name,dir_subj_dropbox);
%                 end
%             end
%         end
%         % disp('')        
end
