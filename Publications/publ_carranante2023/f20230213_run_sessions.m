function f20230213_run_sessions(subjectname,Contrast,bCheck_dropbox)
% function f20230213_run_sessions(modelname,Contrast,bCheck_dropbox)
%
% Author: Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc

if nargin < 3
    bCheck_dropbox = 1;
end
global global_vars

if nargin == 0
    subjectname = input('Enter the ID of the subject to be tested (e.g., ''S09''): '); % 'SAO'; %'SAO'; % modelname = 'king2019';
end

Contrasts2test = {'abda','adga','abpa','adta','apta'};
if nargin < 2
    Show_cell(Contrasts2test);
    bInput = input('Choose the contrast you want to test: ');
    Contrast = Contrasts2test{bInput};    
elseif isempty(Contrast)
    Show_cell(Contrasts2test);
    bInput = input('Choose the contrast you want to test: ');
    Contrast = Contrasts2test{bInput};
else 
    % Nothing to do
end

noise_type = 'bumpv1p2_10dB';
experiment = ['speechACI_Logatome-' Contrast '-S43M']; 

N_conditions = 1; % only one noise
%%% Global variables, change to use only 4000 trials:
N = 4000;
N_target = 2;
N_presentation = N/N_target;
global_vars.N_presentation = N_presentation;
%%%

if bCheck_dropbox
    % Checks whether the participant has completed any condition in another
    % computer:
    try
        il_publ_carranante2023a_utils(subjectname,[],'Copy_results_from_Dropbox',experiment);
    catch
        warning('Check what happened in publ_carranante2023a_utils.m...');
    end
end

available_hardware = {'Petite-Cabine', 'Grande-Cabine'};
% if ~ismac
%     available_hardware{end+1} = 'Sony-MDR-Alejandro';
%     available_hardware{end+1} = 'Sony-WH-Alejandro';
% end
available_hardware{end+1} = 'Default';

Show_cell(available_hardware);                    
bInput = input('Choose your hw config (if you don''t know, use ''Default''): ');
hardware_cfg = available_hardware{bInput};

lvl_target   = 79; % -18 dBFS peak level => -21 dBFS rms 
switch hardware_cfg
    % case 'Sony-MDR-Alejandro'
    %     % My own old headphones
    %     lvl_from_SLM = 90.6; % for the left headphone
    % 
    % case 'Sony-WH-Alejandro' 
    %     % NR headphones
    %     lvl_from_SLM = 91.7; % for the left headphone
        
    case 'Petite-Cabine'
        lvl_from_SLM = 82.3; % for the left headphone
        
    case 'Grande-Cabine'
        lvl_from_SLM = 82.1; % for  the left headphone, preamp zc 0032, Id No 23156 
        
    case 'Default'
        lvl_from_SLM = lvl_target; % assumes dBFS = 100
end
dBFS = 100+(lvl_target-lvl_from_SLM);

Cond_name2store = [fastACI_paths('dir_output') 'Conditions-for-' subjectname '+' experiment '.mat'];

N_per_session = 400;
cfg_game.N = 4000; % a priori knowledge
%%%
N_sessions = ceil(cfg_game.N / N_per_session); % a priori knowledge
Conditions_nr = repmat(1:N_conditions,1,N_sessions);
    
bInit_order_sessions = ~exist(Cond_name2store,'file');
if bInit_order_sessions == 0
    i_current_all = [];
    load(Cond_name2store,'i_current_all');
    
    if max(i_current_all) == 1
        bInit_order_sessions = 1;
    end
end

bIs_first_session = bInit_order_sessions;
Language = 'FR';

global_vars.dBFS = dBFS;
if bIs_first_session
    % Only important for the first session. This choice will be thereafter
    %   stored in cfg_game.
    global_vars.Language = Language;
end

% bInit_order_sessions = 1; warning('temporal')
if bInit_order_sessions == 1
    %%% Ensuring complete data collection up to session 1 and 7 (and then 10):
    N_sessions_to_complete = [1 7];
    N_total_sessions_to_complete = N_conditions * N_sessions_to_complete;

    idx = [];
    idx_i = 1;
    Sessions_to_fill = N_sessions_to_complete(1);
    for j = 1:length(N_sessions_to_complete)
        % idx_f = N_total_sessions_to_complete(j);
        Cond_here = Sessions_to_fill * N_conditions;
        idx = [idx randperm(Cond_here)+idx_i-1];
        
        idx_i = length(idx)+1; % initial index for the following permutation
        if j ~= length(N_sessions_to_complete)
            Sessions_to_fill = N_sessions_to_complete(j+1)-N_sessions_to_complete(j);
        else
            Sessions_to_fill = N_sessions - N_sessions_to_complete(j);
        end
    end
    Cond_here = Sessions_to_fill * N_conditions;
    idx = [idx randperm(Cond_here)+idx_i-1];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    i_current_all = ones([1 1]);
    idx_count = 1; % starting from the very first session...
    
    save(Cond_name2store,'idx','i_current_all','idx_count');
else
    fprintf('idx for Conditions found on file...')
    idx = []; % to be loaded in the next line
    load(Cond_name2store,'idx');
    
    idx_count = [];
    load(Cond_name2store,'idx_count');
    
    hardware_cfg_per_session = [];
    load(Cond_name2store,'hardware_cfg_per_session');
end

Conditions_nr = Conditions_nr(idx);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = idx_count:length(Conditions_nr)
    
    hardware_cfg_per_session{i} = hardware_cfg;
    if i == idx_count
        fprintf('\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
        % This is the first session of today:
        fprintf('Starting the first session of today\n');
        fprintf('\tSubject_ID=%s\n',subjectname);
        fprintf('\tHeadphones=%s\n',hardware_cfg);
        fprintf('\n\tInformed consent signed already?\n');
        fprintf('\tConditions being loaded from: %s\n',Cond_name2store);
        fprintf('\tThis is session %.0f of %.0f, is this correct?\n',i,length(Conditions_nr));
        disp('%-------------------------------------------------------------------------%')
        fprintf('If all the information above is correct, the session is ready to start.\n');
        fprintf('Once the participant is seating inside, press any button to start the session\n\t(press ctrl+c to abort)\n');
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        pause();
    end
    
    idx_condition = Conditions_nr(i);
    
    load(Cond_name2store,'i_current_all');
    data_passation.i_current = i_current_all(idx_condition); 
    
    if data_passation.i_current < cfg_game.N
        [cfg_game,data_passation] = fastACI_experiment(experiment,subjectname,noise_type);
        
        i_current_all(idx_condition) = data_passation.i_current;
        
        if mod(i_current_all(idx_condition),N_per_session) == 0
            idx_count = idx_count + 1;
            bCompleted_this_block = 1;
        else
            % This means that the session was paused: the counter is not incremented
            bCompleted_this_block = 0;
        end
        
        save(Cond_name2store,'idx','i_current_all','idx_count','hardware_cfg_per_session'); % updates
        
        if bCheck_dropbox
            try
                il_publ_carranante2023a_utils(subjectname,[],'Copy_results_to_Dropbox',experiment);
            catch
                warning('Check what happened in publ_carranante2023a_utils.m...');
            end
        end
        
        if bCompleted_this_block == 0
            % Stopping this script if the participant chose to take a break
            return;
        end
    end
    
    disp('')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function outs = il_publ_carranante2023a_utils(Subject_ID, noise_type, type_action, experiment)
% function outs = publ_osses2022b_utils(Subject_ID, noise_type, type_action, experiment)
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

Conditions = {noise_type}; % converting to a cell array if noise_type is a char
% experiment = 'speechACI_Logatome-abda-S43M';

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
    type_actions = [];
    if bDropbox_found
        type_actions{end+1} = 'Copy_results_to_Dropbox';
        type_actions{end+1} = 'Copy_results_from_Dropbox';
    end
    Show_cell(type_actions);
    bInput = input('Choose your action: ');

    type_action = type_actions{bInput};
end

switch type_action
    case 'Run_experiment_for_a_participant'
        f20230213_run_sessions(Subject_ID,bDropbox_found);
        
    % case 'Get_participant_ACI'
    %     flags = {'no_plot'};
    %     publ_osses2022b_preregistration_2_participants_ACIs('subjectname',Subject_ID,flags{:});
        
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
                fprintf('\tTrying to overwrite a new file dated: %s\n',info2.date);
                fprintf('\tby another dated : %s\n',info1.date);
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
