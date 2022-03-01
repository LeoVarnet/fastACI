function f20220119_all_sessions_latin_square(subjectname,bOnly_init)
% function f20220119_all_sessions_latin_square(modelname,bOnly_init)
%
% Author: Alejandro Osses
%
% Original name: g20210924_all_sessions_latin_square
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc

global global_vars

if nargin < 2
    bOnly_init = 0;
end
if nargin == 0
    subjectname = input('Enter the ID of the subject to be tested (e.g., ''S09''): '); % 'SAO'; %'SAO'; % modelname = 'king2019';
end
Conditions = {'white','bumpv1p2_10dB','sMPSv1p3'};
experiment = 'speechACI_Logatome-abda-S43M'; 

N_conditions = length(Conditions);
%%% Global variables, change to use only 4000 trials:
N = 4000;
N_target = 2;
N_presentation = N/N_target;
global_vars.N_presentation = N_presentation;
%%%

if bOnly_init
    for i = 1:N_conditions
        noise_type = Conditions{i};
        fastACI_experiment_init(experiment,subjectname, noise_type);
    end
    return;
else
    % Checks whether the participant has completed any condition in another
    % computer:
    try
        publ_osses2022b_utils(subjectname,[],'Copy_results_from_Dropbox');
    catch
        warning('Check what happened in publ_osses2022b_utils.m...');
    end
end

available_hardware = {'Petite-Cabine', 'Grande-Cabine'};
if ~ismac
    available_hardware{end+1} = 'Sony-MDR-Alejandro';
    available_hardware{end+1} = 'Sony-WH-Alejandro';
end
Show_cell(available_hardware);                    
bInput = input('Choose your hw config: ');
hardware_cfg = available_hardware{bInput};

lvl_target   = 79; % -18 dBFS peak level => -21 dBFS rms 
switch hardware_cfg
    case 'Sony-MDR-Alejandro'
        % My own old headphones
        lvl_from_SLM = 90.6; % for the left headphone
        
    case 'Sony-WH-Alejandro' 
        % NR headphones
        lvl_from_SLM = 91.7; % for the left headphone
        
    case 'Petite-Cabine'
        lvl_from_SLM = 82.3; % for the left headphone
        
    case 'Grande-Cabine'
        lvl_from_SLM = 82.1; % for  the left headphone, preamp zc 0032, Id No 23156 
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
if bIs_first_session
    Languages = {'FR','EN'};
    Show_cell(Languages);
    bInput = input('This is the first session for this participant, please choose the language: ');
    Language = Languages{bInput};
end
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
    i_current_all = ones(1,length(Conditions));
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
    noise_type = Conditions{idx_condition};
    
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
        
        try
            publ_osses2022b_utils(subjectname,[],'Copy_results_to_Dropbox');
        catch
            warning('Check what happened in publ_osses2022b_utils.m...');
        end
        
        if bCompleted_this_block == 0
            % Stopping this script if the participant chose to take a break
            return;
        end
    end
    
    disp('')
end
