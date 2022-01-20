function g20220119_all_sessions_latin_square(modelname,bOnly_init)
% function g20220119_all_sessions_latin_square(modelname,bOnly_init)
%
% Author: Alejandro Osses
%
% Original name: g20210924_all_sessions_latin_square
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global global_vars

if nargin == 0
    modelname = input('Enter the ID of the subject to be tested (e.g., ''S09''): '); % 'SAO'; %'SAO'; % modelname = 'king2019';
end
Conditions = {'white','bumpv1p2_10dB','sMPSv1p3'};
experiment = 'speechACI_Logatome-abda-S43M'; 

N_conditions = length(Conditions);

if bOnly_init
    for i = 1:N_conditions
        noise_type = Conditions{i};
        fastACI_experiment_init(experiment,modelname, noise_type);
    end
    return;
end

available_hardware = {  'Sony-MDR-Alejandro','Sony-WH-Alejandro', ...
                        'Petite-Cabine', 'Grande-Cabine'};
Show_cell(available_hardware);                    
bInput = 2; warning('temporal'); % input('Choose your hw config: ');
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
        error('Not yet enabled')
        lvl_from_SLM = 82.3; % for the left headphone
        
    case 'Grande-Cabine'
        error('Not yet enabled')
        lvl_from_SLM = 90.6; % for the left headphone
end
dBFS = 100+(lvl_target-lvl_from_SLM);
global_vars.dBFS = dBFS;
global_vars.Language = 'EN';

Cond_name2store = [fastACI_paths('dir_output') 'Conditions-for-' modelname '+' experiment '.mat'];

%%% Global variables:
N = 4000;
N_target = 2;
N_presentation = N/N_target;
global_vars.N_presentation = N_presentation;

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

bInit_order_sessions = 1; warning('temporal')
if bInit_order_sessions == 1
    %%% Ensuring complete data collection up to session 1 and 7 (and then 10):
    N_sessions_to_complete = [1 7];
    N_total_sessions_to_complete = N_conditions * N_sessions_to_complete;

    idx = [];
    idx_i = 1;
    Sessions_to_fill = N_sessions_to_complete(1);
    for j = 1:length(N_sessions_to_complete)
        idx_f = N_total_sessions_to_complete(j);
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
end

Conditions_nr = Conditions_nr(idx);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% switch modelname
%     case 'osses2021'        
%         templ_num = 100;
%         
%         dir_here = [fastACI_basepath 'Simulations' filesep];
%         fname_cfg =  [dir_here modelname '_cfg.m'];
% 
%         run_str = 'run-4'; % relanoiborra2019
% 
%         p = il_get_model_config_DAGA(run_str,modelname,templ_num);
%         text_to_write = readfile_replace('model_cfg_replace.txt',p);
% 
%         if exist(fname_cfg,'file')
%             fprintf('----------------------------------------------------------------------------\n')
%             fprintf('file %s exists, \npress any key to continue (will overwrite) or press ctrl+C to cancel \n',fname_cfg);
%             fprintf('----------------------------------------------------------------------------\n')
%             pause
%         end
% 
%         fid = fopen(fname_cfg, 'w');
%         fwrite(fid, text_to_write);
%         fclose(fid);
% end

for i = idx_count:length(Conditions_nr)
    idx_condition = Conditions_nr(i);
    noise_type = Conditions{idx_condition};
    
    load(Cond_name2store,'i_current_all');
    data_passation.i_current = i_current_all(idx_condition); 
    
    if data_passation.i_current < cfg_game.N
        [cfg_game,data_passation] = fastACI_experiment(experiment,modelname,noise_type);
        
        i_current_all(idx_condition) = data_passation.i_current;
        idx_count = idx_count + 1;
        save(Cond_name2store,'idx','i_current_all','idx_count'); % updates
    end
    
    % if data_passation.i_current >= cfg_game.N
    %     folder_src = cfg_game.dir_results;
    %     folder_new = [fileparts(cfg_game.dir_results(1:end-1)) filesep 'Results-multisession-' run_str filesep];
    % 
    %     movefile(folder_src,folder_new);
    % end
    
    disp('')
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function p = il_get_model_config_DAGA(run_str,modelname,templ_num)
% 
% p = [];
% p.modelname = modelname;
% 
% p.templ_num = templ_num;
% p.bStore_template = 1;
% 
% switch run_str
%     case 'run-1'
%         p.decision_script = 'aci_detect';
%         p.template_script = 'model_template';
%         p.template_every_trial = 0;
%         p.det_lev = -6;
%         p.type_decision = 'optimal_detector';
%         p.thres_for_bias = 0;
%         
%     case 'run-3-m1p55'
%         p.decision_script = 'aci_detect';
%         p.template_script = 'model_template';
%         p.template_every_trial = 0;
%         p.det_lev = -6;
%         p.type_decision = 'optimal_detector';
%         p.thres_for_bias = -1.55;
%         
%     case 'run-3-p0p39'
%         p.decision_script = 'aci_detect';
%         p.template_script = 'model_template';
%         p.template_every_trial = 0;
%         p.det_lev = -6;
%         p.type_decision = 'optimal_detector';
%         p.thres_for_bias = 0.39;
%         
%     case 'run-3-p0p78'
%         p.decision_script = 'aci_detect';
%         p.template_script = 'model_template';
%         p.template_every_trial = 0;
%         p.det_lev = -6;
%         p.type_decision = 'optimal_detector';
%         p.thres_for_bias = 0.78;
%     case 'run-4'
%         p.decision_script = 'aci_detect';
%         p.template_script = 'model_template';
%         p.template_every_trial = 0;
%         p.det_lev = -6;
%         p.type_decision = 'relanoiborra2019_decision';
%         
%     otherwise
%         error('Run condition not recognised')
% end
% 
% if ~isfield(p,'thres_for_bias')
%     p.thres_for_bias = nan;
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [idx_out] = il_Check_completed_permutations(idx,N_sessions_to_complete, Conditions_nr);
% 

% 
% % Assign = 
% idx_out = [];
% Cond_out = [];
% idx_remaining = idx;
% current_idx = 1;
% i = 1;
% 
% last_idx = 1;
% for i = 1:length(N_total_sessions_to_complete)
%     idx_here = idx(last_idx:N_total_sessions_to_complete(i));
%     Conds_here = Conditions_nr(idx_here);
%     if length(unique(Conds_here)) == N_conditions
%         bReassign = 0;
%     else
%         bReassign = 1;
%     end
%     if bReassign
%         N_repeated = [];
%         N_missing  = [];
%         for j=1:N_conditions
%             if sum(Conds_here == j)>1
%                 N_repeated(end+1) = j;
%             elseif sum(Conds_here == j) == 0
%                 N_missing(end+1) = j;
%             end
%         end
%         % Looks for the index of the repeated one:
%         idx_here = find(
%         disp('')
%     end
%     disp('')
% end

% % for i = 1:length(N_sessions_to_complete)
%     %%% first sample:
%     idx_out(i) = idx(i);
%     Cond_out(i) = Conditions_nr(idx(i));
%     idx_remaining(i) = nan;
%     i = i+1;
%     
%     N_current = length(Cond_out);
%     % Look for the following non repeated sample:
%     while length(Cond_out) == N_current
%         if Cond_idx(i)
%     end
%     disp('')
%     %%% second
%     
% end
