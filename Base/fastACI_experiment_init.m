function cfg_crea = fastACI_experiment_init(experiment_full,Subject_ID, Condition)
% function cfg_crea = fastACI_experiment_init(experiment_full,Subject_ID, Condition)
%
% Description:
%       It creates a new participant file.
%       It first runs:
%           'experiment'_set.m
%
% Changes by AO:
%   - cfg_crea.color_noise changed by cfg_crea.noise_type
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc, close all
if nargin < 3
    Condition = [];
end
if nargin < 2
    Subject_ID = input('Enter the Subject ID (e.g., ''S01''): ');
end

if ~isempty(Condition)
    str_cond = sprintf(', condition=%s',Condition);
else
    str_cond = '';
end
fprintf('Running %s.m (experiment=%s, subject=%s%s)\n',mfilename,experiment_full,Subject_ID,str_cond);

% [dir_results, dir_results_completed] = Check_local_dir_data(experiment,Subject_ID);
dir_results = Check_local_dir_data(experiment_full,Subject_ID);

%%% Checking whether there are separable conditions (separator='-')
if sum(experiment_full=='-') % Checking whether it contains an hyphen
    experiment   = strsplit(experiment_full,'-');
    N_Cond_extra = length(experiment)-1;
    Cond_extra   = experiment(2:end);
    experiment = experiment{1};
else
    N_Cond_extra = 0;
    Cond_extra = [];
    experiment = experiment_full;
end

cfg_crea = [];
cfg_crea.experiment  = experiment;
cfg_crea.experiment_full = experiment_full;
for i = 1:N_Cond_extra
    exp2eval = sprintf('cfg_crea.Cond_extra_%.0f = Cond_extra{%.0f};',i,i);
    eval(exp2eval);
end
cfg_crea.Subject_ID  = Subject_ID;
if ~isempty(Condition)
    cfg_crea.Condition = Condition;
end

%%%
% Checking existing crea _files:
files_prev = Get_filenames(dir_results,['cfgcrea_*' Condition '.mat']);
if ~isempty(files_prev)
    fprintf('%s: %.0f create file(s) (cfgcrea) was/were found on disk. No creation file will be created.\n',upper(mfilename),length(files_prev));
    fprintf('\tIf you really want to create a new create file, remove the existing file: %s\n',files_prev{1});
    pause_time = 3;
    fprintf('\t Pausing for %.0f seconds... (press ctrl+c to cancel)\n',pause_time);
    pause(pause_time);
    
    if length(files_prev) > 1
        % So, one creafile was found only
        fprintf('\t Skipping the waveform generation because there is more than 1 crea file on disk)\n');
    elseif length(files_prev) == 1
        % We will try to generate the waveforms again. They will be generated
        %    if the waveform folders do not exist...
        cfg_crea = []; % loaded again in the next line
        cfg_crea = load([dir_results files_prev{1}],'cfg_crea');
        %%% Checking if an *.init file is found on disc:
        script_name = sprintf('%s_init',experiment);
        if exist([script_name '.m'],'file')
            fprintf('\tScript %s.m found on disc...\n',script_name);
            exp2eval = sprintf('cfg_crea=%s(cfg_crea);',script_name);
            eval(exp2eval);
        end
    end
    fprintf('\t Pausing for %.0f seconds... (press ctrl+c to cancel)\n',pause_time);
    pause(pause_time);
    return
end
%%%

exp2eval = sprintf('cfg_crea=%s_set(cfg_crea);',experiment);
eval(exp2eval);

exp2eval = sprintf('cfg_crea=%s_cfg(cfg_crea);',experiment);
eval(exp2eval);

if isfield(cfg_crea,'N_signal')
    error('Error... (temporal, remove soon - AO on 17/03/2021)')
end
if isfield(cfg_crea,'N_noise')
    error('Error... (temporal, remove soon - AO on 17/03/2021)')
end

if ~isfield(cfg_crea,'N')
    cfg_crea.N = cfg_crea.N_target*cfg_crea.N_presentation;
end

%%% Checking if an *.init file is found on disc:
script_name = sprintf('%s_init',experiment);
if exist([script_name '.m'],'file')
    fprintf('\tScript %s.m found on disc...\n',script_name);
    exp2eval = sprintf('cfg_crea=%s(cfg_crea);',script_name);
    eval(exp2eval);
else
    fprintf('\tNo initialisation file found for experiment %s...\n',experiment);
end

if ~isfield(cfg_crea,'response_correct_target')
    if cfg_crea.N_target > 2
        if ~isfield(cfg_crea,'response_correct_target')
            error('A maximum of two alternatives ''1'' and ''2'' can be asked to the participants. Your experiment seems to use %.0f sounds, but they should be mapped to only two responses. Specify the field ''response_correct_target''.',cfg_game.N_target);
        end
    elseif cfg_crea.N_target == 2
        cfg_crea.response_correct_target = [1 2];
    end
end

%%% Save parameters
[clock_str, clock_now] = Get_date_and_time_str;
cfg_crea.date = clock_now;

filter2use = [Subject_ID '_' experiment_full];
if ~isempty(Condition)
    filter2use = [filter2use '_' Condition];
end
savename = ['cfgcrea_' clock_str '_' filter2use];

% dir_results_subj = [dir_results Subject_ID filesep];
mkdir(dir_results); % makes sure the folder for the participant exists...
info_toolbox = Get_toolbox_info(mfilename);
save([dir_results savename], 'cfg_crea','info_toolbox');
fprintf(['cfg file saved: ' savename '.mat\n\n']);
