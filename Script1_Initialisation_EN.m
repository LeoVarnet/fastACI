function cfg_crea = Script1_Initialisation_EN(experiment,Subject_ID)
% function cfg_crea = Script1_Initialisation_EN(experiment,Subject_ID)
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
if nargin < 2
    Subject_ID = input('Enter the Subject ID (e.g., ''S01''): ');
end
if nargin == 0
    experiments = {'modulationACI','modulationACI_seeds','speechACI_varnet2015'};
    Show_cell(experiments);
    bExp = input('Choose the experiment that you want to run from the list above: ');
    
    experiment = experiments{bExp};
end

fprintf('Running %s.m (experiment=%s, subject=%s)\n',mfilename,experiment,Subject_ID);

cfg_crea = [];
exp2eval = sprintf('cfg_crea=%s_set;',experiment);
eval(exp2eval);

cfg_crea.experiment  = experiment;
cfg_crea.Subject_ID  = Subject_ID;
if ~isfield(cfg_crea,'N')
    cfg_crea.N = cfg_crea.N_signal*cfg_crea.N_noise;
end

[path,name,ext]=fileparts(which(mfilename)); % path will be the folder where this file is located...
dir_main = [path filesep];    %'C:\Users\Varnet Leo\Dropbox\Professionnel\Matlab\MyScripts\modulationACI\AM';
dir_results = [dir_main 'Interim_results' filesep];
if ~isfolder(dir_results)
    if isfolder(dir_main)
        % If folder Interim_results does not exist, it is created, but this 
        %   is only done if dir_main already exists...
        mkdir(dir_results);
    end
end
    
%%% Checking if an *.init file is found on disc:
script_name = sprintf('%s_init',experiment);
if exist([script_name '.m'],'file')
    fprintf('\tScript %.m found on disc...\n',script_name);
    exp2eval = sprintf('cfg_crea=%s(cfg_crea);',script_name);
    eval(exp2eval);
else
    fprintf('\tNo initialisation file found for experiment %s...\n',experiment);
end

if ~isfield(cfg_crea,'response_correct_target')
    if cfg_crea.N_signal > 2
        if ~isfield(cfg_crea,'response_correct_target')
            error('A maximum of two alternatives ''1'' and ''2'' can be asked to the participants. Your experiment seems to use %.0f sounds, but they should be mapped to only two responses. Specify the field ''response_correct_target''.',cfg_game.N_signal);
        end
    elseif cfg_crea.N_signal == 2
        cfg_crea.response_correct_target = [1 2];
    end
end

%%% Save parameters
[clock_str, clock_now] = Get_date_and_time_str;
cfg_crea.date = clock_now;

savename = ['cfgcrea_' clock_str '_' Subject_ID '_' cfg_crea.experiment];
save([dir_results savename], 'cfg_crea');
fprintf(['cfg file saved: ' savename '.mat\n\n']);
