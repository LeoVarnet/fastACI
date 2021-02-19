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
dir_main = [path filesep];    %'C:\Users\Varnet L�o\Dropbox\Professionnel\Matlab\MyScripts\modulationACI\AM';
dir_results = [dir_main 'Interim_results' filesep];

%%% Checking if an *.init file is found on disc:
script_name = sprintf('%s_init',experiment);
if exist([script_name '.m'],'file')
    fprintf('\tScript %.m found on disc...\n',script_name);
    exp2eval = sprintf('cfg_crea=%s(cfg_crea);',script_name);
    eval(exp2eval);
else
    fprintf('\tNo initialisation file found for experiment %s...\n',experiment);
end

%%% Save parameters
[clock_str, clock_now] = Get_date_and_time_str;
cfg_crea.date = clock_now;

savename = ['cfgcrea_' clock_str '_' Subject_ID '_' cfg_crea.experiment];
save([dir_results savename], 'cfg_crea');
fprintf(['cfg file saved: ' savename '.mat\n\n']);
