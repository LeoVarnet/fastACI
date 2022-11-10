function [h,hname,outs] = pres_osses2022_02_AABBA_utils(experiment, subject, masker, varargin)
% function [h,hname,outs] = pres_osses2022_02_AABBA_utils(experiment, subject, masker, varargin)
%
% My comments: 
%       - 'lasso' with 'lassoslow' should be the same
%       
% Author: Leo Varnet, adapted by Alejandro
% Original name: g20220202_AABBA_the_other_plots_and_pipeline_white.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

definput.keyvals.dir_out=[];
[flags,keyvals]  = ltfatarghelper({},definput,varargin);
dir_out_ACI = keyvals.dir_out;

if nargin == 0
    clc, close all   
end

h = []; hname = [];
outs = [];

if ~strcmp(experiment,'speechACI_Logatome-abda-S43M')
    error('This script is only valid for experiment ''speechACI_Logatome-abda-S43M''');
end

N_subjects = 1;

if ~strcmp(masker,'white')
    error('This script is only valid for ''white'' backgrounds');
end
N_maskers = 1;

if N_maskers <= 3
    maskercolors = {'b','r',rgb('Green')};
end

TF_type = 'gammatone';
expvar_after_reversal = 4; 
% N_trials = 4000; idx_trialselect = 1:N_trials; warning('Temporal')
glmfct = 'l1glm'; % former 'lassoglmslow'

switch glmfct
    case 'l1glm'
        N_lambda = 30;
        Lambdas = logspace(-4, -1, N_lambda);
        idx = find(Lambdas >= 10^-3);
        Lambdas = Lambdas(idx);
end

dir_subject = [fastACI_dir_data experiment filesep subject filesep];

dir_suff = '-AABBA-2022-01';
switch subject
    case 'SLV'
        fname_results = 'savegame_2021_11_19_12_37_SLV_speechACI_Logatome-abda-S43M_white.mat';
    % case 'dau1997'
    %     fname_results = 'savegame_2022_01_16_17_30_dau1997_speechACI_Logatome-abda-S43M_white.mat';
    %     % dir_suff = 'Run-1'; % looks for Run-1
    % case 'king2019'
    %     fname_results = 'savegame_2022_01_17_19_37_king2019_speechACI_Logatome-abda-S43M_white.mat';
    %     % dir_suff = 'Run-1';
    % case 'relanoiborra2019'
    %     fname_results = 'savegame_2022_01_17_02_55_relanoiborra2019_speechACI_Logatome-abda-S43M_white.mat';
    %     % dir_suff = 'Run-1';
    % case 'osses2021'
    %     fname_results = 'savegame_2021_12_09_00_11_osses2021_speechACI_Logatome-abda-S43M_white.mat';
    %     % dir_suff = 'Run-6-opt-calibrated';
    % case 'osses2022a'
    %     fname_results = 'savegame_2021_12_13_17_04_osses2022a_speechACI_Logatome-abda-S43M_white.mat';
    %     % dir_suff = 'Run-1';
end

% loading data
switch subject
    case {'osses2021','osses2022a','dau1997','king2019','maxwell2020','relanoiborra2019'}
        % Case for simulations:
        % dirs{1} = [dir_subject dir_suff filesep];
        dirs{1} = [dir_subject filesep];
        dir_results = [dirs{1} 'Results' dir_suff filesep];
        
        fname_results = Get_filenames(dir_results,'savegame*white*.mat');
        if length(fname_results) == 1
            fname_results = fname_results{1};
        else
            error('Check what happened')
        end

    otherwise
        %%% Assuming that the participant is a real listener:
        dir_results = [dir_subject 'Results' filesep];

end

[cfg_game,data_passation] = Convert_ACI_data_type([dir_results fname_results]);
cd(dir_results)

%%%
subj = Get_subjectname_from_dirname(cfg_game.dir_noise);
if strcmp(subj,'SLeo')
    error('Continue debugging here...')
    % [xx,cfg_game.dir_noise] = Get_subjectname_from_dirname(cfg_game.dir_noise,'SLV');
    % [xx,cfg_game.dir_target] = Get_subjectname_from_dirname(cfg_game.dir_target,'SLV');
    % if exist(cfg_game.dir_noise,'dir')
    %     copyfile([dir_results fname_results],[dir_results fname_results(1:end-4) '-old-dirname.mat']);
    %     save([dir_results fname_results],'cfg_game','data_passation');
    % else
    %     % Then, the noises were not originally generated in this computer:
    %     cfg_game = Check_cfg_crea_dirs(cfg_game);
    % 
    %     % Now we try again to store the new location:
    %     if exist(cfg_game.dir_noise,'dir')
    %         copyfile([dir_results fname_results],[dir_results fname_results(1:end-4) '-old-dirname.mat']);
    %         save([dir_results fname_results],'cfg_game','data_passation');
    %     end
    % end
    % fprintf('The subject name of the target/noise directories was changed from SLeo to SLV...\n');
end
%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare ACI analysis
cfg_game = Check_cfg_crea_dirs(cfg_game); % Updates dir_target/dir_noise to local folders
Data_matrix = [];

flags_for_input = {TF_type, ...
    'trialtype_analysis', 'total',...
    'N_folds', 10, ...
    'dir_noise',cfg_game.dir_noise, ...
    'dir_target',cfg_game.dir_target, ...
    'add_signal',0, ...
    'apply_SNR',0, ...
    'skip_if_on_disk',1, ...
    'no_permutation', ...
    'no_bias', ...
    'no_plot', ... % 'idx_trialselect',idx_trialselect, ...
    'expvar_after_reversal',expvar_after_reversal, ...
    'lambda',Lambdas, ...
    'pyramid_script','imresize', ...
    'pyramid_shape',0, ...
    'dir_out',dir_out_ACI ...
    };

[ACI,cfg_ACI,results, Data_matrix] = fastACI_getACI(fname_results, glmfct, flags_for_input{:}, 'Data_matrix', Data_matrix);

outs.ACI = ACI;
outs.cfg_ACI = cfg_ACI;
outs.results = results;
outs.cfg_game = cfg_game;
outs.data_passation = data_passation;
%%% End ACI analysis