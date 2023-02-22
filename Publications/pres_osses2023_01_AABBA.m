function pres_osses2023_01_AABBA(varargin)
% functiong pres_osses2023_01_AABBA(varargin)
%
% Generates some figures from the presentation titled 'Paris II - Simulating 
%   the perception of soundscapes, speech-, AM-, and FM-sounds' given in the
%   context of the AABBA meeting celebrated in Vienna on 3 and 4 February 2022.
%
% At this moment, this script requires that you already have the simulation
%    results stored locally in your computer.
%
% % To display the figure in page 11 of the final presentation (PDF) use :::
%     % ACI of participant S2 from Osses & Varnet (2021, DAGA):
%     pres_osses2022_02_AABBA('fig_page11');
%
% % To display the top figure in page 13 use :::
%     pres_osses2022_02_AABBA('fig_page13_top'); % ACIs for models and the human listener
%
% % To display the bottom figure in page 13 use :::
%     pres_osses2022_02_AABBA('fig_page13_bottom'); % thresholds and correlation plot
%
% Author: Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if nargin == 0
%     help pres_osses2023_01_AABBA;
%     return
% end

bDo1 = 0;
bDo2 = 0; % Run simulations
bDo3 = 1; % Obtain ACI from simulations

%%% Common parameters:
experiment = 'speechACI_Logatome-apta-S43M';
Condition = 'white'; % default

if bDo1
    %%% Step 1:
    Subject_ID = 'SAABBA_1'; % Human data
    fastACI_experiment(experiment,Subject_ID,Condition);
    %%%
end

if bDo2
    Subject_ID = 'osses2022a'; % Simulation data
    
    in_std = 0; % default
    maxtimelag_ms = 0; % default
    fname_template_suffix = 'white-2023-01-17';
    fg_bias_each_session = 'bias_global'; % default
    
    Calibration_completed = 0;
    
    for i = 1:10
        if Calibration_completed == 0 
            thres_for_bias = [];
            % Calibration starts
            flags_here = {'in_std',in_std,'thres_for_bias',thres_for_bias, ...
                'fname_template_suffix',fname_template_suffix,'maxtimelag_ms',maxtimelag_ms,'bias_global'};
            [~,~,kv] = fastACI_experiment(experiment,Subject_ID,Condition,flags_here{:});
            thres_for_bias = kv.thres_for_bias;
            
            Calibration_completed = 1;
        end

        flags_here = {'in_std',in_std,'thres_for_bias',thres_for_bias, ... % 'thres_for_bias_each_session',thres_for_bias_each_session, ...
            'fname_template_suffix',fname_template_suffix,'maxtimelag_ms',maxtimelag_ms,fg_bias_each_session};
        [cfg_game,data_passation] = fastACI_experiment(experiment,Subject_ID,Condition,flags_here{:});
    end
end

if bDo3
    Subject_ID = 'osses2022a'; % Simulation data
    
    dir_exp = [fastACI_dir_data experiment filesep Subject_ID filesep 'Results' filesep];
    fname_results = Get_filenames(dir_exp,['savegame_*_' Subject_ID '_*' Condition '.mat']);
    fname_results = [dir_exp fname_results{1}];
    
    N_lambda = 30;
    Lambdas = logspace(-4, -1, N_lambda);
    idx = find(Lambdas >= 10^-3);
    Lambdas = Lambdas(idx);
    
    glmfct = 'l1glm'; % L1 ( lasso ) GLM function
    TF_type = 'gammatone'; % Type of T-F conversion
    flags_in = {    'trialtype_analysis','total', ...
                    'N_folds', 10, ...
                    'no_permutation', ...
                    'no_bias', ...
                    'plot', ... % or set to 'no_plot '
                    'pyramid_script','imresize', ...
                    'pyramid_shape',0, ...
                    'lambda', Lambdas};
    [ACI, cfg_ACI, results, Data_matrix] = fastACI_getACI(fname_results,TF_type,glmfct,flags_in{:}) ;
end

%%% Extra simulation:
%% 1. Creates waveforms
% Subject_ID = 'S03';
% experiment = 'localisationILD';
% fastACI_experiment(experiment, Subject_ID);