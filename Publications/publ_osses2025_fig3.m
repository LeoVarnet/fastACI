function [] = publ_osses2025_fig3()
% function publ_osses2025_fig3()
%%%%% Figure 3 from fastACI paper %%%%%
%
% % To display Fig. 3 of Osses et al use
%     publ_osses2025_fig3();
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear all 

SNRprctile = [5 95]; 

h = [];
hname = [];

% cho ose the data you want to analyse
% 
% experiment = 'speechACI_varnet2013';
% participant = 'SLeo';%'osses2021';%
% opts_ACI.DimCI = 'gammatone'; %'gammatone';% 'lyon';%
% masker = 'white';
% dir_participant = [fastACI_paths('dir_data') experiment filesep participant filesep];
% Ntrials2load = 10000;
experiment = 'speechACI_Logatome-abda-S43M';
participant = 'S04';%'osses2021';%
opts_ACI.DimCI = 'gammatone'; %'gammatone';% 'lyon';%
masker = 'bumpv1p2_10dB';
dir_participant = [fastACI_paths('dir_data') experiment filesep participant filesep];
Ntrials2load = 10000;


cd([dir_participant 'Results' filesep])
D = dir([dir_participant 'Results' filesep 'savegame_*' masker '.mat']);% fname_results = D(end).name;
fname_results = D(end).name; folder_results = D(end).folder;
load([folder_results filesep fname_results])

if exist('data_passation')
    SNR = data_passation.expvar;
elseif exist('ListStim')
    SNR = vertcat(ListStim.RSB);
end
clear D

% cfg_game.dir_target = [dir_participant 'speech-samples' filesep];
% cfg_game.dir_noise = [dir_participant 'NoiseStim' filesep];

Data_matrix = [];

%% Options for ACI estimation

flags_for_input = {...
    'trialtype_analysis','total', ...
    'no_permutation', ...
    'no_bias', ...
    'no_plot'...
    };
  %  'N_folds', 10, ...
  %  'add_signal',0, ...
  %  'apply_SNR',0, ...
%'idx_trialselect', 1:Ntrials2load, ...
    %'dir_noise',cfg_game.dir_noise, ...
    %'dir_target',cfg_game.dir_target, ...
% 'perc', [5 95]
%'expvar_limits',[prctile(SNR,SNRprctile(1)) prctile(SNR,SNRprctile(2))], ... %'dir_out','/Users/leovarnet/Documents/ACIs/', ...
    %
    %     'f_limits', [0 8000],...[0 4050], ...
%     'spect_Nwindow', 512, ...
%     'spect_overlap', 1/2, ...
%     'spect_NFFT', 512/2,...3*512, ...
%     'spect_unit', 'linear', ...'dB', ... 

%% Calculate and plot revcorr ACI - classic_revcorr

%flags_for_input(1:2) = {'idx_trialselect', 1:10000};

type = 'correlation';
[ACI,cfg_ACI,results,Data_matrix] = fastACI_getACI(fname_results, opts_ACI.DimCI, type, flags_for_input{:},'Data_matrix',Data_matrix);

% display ACI
figure('Position', [100 100 400 300])
affichage_tf(ACI,'CI', 'cfg', cfg_ACI); hold on
%affichage_tf_add_Praat_metrics(cfg_game.dir_target,cfg_ACI,[], {'-','-.'},{[0.6,0,0],[0,0,0.6]},1.5);
title(['correlation'], 'interpreter','none')

%% Calculate and plot revcorr ACI - weighted sum

flags_for_input(1:2) = {'idx_trialselect', 1:10000};

type = 'weighted_sum';
[ACI,cfg_ACI,results,Data_matrix] = fastACI_getACI(fname_results, opts_ACI.DimCI, type, flags_for_input{:},'Data_matrix',Data_matrix);

% display ACI
figure('Position', [100 100 400 300])
affichage_tf(ACI,'CI', 'cfg', cfg_ACI); hold on
%affichage_tf_add_Praat_metrics(cfg_game.dir_target,cfg_ACI,[], {'-','-.'},{[0.6,0,0],[0,0,0.6]},1.5);
title(['weighted sum'], 'interpreter','none')
% 
% %% Calculate and plot revcorr ACI - weighted sum
% 
% flags_for_input(1:2) = {'idx_trialselect', 1:10000};
% 
% type = 'glm';
% [ACI,cfg_ACI,results,Data_matrix] = fastACI_getACI(fname_results, opts_ACI.DimCI, type, flags_for_input{:},'Data_matrix',Data_matrix);
% 
% % display ACI
% figure('Position', [100 100 400 300])
% affichage_tf(ACI,'CI', 'cfg', cfg_ACI); hold on
% %affichage_tf_add_Praat_metrics(cfg_game.dir_target,cfg_ACI,[], {'-','-.'},{[0.6,0,0],[0,0,0.6]},1.5);
% title(['glm'], 'interpreter','none')

%% Calculate and plot revcorr ACI - weighted sum without zscore

type = 'weighted_sum';
[ACI,cfg_ACI,results,Data_matrix] = fastACI_getACI(fname_results, opts_ACI.DimCI, type, flags_for_input{:},'zscore',0,'Data_matrix',Data_matrix,'skip_if_on_disk',0);

% display ACI
figure('Position', [100 100 400 300])
affichage_tf(ACI,'CI', 'cfg', cfg_ACI); hold on
%affichage_tf_add_Praat_metrics(cfg_game.dir_target,cfg_ACI,[], {'-','-.'},{[0.6,0,0],[0,0,0.6]},1.5);
title(['weighted sum, no z-scoring'], 'interpreter','none')

%% Calculate and plot revcorr ACI - glmfitqp

type = 'glmfitqp';
[ACI,cfg_ACI,results,Data_matrix] = fastACI_getACI(fname_results, opts_ACI.DimCI, type, flags_for_input{:},'Data_matrix',Data_matrix);

% display ACI
figure('Position', [100 100 400 300])
affichage_tf(ACI,'CI', 'cfg', cfg_ACI); hold on
%affichage_tf_add_Praat_metrics(cfg_game.dir_target,cfg_ACI,[], {'-','-.'},{[0.6,0,0],[0,0,0.6]},1.5);
title(['ACI ' type], 'interpreter','none')

%% Calculate and plot revcorr ACI - lasso

type = 'lasso';
[ACI,cfg_ACI,results,Data_matrix] = fastACI_getACI(fname_results, opts_ACI.DimCI, type, flags_for_input{:},'Data_matrix',Data_matrix);

% display ACI
figure('Position', [100 100 400 300])
affichage_tf(ACI,'CI', 'cfg', cfg_ACI); hold on
%affichage_tf_add_Praat_metrics(cfg_game.dir_target,cfg_ACI,[], {'-','-.'},{[0.6,0,0],[0,0,0.6]},1.5);
title(['ACI ' type], 'interpreter','none')
end