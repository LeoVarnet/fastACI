function [str_stim,data_passation] = segmentation_user(cfg,data_passation)
% function [str_stim,data_passation] = speechLAMIv2_user(cfg,data_passation)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

str_stim = [];

i_current = data_passation.i_current;
% SNR       = data_passation.expvar(i_current);

n_stim1    = data_passation.n_stim(i_current);
n_stim2    = data_passation.n_stim(i_current)+cfg.N/2;
% n_signal  = cfg.n_targets_sorted(n_stim); % Not needed

fname_noise1 = [cfg.dir_noise cfg.ListStim(n_stim1).name];
fname_noise2 = [cfg.dir_noise cfg.ListStim(n_stim2).name];

noise1 = audioread(fname_noise1);
noise2 = audioread(fname_noise2);
  
str_stim.tuser            = noise1;
str_stim.tref             = noise2;
str_stim.stim_tone_alone  = str_stim.tuser;
