function [str_stim,data_passation] = segmentation_user(cfg,data_passation)
% function [str_stim,data_passation] = speechLAMIv2_user(cfg,data_passation)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

str_stim = [];

i_current = data_passation.i_current;
% SNR       = data_passation.expvar(i_current);

n_stim    = data_passation.n_stim(i_current);
% n_signal  = cfg.n_targets_sorted(n_stim); % Not needed

if isfield(cfg,'Rove_level')
    % Rove_range has to be adjusted in the _set.m script.
    data_passation.Rove_level(i_current) = cfg.Rove_level(n_stim);
    presentation_gain = 10^(data_passation.Rove_level(i_current)/20);
else
    presentation_gain = 1;
end

fname_noise = [cfg.dir_noise cfg.ListStim(n_stim).name];
noise = audioread(fname_noise);
 
tuser_cal = noise;
 
str_stim.tuser            = presentation_gain*tuser_cal;
str_stim.stim_tone_alone  = str_stim.tuser;
