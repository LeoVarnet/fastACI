function [str_stim,data_passation] = replication_ahumada1975_user(cfg,data_passation)
% function [str_stim,data_passation] = replication_ahumada1975_user(cfg,data_passation)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

str_stim = [];
%%% Reading info from data_passation:
i_current = data_passation.i_current;
n_stim    = data_passation.n_stim(i_current);
expvar    = data_passation.expvar(i_current); % It is SNR here
 
%%% Defining whether the current trial is 'target1' or 'target2'
n_signal  = cfg.n_targets_sorted(n_stim); 

if n_signal == 1
    signal = cfg.signal;
else
    signal = 0*cfg.signal;
end
%%%
fname_noise = [cfg.dir_noise cfg.ListStim(n_stim).name];
 
noise = audioread(fname_noise);
 
% Filtering between 20 Hz and 4000 Hz -- Kron-Hite filter replaced by a 4th order Butterworth
[B,A] = butter(4,[20 4000]/(cfg.fs/2));
noise = filter(B,A,noise);

gain_snr = 10^(expvar/20);
signal = gain_snr * signal;

tuser_cal = noise+signal;
 
str_stim.tuser            = tuser_cal;
str_stim.stim_tone_alone  = signal;
