function cfg_inout = toneinnoise_ahumada1975_set(cfg_inout)
% function cfg_inout = toneinnoise_ahumada1975_set(cfg_inout)
%
% Function comparable to *_set.m functions from AFC toolbox
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    cfg_inout = [];
end
 
cfg = [];
  
cfg.filename_target = {'tone-500-Hz.wav'};
   
%%% Parameters to create the targets:
%     (c) = compulsory for this experiment
cfg.N_target   = 2; % target present, target absent
cfg.N_presentation = 400;  % (c) number of stimuli / condition
 
%%% Sampling frequency should be compulsory if targets are not on disk:
cfg.fs         = 10000; % (c) Hz, sampling frequency 
cfg.f_target   = 500; % Hz
cfg.dur_target = 100e-3; 
cfg.dur_noise  = 500e-3;
cfg.ramp_dur_target = 10e-3; % s, up and down ramp ('fade in and fade out')
cfg.ramp_dur_noise  = 50e-3; % s, up and down ramp ('fade in and fade out')
cfg.lvl_noise  = 70;  % dB, level for the background noises, 85 dB is the level used by Ahumada
cfg.lvl_target = 70;  % 
cfg.dBFS       = 100; % (c) 93.61; % Calibration with Headphones HD 650
  
t = (1:round(cfg.dur_target*cfg.fs))/cfg.fs;
signal = sin(2*pi*cfg.f_target*t(:));
signal = scaletodbspl(signal,cfg.lvl_target,cfg.dBFS);

rp = cos_ramp(length(signal),cfg.fs,cfg.ramp_dur_target*1000,cfg.ramp_dur_target*1000);
signal = signal.* rp(:);

N_sil = round(cfg.fs*(cfg.dur_noise-cfg.dur_target)/2);
sil = zeros([N_sil 1]);

cfg.signal = [sil; signal; sil];

cfg_inout = Merge_structs(cfg,cfg_inout);
