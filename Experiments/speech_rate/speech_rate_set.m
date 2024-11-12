function cfg_inout = speech_rate_set(cfg_inout)
% function cfg_out = speechLAMIv2_set(cfg_inout)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global global_vars

if nargin == 0
    cfg_inout = [];
end

dir_data_experiment = [fastACI_paths('dir_data') cfg_inout.experiment_full filesep];

% dir_logatome_src = '/media/alejandro/My Passport/Databases/data-Speech/french/Logatome/'; % Logatome-average-power-speaker-S46M_FR.mat
dir_target = [dir_data_experiment cfg_inout.Subject_ID filesep 'speech-samples' filesep];
dir_noise  = [dir_data_experiment cfg_inout.Subject_ID filesep 'Stim-processed' filesep];

%%% Parameters to create the targets:
cfg.fs        = 48000; % Hz, sampling frequency
cfg.dur_ramp  = 1e-1;%75e-3; % cosine ramp
cfg.SPL       = 65; % target level, by default level of the noise (the speech 
                    % level is adapted)
%%%
dBFS = 100; % full scale value to store the waveforms
cfg.dBFS = dBFS;
%%%

%cfg.bRove_level = 0; % New option as of 16/04/2021
%cfg.Rove_range  = 2.5; % plus/minus this value, changed from 4 to 2.5 dB on 26/05/2021

if isfield(global_vars,'Language')
    Language = global_vars.Language;
else
    Language = 'FR'; % French by default
end
cfg.Language = Language; % or 'EN'

cfg.N_presentation = 300; % number of stimuli / condition
cfg.N_target  = 1;     % Number of conditions
cfg.N         = cfg.N_target*cfg.N_presentation;

cfg_inout.dir_data_experiment = dir_data_experiment;

% Change the following names:
if isfield(cfg_inout,'dir_target')
    if ~exist(cfg_inout.dir_target,'dir')
        if iswindows
            warning('Leo: I had to re-enable this option here to have the possibility to use noises from other subjects, but I remember that this was an error source for you...')
            disp('Pausing for 10 s... (this is a temporal message)')
            pause(10)
        end
        cfg_inout.dir_target = dir_target;
    end
else
    cfg_inout.dir_target = dir_target;
end

if isfield(cfg_inout,'dir_noise')
    if ~exist(cfg_inout.dir_noise,'dir')
        cfg_inout.dir_noise  = dir_noise;
    end
else
    cfg_inout.dir_noise  = dir_noise;
end
cfg_inout.noise_type = '';
 
cfg_inout = Merge_structs(cfg,cfg_inout);
