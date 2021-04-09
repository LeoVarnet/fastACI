function cfg_inout = speechACI_Logatome_set(cfg_inout)
% function cfg_out = speechACI_Logatome_set(cfg_inout)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    cfg_inout = [];
end
 
dir_main = '/home/alejandro/Documents/Databases/data/fastACI/speechACI_Logatome/';
dir_logatome_src = '/media/alejandro/My Passport/Databases/data-Speech/french/Logatome/'; % Logatome-average-power-speaker-S46M_FR.mat

dir_speech = [dir_main cfg_inout.Subject_ID filesep 'speech-samples' filesep];
noise_type = 'SSN';
% noise_type = 'white';

switch noise_type
    case 'SSN'
        dir_noise  = [dir_main cfg_inout.Subject_ID filesep 'NoiseStim-SSN'    filesep];
    case 'white'
        dir_noise  = [dir_main cfg_inout.Subject_ID filesep 'NoiseStim-white'  filesep];
end

cfg_inout.Condition = noise_type;

dBFS = 100;
 
%%% Parameters to create the targets:
cfg.fs        = 16000; % Hz, sampling frequency
cfg.dur_ramp  = 75e-3; % cosine ramp
cfg.SPL       = 65; % target level, by default level of the noise (the speech 
                    % level is adapted)
cfg.dBFS      = dBFS;

% Change the following names:
cfg.N_presentation = 2500;  % number of stimuli / condition
cfg.N_target  = 2;     % Number of conditions
cfg.N         = cfg.N_target*cfg.N_presentation;

cfg_inout.dir_main   = dir_main;
cfg_inout.dir_logatome_src = dir_logatome_src;

cfg_inout.dir_speech = dir_speech;
cfg_inout.dir_noise  = dir_noise;
cfg_inout.noise_type = noise_type;
 
cfg_inout = Merge_structs(cfg,cfg_inout);
