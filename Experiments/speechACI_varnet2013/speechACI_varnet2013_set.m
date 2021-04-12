function cfg_inout = speechACI_varnet2013_set(cfg_inout)
% function cfg_out = speechACI_varnet2013_set(cfg_inout)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    cfg_inout = [];
end
 
dir_main = '/home/alejandro/Documents/Databases/data/fastACI/speechACI_varnet2013/';
dir_speech = [dir_main cfg_inout.Subject_ID filesep 'speech-samples' filesep];
if ~isfield(cfg_inout,'Condition')
    dir_name_noise = 'NoiseStim';
else
    switch lower(cfg_inout.Condition) % lower case
        case 'white'
            dir_name_noise = 'NoiseStim';
        case 'ssn'
            dir_name_noise = 'NoiseStim-SSN';
        otherwise
            error('%s: Condition not recognised',upper(mfilename));
    end
end
dir_noise  = [dir_main cfg_inout.Subject_ID filesep dir_name_noise filesep];

dBFS = 100;

%%% Parameters to create the targets:
cfg.fs        = 44100; % Hz, sampling frequency
cfg.dur_ramp  = 75e-3; % cosine ramp
cfg.SPL       = 65; % target level, by default level of the noise (the speech 
                    % level is adapted)
cfg.dBFS      = dBFS;

% Change the following names:
cfg.N_presentation = 2500; % 5000;  % number of stimuli / condition
cfg.N_target  = 2;     % Number of conditions
cfg.N         = cfg.N_target*cfg.N_presentation;

cfg_inout.dir_speech = dir_speech;
cfg_inout.dir_noise  = dir_noise;

% TODO: change this name:
% cfg_inout.dir_stim  = dir_noise;
 
cfg_inout = Merge_structs(cfg,cfg_inout);
