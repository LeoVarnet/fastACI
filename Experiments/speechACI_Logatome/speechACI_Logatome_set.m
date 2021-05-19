function cfg_inout = speechACI_Logatome_set(cfg_inout)
% function cfg_out = speechACI_Logatome_set(cfg_inout)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    cfg_inout = [];
end

if ismac % lab's computer
    dir_main = '/Users/leovarnet/ownCloud/Data/Projet fastACI/';
    % dir_main = 'C:\Users\Varnet Léo\Dropbox\Data\Projet fastACI\';
elseif isunix % Alejandro's computer
    dir_main = '/home/alejandro/Documents/Databases/data/fastACI/speechACI_Logatome/';
elseif ispc % Leo's computer
    dir_main = 'C:\Users\Léo\ownCloud\Data\Projet fastACI/';%'C:\Users\Varnet Lï¿½o\Dropbox\Data\Projet fastACI/';
    % dir_main = 'C:\Users\Varnet Léo\Dropbox\Data\Projet fastACI\';
end

dir_logatome_src = '/media/alejandro/My Passport/Databases/data-Speech/french/Logatome/'; % Logatome-average-power-speaker-S46M_FR.mat

dir_speech = [dir_main cfg_inout.Subject_ID filesep 'speech-samples' filesep];

if ~isfield(cfg_inout,'Condition')
    cfg_inout.Condition = 'SSN'; % by default it is speech-shaped noise
end
    
% cfg_inout.Condition = noise_type;

switch lower(cfg_inout.Condition) % lower case
    case 'white'
        dir_name_noise = 'NoiseStim-white';
        noise_type = 'white';
    case 'pink'
        dir_name_noise = 'NoiseStim-pink';
        noise_type = 'pink';
    case 'ssn'
        dir_name_noise = 'NoiseStim-SSN';
        noise_type = 'SSN';
    otherwise
        error('%s: Condition not recognised',upper(mfilename));
end
dir_noise  = [dir_main cfg_inout.Subject_ID filesep dir_name_noise filesep];

dBFS = 100;
 
%%% Parameters to create the targets:
cfg.fs        = 16000; % Hz, sampling frequency
cfg.dur_ramp  = 75e-3; % cosine ramp
cfg.SPL       = 65; % target level, by default level of the noise (the speech 
                    % level is adapted)
cfg.dBFS      = dBFS;

cfg.bRove_level = 1; % New option as of 16/04/2021
cfg.Rove_range  = 4; % plus/minus this value

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
