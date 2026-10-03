function cfg_out = major2026_cfg(cfg_in)
% function str_stim = major2026_cfg(cfg_in)
%
% Function comparable to *_cfg.m functions from AFC toolbox
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    cfg_in = [];
end

cfg = [];
cfg_out = cfg_in; % copying input to output struct

%cfg.filename_target = Get_filenames(cfg_in.dir_target,[cfg_in.Cond_extra_2 '*.wav']);
cfg.filename_target  = {'v1_close_chi.wav','v1_close_ci.wav','v1_far_chi.wav','v1_far_ci.wav','v2_close_chi.wav','v2_close_ci.wav','v2_far_chi.wav','v2_far_ci.wav'};
cfg.response_names = {'bifaci','bifachi'};
cfg.response_correct_target = [2,1,2,1,2,1,2,1];
cfg.correctness_matrix = [1,2];

cfg.warmup         = 1; % 'oui'

bDebug = 1; 
cfg.displayN       = bDebug; % 'oui'
cfg.feedback       = 1;

cfg.sessionsN      = 400; 
cfg_out.randorder  = 1;

cfg.startvar = 0;  % old name 'm_start'
cfg.adapt = 'weighted-up-down'; % or 2
if isfield(cfg_in,'Condition')
    cfg.expvar_description = 'SNR (dB)';
    cfg.step_resolution = 'linear';
    cfg.start_stepsize = 2; % dB
    cfg.maxvar = 10;
end

cfg.adapt_stepsize = 50/100;

% Staircase algorithm parameters
switch cfg.adapt
    case {1, 'transformed-up-down'}
        cfg.rule = [1 2]; % [up down]-rule: [1 2] = 1-up 2-down   
        cfg.target_score = .707; % fixed value here
        cfg.step_up    = 1;
        cfg.step_down  = 1;
        cfg.min_stepsize = 1; % dB

    case {2, 'weighted-up-down'}
        cfg.rule = [1 1]; 
        cfg.target_score = .707;
        cfg.step_down  = 1;
        cfg.step_up    = cfg.step_down*cfg.target_score/(1-cfg.target_score); % Kaernbach1991, Eq. 1
        cfg.min_stepsize = 1/cfg.step_up;

    otherwise
        error('Not validated yet...')
end

cfg_out = Merge_structs(cfg,cfg_out);
