function cfg_inout = toneinnoise_ahumada1975_cfg(cfg_inout)
% function str_stim = toneinnoise_ahumada1975_cfg(cfg_in)
%
% Function comparable to *_cfg.m functions from AFC toolbox
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    cfg_inout = [];
end

%% Uncomment/comment the following fields to customise this experiment:
cfg_inout.response_names = {'ton dans le bruit','bruit seul'};%{'Tone present', 'noise only'}; % arbitrary string, make sure it coincides with filename_target
cfg_inout.Language = 'FR';%'EN'; % 
cfg_inout.feedback = 1; % 0
cfg_inout.sessionsN = 400; % arbitrary number
cfg_inout.warmup         = 0; % 'oui', CAUTION: Overwritten in the case of simulation

% cfg_in.randorder = 1; 1 is the default

%% Compulsory options:
cfg_inout.adapt = 'weighted-up-down'; % 'weighted-up-down'=2 or 1='transformed up-down'
% 1. Compulsory fields for cfg_in.adapt:
cfg_inout.startvar = 0;  % 0 dB SNR
cfg_inout.expvar_description = 'SNR (dB)'; % Optional, but related to startvar

% 2. Optional fields (if omitted, the defaults will be assigned later):
if cfg_inout.bCrea == 1
    % Only loaded during creation:
    [cfg_inout,txt2show] = staircase_defaults(cfg_inout);
    if ~isempty(txt2show)
        fprintf('\t%s.m: List of assigned fields for the staircase procedure:',mfilename);
        Show_cell(txt2show);
        disp('Pausing for 5 seconds...');
        pause(5);
    end
end
 
%% Otherwise, the defaults will be adopted:
if ~isfield(cfg_inout,'response_names')
    cfg_inout.response_names = cfg_inout.filename_target;
end

if ~isfield(cfg_inout,'Language')
    cfg_inout.Language = cfg_inout.keyvals.Language;
end

if ~isfield(cfg_inout,'feedback')
    cfg_inout.feedback = cfg_inout.keyvals.feedback;
end
if ~isfield(cfg_inout,'sessionsN')
    cfg_inout.sessionsN = cfg_inout.keyvals.sessionsN; % CAUTION: Overwritten in the case of simulation
end

if ~isfield(cfg_inout,'randorder')
    cfg_inout.randorder = cfg_inout.keyvals.randorder;
end
