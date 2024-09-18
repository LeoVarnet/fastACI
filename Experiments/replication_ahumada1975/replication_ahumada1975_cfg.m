function cfg_inout = replication_ahumada1975_cfg(cfg_inout)
% function str_stim = replication_ahumada1975_cfg(cfg_in)
%
% Function comparable to *_cfg.m functions from AFC toolbox
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    cfg_inout = [];
end

cfg_inout.Language = 'FR';%'EN'; % 
switch cfg_inout.Language
    case 'FR'
        cfg_inout.response_names = {'ton dans le bruit','bruit seul'}; % arbitrary string, make sure it coincides with filename_target
    case 'EN'
        cfg_inout.response_names = {'Tone present', 'Tone absent'}; % arbitrary string, make sure it coincides with filename_target
end

cfg_inout.feedback = 1;
cfg_inout.sessionsN = 400; % arbitrary number

cfg_inout.adapt = 0; % No adaptive staircase
cfg_inout.probe_periodicity = 10; % insert a probe (easy) stimulus every Nth trial

cfg_inout.warmup = 0; % No warm-up phase for this particular experiment

cfg_inout.startvar = -21.19;% -15;% so that E/N_0 = 11.8 dB % % dB SNR
cfg_inout.expvar_description = 'SNR (dB)'; % Optional, but related to startvar
 
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
