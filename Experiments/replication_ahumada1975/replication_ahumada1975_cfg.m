function cfg_inout = replication_ahumada1975_cfg(cfg_inout)
% function str_stim = replication_ahumada1975_cfg(cfg_in)
%
% Function comparable to *_cfg.m functions from AFC toolbox
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    cfg_inout = [];
end

cfg_inout.response_names = {'ton dans le bruit','bruit seul'};%{'Tone certainly not present', 'Tone probably not present', 'Tone probably present', 'Tone certainly present'}; % arbitrary string, make sure it coincides with filename_target
%cfg_inout.correctness_matrix = [2 2 1 1]; % correspondence between Likert scale and actual targets
%cfg_inout.response_names = {'Tone present', 'Tone absent'}; % arbitrary string, make sure it coincides with filename_target
cfg_inout.Language = 'FR';%'EN'; % 
cfg_inout.feedback = 0;
cfg_inout.sessionsN = 400; % arbitrary number
% cfg_in.randorder = 1; 1 is the default

cfg_inout.adapt = 0; % No adaptive staircase
cfg_inout.probe_periodicity = 10; % insert a probe (easy) stimulus every Nth trial

%%
cfg_inout.warmup = 0; % No warm-up phase for this particular experiment

cfg_inout.startvar = -21.19;% so that E/N_0 = 11.8 dB % -15;  % dB SNR
cfg_inout.expvar_description = 'SNR (dB)'; % Optional, but related to startvar

% % 2. Optional fields (if omitted, the defaults will be assigned later):
% if cfg_inout.bCrea == 1
%     % Only loaded during creation:
%     [cfg_inout,txt2show] = staircase_defaults(cfg_inout);
%     if ~isempty(txt2show)
%         fprintf('\t%s.m: List of assigned fields for the staircase procedure:',mfilename);
%         Show_cell(txt2show);
%         disp('Pausing for 5 seconds...');
%         pause(5);
%     end
% end
 
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
