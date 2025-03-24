function cfg_out = speech_rate_cfg(cfg_in)
% function str_stim = speechLAMIv2_cfg(cfg_in)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    cfg_in = [];
end

cfg = [];
cfg_out = cfg_in; % copying input to output struct

cfg.filename_target = Get_filenames(cfg_in.dir_target,'*.wav'); %Get_filenames(cfg_in.dir_target,[cfg_in.Cond_extra_2 '*.wav']);

if ~isempty(cfg.filename_target)
    for i = 1:length(cfg.filename_target)
        ftarget = strsplit(cfg.filename_target{i}(1:end-4),'_'); % excludes the '.wav'
        cfg.response_names{i} = ftarget{end}; % second syllable
        % cfg.response_names   = {'pa', 'ta'}; 
    end
end
cfg.response_correct_target = [1,2];
cfg.warmup         = 1; % 'oui', CAUTION: Overwritten in the case of simulation
 
% cfg.Language       = 'FR'; % or 'EN'

bDebug = 1; 
cfg.displayN       = 1; % 'oui'
cfg.feedback       = 0;
 
cfg.sessionsN      = 400; % CAUTION: Overwritten in the case of simulation
cfg_out.randorder  = 1;
  
cfg.startvar = 0;  % old name 'm_start'

%%% IDLE variables:
cfg.adapt    = 2; % warning('default is 1 - now we are only testing')    %  1; warning('adapt type temporarily set to 1\n');%      
cfg.expvar_description = 'Idle expvar, this experimental variable is not used during the experiment';
cfg.step_resolution = 'linear';
cfg.start_stepsize = 2; % dB
cfg.min_stepsize = 1; % dB
cfg.adapt_stepsize = 50/100;

% Staircase algorithm parameters
switch cfg.adapt
    case {2, 'weighted-up-down'}
        % IDLE
        cfg.rule = [1 1]; 
        target_score = .707;
        cfg.step_down  = 1;
        cfg.step_up    = cfg.step_down*target_score/(1-target_score); % Kaernbach1991, Eq. 1
        cfg.min_stepsize = 1/cfg.step_up;
end
%%% End IDLE variables

cfg_out = Merge_structs(cfg,cfg_out);
