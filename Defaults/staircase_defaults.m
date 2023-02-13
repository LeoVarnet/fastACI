function [cfg, txt_added_fields] = staircase_defaults(cfg)
% function [cfg, txt_added_fields] = staircase_defaults(cfg)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% debut_i = disp(''); % str_inout.debut_i;

txt_added_fields = [];
if ~isfield(cfg,'start_stepsize')
    cfg.start_stepsize = 1; % dB
    txt_added_fields{end+1} = sprintf('Default: start_stepsize=%.1f',cfg.start_stepsize);
end

if ~isfield(cfg,'target_score')
    cfg.target_score = .707;
    txt_added_fields{end+1} = sprintf('Default: target_score=%.3f',cfg.target_score);
end
target_score = cfg.target_score;

switch cfg.adapt 
    case {1, 'transformed-up-down'}
        % cfg.step_resolution = 'linear';
        switch target_score
            case .707
                rule = [1 2];
            otherwise
                error('Add manually to %s the new target score with its underlying rule',mfilename);
        end
        step_up = 1;
        step_down = 1;
        
        cfg.step_down = 1;
        adapt_stepsize = 0.9; % factor. 1=means no adapted step size
        min_stepsize = 0.5; % dB
        
    case {2, 'weighted-up-down'}
        rule = [1 1]; 
        step_down = 1;
        step_up = step_down*target_score/(1-target_score); % Kaernbach1991, Eq. 1
        
        adapt_stepsize = 0.5;
        min_stepsize = 1/step_up;
    otherwise
        error('Not validated yet...')
end

if ~isfield(cfg,'rule')
    % Assigning the new rule
    cfg.rule = rule; % [up down]-rule: [1 2] = 1-up 2-down   
    txt_added_fields{end+1} = sprintf('Default: rule=[%.1f %.1f]',cfg.rule);
else
    % Otherwise it checks whether the rule coincides with the target_score:
    if cfg.rule(1) ~= rule(1)
        error('a ''rule'' is specified that does not coincide with the requested target_score');
    end
    if cfg.rule(2) ~= rule(2)
        error('a ''rule'' is specified that does not coincide with the requested target_score');
    end
end

if ~isfield(cfg,'step_up')
    cfg.step_up = step_up;
    txt_added_fields{end+1} = sprintf('Default: step_up=%.2f',cfg.step_up);
else
    if cfg.step_up ~= step_up
        error('Wrong step up');
    end
end

if ~isfield(cfg,'step_down')
    cfg.step_down = step_down;
    txt_added_fields{end+1} = sprintf('Default: step_down=%.2f',cfg.step_down);
else
    if cfg.step_down ~= step_down
        error('Wrong step down');
    end
end

if ~isfield(cfg,'adapt_stepsize')
    cfg.adapt_stepsize = adapt_stepsize; % factor. 1=means no adapted step size
    txt_added_fields{end+1} = sprintf('Default: adapt_stepsize=%.2f',cfg.adapt_stepsize);
end

if ~isfield(cfg,'min_stepsize')
    cfg.min_stepsize = min_stepsize; 
    txt_added_fields{end+1} = sprintf('Default: min_stepsize=%.2f',cfg.min_stepsize);
end