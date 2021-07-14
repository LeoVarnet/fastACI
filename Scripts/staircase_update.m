function [str_inout,cfg] = staircase_update(str_inout,cfg)
% function [str_inout,cfg] = staircase_update(str_inout,cfg)
%
% staircase update phase
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iscorrect       = str_inout.iscorrect;
n_correctinarow = str_inout.n_correctinarow;
stepsize        = str_inout.stepsize;
expvar          = str_inout.expvar;
reversal_current = str_inout.reversal_current;
if isfield(str_inout,'staircase_direction')
    staircase_direction = str_inout.staircase_direction;
else
    staircase_direction = [];
end

n_up   = cfg.rule(1);
n_down = cfg.rule(2);

if n_up ~= 1
    error('%s: Only one-up ''n''-down adaptation rules have been validated...',upper(mfilename))
end
% 2-down 1-up staircase method on m
if ~iscorrect
    switch cfg.step_resolution
        case 'linear'
            expvar = expvar + stepsize*cfg.step_up;
        case 'octave'
            if expvar == 0
                expvar = -1;
            else
                expvar = sign(expvar)*abs(expvar)*2^(-cfg.step_up);
            end
    end
    
    if ~isempty(staircase_direction)
        if strcmp(staircase_direction,'down') 
            % Then this is a turning point from correct to incorrect
            reversal_current = reversal_current+1;
            staircase_direction = 'up';
        end
    else
        % Then we define the direction
        staircase_direction = 'up';
    end
elseif n_correctinarow == n_down
    switch cfg.step_resolution
        case 'linear'
            expvar = expvar - stepsize*cfg.step_down;
        
        case 'octave'
            if expvar == 0
                expvar = -1;
            else
                expvar = sign(expvar)*abs(expvar)*2^(cfg.step_down);
            end
    end
    n_correctinarow=0;
    
    if reversal_current == 0
        % Then we define the direction
        staircase_direction = 'down';
    else
        if strcmp(staircase_direction,'up') 
            % Then this is a turning point from incorrect to correct
            reversal_current = reversal_current+1;
            staircase_direction = 'down';
        end
    end
elseif n_correctinarow == 1
    % Nothing to do, except to assign the staircase_direction for the first trial
    if reversal_current == 0
        % Then we define the direction
        staircase_direction = 'down';
    end
end

% if expvar>1
%     error('This line should be removed')
%     % expvar=1;
% end

if reversal_current ~= 0 && mod(reversal_current,2) == 0
    % if n_correctinarow == 0 % change stepsize
    if stepsize>cfg.min_stepsize
        stepsize = stepsize*cfg.adapt_stepsize;
        stepsize = max(stepsize,cfg.min_stepsize);
    else
        stepsize = cfg.min_stepsize;
    end
end
str_inout.n_correctinarow = n_correctinarow;
str_inout.stepsize        = stepsize;
str_inout.expvar          = expvar;

str_inout.reversal_current = reversal_current;
str_inout.staircase_direction = staircase_direction;
