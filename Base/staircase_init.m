function [str_inout,cfg] = staircase_init(str_inout,cfg)
% function [str_inout,cfg] = staircase_init(str_inout,cfg)
%
% Initialise a new staircase
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
debut_i  = str_inout.debut_i;

response=[];
n_correctinarow = 0;
expvar = cfg.startvar;

i_current = debut_i;
if isfield(cfg,'start_stepsize')
    stepsize = cfg.start_stepsize;
else
    % This is the case when a constant stimulus procedure is being run
end
isbreak = 0;

if ~isfield(cfg,'step_resolution')
    cfg.step_resolution = 'linear';
end

str_inout.response = response;
str_inout.n_correctinarow = n_correctinarow;
str_inout.expvar = expvar;
str_inout.i_current = i_current;
if isfield(cfg,'start_stepsize')
    str_inout.stepsize = stepsize;
end
str_inout.isbreak = isbreak;

str_inout.reversal_current = 0;
