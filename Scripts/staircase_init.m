function [str_inout,cfg] = staircase_init(str_inout,cfg)
% function [str_inout,cfg] = staircase_init(str_inout,cfg)
%
% Initialise a new staircase
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iswarmup = cfg.warmup;
debut_i  = str_inout.debut_i;

response=[];
n_correctinarow = 0;
if iswarmup
    if isfield(cfg,'maxvar')
        expvar = cfg.maxvar;
    else
        expvar = cfg.startvar;
    end
else
    expvar = cfg.startvar;
end
i_current = debut_i;
stepsize = cfg.start_stepsize;
isbreak = 0;

if ~isfield(cfg,'step_resolution')
    cfg.step_resolution = 'linear';
end

str_inout.response = response;
str_inout.n_correctinarow = n_correctinarow;
str_inout.expvar = expvar;
str_inout.i_current = i_current;
str_inout.stepsize = stepsize;
str_inout.isbreak = isbreak;

str_inout.reversal_current = 0;
