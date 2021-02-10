function str_inout = staircase_update(str_inout,cfg)
% function str_inout = staircase_update(str_inout,cfg)
%
% staircase update phase
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iscorrect       = str_inout.iscorrect;
n_correctinarow = str_inout.n_correctinarow;
stepsize        = str_inout.stepsize;
expvar          = str_inout.expvar;

% 2-down 1-up staircase method on m
if ~iscorrect
    expvar = expvar + stepsize;
elseif n_correctinarow==2
    expvar = expvar - stepsize;
    n_correctinarow=0;
end
if expvar>1
    expvar=1;
end
if n_correctinarow == 0 % change stepsize
    if stepsize>cfg.min_stepsize
        stepsize = stepsize*cfg.adapt_stepsize;
    else
        stepsize = cfg.min_stepsize;
    end
end
str_inout.n_correctinarow = n_correctinarow;
str_inout.stepsize        = stepsize;
str_inout.expvar          = expvar;
