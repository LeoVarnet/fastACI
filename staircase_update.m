function str_inout = staircase_update(str_inout,cfg)
% function str_inout = staircase_update(str_inout,cfg)
%
% staircase update phase
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iscorrect       = str_inout.iscorrect;
n_correctinarow = str_inout.n_correctinarow;
stepsize        = str_inout.stepsize;
m               = str_inout.m;

% 2-down 1-up staircase method on m
if ~iscorrect
    m = m + stepsize;
elseif n_correctinarow==2
    m = m - stepsize;
    n_correctinarow=0;
end
if m>1
    m=1;
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
str_inout.m               = m;
