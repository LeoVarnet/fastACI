function str_inout = staircase_init(str_inout,cfg)
% function str_inout = staircase_init(str_inout,cfg)
%
% Initialise a new staircase
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iswarmup = cfg.warmup;
debut_i  = str_inout.debut_i;

response=[];
n_correctinarow = 0;
if iswarmup
    m = 0;
else
    m = cfg.m_start;
end
i = debut_i;
stepsize = cfg.start_stepsize;
isbreak = 0;

str_inout.response = response;
str_inout.n_correctinarow = n_correctinarow;
str_inout.m = m;
str_inout.i = i;
str_inout.stepsize = stepsize;
str_inout.isbreak = isbreak;
