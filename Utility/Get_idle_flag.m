function flags_here = Get_idle_flag
% function flags_here = Get_idle_flag
%
% Important if MATLAB is R2012b (or maybe also older versions)

flags_here.idle = {'idle'}; 
flags_here.do_idle = 1; 
flags_here.do_no_idle = 0; % for compatibility with older MATLAB versions