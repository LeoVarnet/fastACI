function fig_name = name2figname(name)
% function fig_name = name2figname(name)
%
% Replaces the 'non-compatible' characters for displaying names of files as
% it while plotting...
%
% 1. Description:
%   Replaces each '_' character '-'. 
%
% 2. Additional info:
%   Tested cross-platform: Yes
%
% 3. Stand-alone example:
%   name = 'hola_a_todos';
%   fig_name = name2figname(name);
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 2012-2013
% Last update on: 15/05/2014 % Update this date manually
% Last use on   : 15/05/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fig_name        = name;
cont            = find(name == '_');
fig_name(cont)  = '-';

if length(fig_name) > 24
    try
        [tmp1 fig_name] = fileparts(fig_name);
    catch
        fig_name = ['[..]' fig_name(end-24:end)];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
