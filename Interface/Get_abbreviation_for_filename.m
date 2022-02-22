function [str_out,str_in] = Get_abbreviation_for_filename(str_in)
% function [str_out,str_in] = Get_abbreviation_for_filename(str_in)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch str_in % sorted alphabetically
    case 'addsignal'
        str_out = 'sig';
    case 'bias'
        str_out = 'b';
    case 'correct'
        str_out = 'co';
    case 'gammatone'
        str_out = 'gt';
    case 'incorrect'
        str_out = 'inco';
    case 'lasso'
        str_out = 'lasso'; % old (before 22/02/2022)='l1lm'
    case {'l1lm','lassoslow'}
        str_out = 'l1lm'; % old (before 22/02/2022)='l1lms'
    case 'lassoglm' 
        str_out = 'lassoglm'; % old (before 22/02/2022)='l1glm'
    case {'l1glm','lassoglmslow'}
        str_out = 'l1glm'; % old (before 22/02/2022)='l1glms'
    case 'no_bias'
        str_out = 'nob';
    case 'perc'
        str_out = 'pe';
    case 'up-to-trial'
        str_out = 'tr';
    otherwise
        str_out = str_in;
end     
