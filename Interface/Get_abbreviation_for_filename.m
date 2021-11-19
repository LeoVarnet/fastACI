function [str_out,str_in] = Get_abbreviation_for_filename(str_in)
% function [str_out,str_in] = Get_abbreviation_for_filename(str_in)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch str_in
    case 'up-to-trial'
        str_out = 't';
    otherwise
        str_out = str_in;
end     
