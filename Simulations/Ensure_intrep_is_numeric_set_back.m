function R = Ensure_intrep_is_numeric_set_back(intrep,info)
% function R = Ensure_intrep_is_numeric_set_back(intrep,info)
%
% 1. Description:
%
% 2. Stand-alone example:
%       [RMT,xx,info] = Ensure_intrep_is_numeric(RMT);
%       R = Ensure_intrep_is_numeric_set_back(intrep,info);
%       % RMT and R should be exactly the same
% 
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2017
% Created on    : 22/09/2017
% Last update on: 22/09/2017 
% Last use on   : 22/09/2017 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isfield(info,'size_cell')
    N = info.size_each_col;
    M = info.size_each_row;
    R = reshape(intrep,N,M);
else
    for i = 1:info.size_cell(1)
        N = info.size_each_cell(i,1); 
        M = info.size_each_cell(i,2);
        R{i,1} = [];
        for j = 1:M
            R{i,1} = [R{i,1} intrep(1:N)];
            intrep(1:N) = [];
        end
    end
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
