function [intrep, size_intrep, info] = Ensure_intrep_is_numeric(intrep)
% function [intrep, size_intrep, info] = Ensure_intrep_is_numeric(intrep)
%
% 1. Description:
%       size_intrep is the size of the internal representation after being 
%       converted to a numeric value (matrix form).
% 
% 2. Stand-alone example:
%       [RMT,xx,info] = Ensure_intrep_is_numeric(RMT);
%       R = Ensure_intrep_is_numeric_set_back(intrep,info);
%       % RMT and R should be exactly the same
% 
% 3. Additional info:
%       Tested cross-platform: No
%       See also Ensure_intrep_is_numeric_set_back.m, il_ensure_is_numeric in exp_osses2016_3AFC.m
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2017
% Created on    : 14/05/2017
% Last update on: 29/08/2017 
% Last use on   : 29/08/2017 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~iscell(intrep)
    size_intrep = size(intrep);
    if nargout >=3
        info.size_each_col = size(intrep,1);
        info.size_each_row = size(intrep,2);
    end
    intrep = intrep(:);
else
    if nargout >= 3
        count_i = 0;
        info.is_cell = 1;
        info.size_cell = size(intrep);
        for i = 1:length(intrep)
            info.size_each_cell(i,:) = size(intrep{i});
            info.band_i(i) = count_i+1;
            count_i = count_i + info.size_each_cell(i,1)*info.size_each_cell(i,2);
        end
    end
    intrep = cell2mat(transpose(intrep)); 
    size_intrep = size(intrep);
    intrep = intrep(:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

