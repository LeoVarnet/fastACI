function [var_out, bFound] = Get_field_from_cell(cell_in,str_in)
% function [var_out, bFound] = Get_field_from_cell(cell_in,str_in)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

var_out = [];

bFound = 0;

for i = 1:length(cell_in)-1
    switch cell_in{i}
        case str_in
            var_out = cell_in{i+1};
            bFound = 1;
            break;
    end
end
