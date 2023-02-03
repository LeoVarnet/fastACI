function [p,bStatus] = Ensure_field(p, field_name, default_value, bSilent)
% function [p,bStatus] = Ensure_field(p, field_name, default_value, bSilent)
%
% Ensure that a struct field exists, else give it a default value.
% If the field existed in the input p, then the output p is identical.
% Else a new field is created, with the specified default value.
% function p = Ensure_field(p, field_name, default_value);
% 
% Inputs:
%   p:             Parameter struct.
%   field_name:    Name of field (string).
%   default_value: Value to set field to if field does not exist in p.
%
% Outputs:
%   p:             Parameter struct.
%   bStatus:       1 if field was assigned, 0 if not
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Authors       : Brett Swanson
% Edited by     : Alejandro Osses (ale.a.osses@gmail.com), TU/e Eindhoven
% Last edited on: 24/03/2017
% Last used on  : 24/03/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 4
    bSilent = 0;
end

if ~isfield(p, field_name)
	p = setfield(p, field_name, default_value);
    bStatus = 1;
    if nargout == 1 & bSilent == 0
        disp([mfilename '.m: Struct field ''' field_name ''' assigned'])
    end
else
    bStatus = 0;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
