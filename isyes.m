function [ boolean ] = isyes( A )
%ISYES(A) returns true if A is a string containing 'oui' or 'yes' or 'on',
%and false if A is a string containing 'non' or 'no' or 'off'

if ischar(A)
    if strcmp(A, 'yes') || strcmp(A, 'oui') || strcmp(A, 'on');
        boolean = true;
    elseif strcmp(A, 'no') || strcmp(A, 'non') || strcmp(A, 'off');
        boolean = false;
    else
        error('Input argument of isyes is not a valid string (valid strings are : ''yes'', ''oui'', ''on'', ''no'', ''non'', ''off'')')
    end
else
    error('Input argument of isyes must be a string')
end
end

