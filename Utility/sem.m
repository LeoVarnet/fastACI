function y = sem(varargin)
%SEM Standard error of the mean.
%   For vectors, Y = SEM(X) returns the standard error of the mean.  For matrices,
%   Y is a row vector containing the standard error of the mean of each column.  For
%   N-D arrays, SEM operates along the first non-singleton dimension of X.

y = sqrt(var(varargin{:}))/sqrt(length(varargin{1}));
