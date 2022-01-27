function o = optimget_local(options,name,default,flag)
%OPTIMGET_LOCAL Get OPTIM OPTIONS parameters.
%   VAL = OPTIMGET_LOCAL(OPTIONS,'NAME') extracts the value of the named parameter
%   from optimization options structure OPTIONS, returning an empty matrix if
%   the parameter value is not specified in OPTIONS.  It is sufficient to
%   type only the leading characters that uniquely identify the
%   parameter.  Case is ignored for parameter names.  [] is a valid OPTIONS
%   argument.
%
%   VAL = OPTIMGET_LOCAL(OPTIONS,'NAME',DEFAULT) extracts the named parameter as
%   above, but returns DEFAULT if the named parameter is not specified (is [])
%   in OPTIONS.  For example
%
%     val = optimget(opts,'TolX',1e-4);
%
%   returns val = 1e-4 if the TolX parameter is not specified in opts.
%
%   See also OPTIMSET.

%   Copyright 1984-2011 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2012/02/06 17:27:10 $

% undocumented usage for fast access with no error checking
if (nargin == 4) && isequal(flag,'fast')
    o = optimgetfast(options,name,default);
    return
end

if nargin < 2
    error(message('MATLAB:optimget:NotEnoughInputs'));
end
if nargin < 3
    default = [];
end

if ~isempty(options) && ~isa(options,'struct')
    error(message('MATLAB:optimget:Arg1NotStruct'));
end

if isempty(options)
    o = default;
    return;
end

allfields = {'Display'; 'MaxFunEvals';'MaxIter';'TolFun';'TolX';'FunValCheck';'OutputFcn';'PlotFcns'};

% Include specialized options if appropriate
if uselargeoptimstruct
    optimfields = optimoptiongetfields;  
    allfields = [allfields; optimfields];
end

Names = allfields;

name = deblank(name(:)'); % force this to be a row vector
j = find(strncmpi(name,Names,length(name)));
if isempty(j)               % if no matches
    error(message('MATLAB:optimget:InvalidPropName', name));
elseif length(j) > 1            % if more than one match
    % Check for any exact matches (in case any names are subsets of others)
    k = find(strcmpi(name,Names));
    if length(k) == 1
        j = k;
    else
        msg = ['(' Names{j(1),:}];
        for k = j(2:length(j))'
            msg = [msg ', ' Names{k,:}];
        end
        msg = [msg, '.)'];
        error(message('MATLAB:optimget:AmbiguousPropName', name, msg));
    end
end

if any(strcmp(Names,Names{j,:}))
    o = options.(Names{j,:});
    if isempty(o)
        o = default;
    end
else
    o = default;
end

%------------------------------------------------------------------
function value = optimgetfast(options,name,defaultopt)
%OPTIMGETFAST Get OPTIM OPTIONS parameter with no error checking so fast.
%   VAL = OPTIMGETFAST(OPTIONS,FIELDNAME,DEFAULTOPTIONS) will get the
%   value of the FIELDNAME from OPTIONS with no error checking or
%   fieldname completion. If the value is [], it gets the value of the
%   FIELDNAME from DEFAULTOPTIONS, another OPTIONS structure which is
%   probably a subset of the options in OPTIONS.
%

if isempty(options)
     value = defaultopt.(name);
     return;
end
% We need to know if name is a valid field of options, but it is faster to use 
% a try-catch than to test if the field exists and if the field name is
% correct. If the options structure is from an older version of the
% toolbox, it could be missing a newer field.
try
    value = options.(name);
catch ME
    value = [];
end

if isempty(value)
    value = defaultopt.(name);
end


