function B = Script4_Calcul_ACI_modified_impyramid(A, direction)
%MODIFIED_IMPYRAMID 
% My own modified version of MATLAB function impyramid. THe only difference
% is that outputSize = 2*[M N] instead of 2*[M N] - 1
%
% Author: Leo Varnet
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A = permute(A,[2 3 1]);

validateattributes(A, {'numeric', 'logical'}, {}, mfilename, 'A', 1);
direction = validatestring(direction, {'reduce', 'expand'}, ...
    mfilename, 'DIRECTION', 2);

M = size(A,1);
N = size(A,2);

if strcmp(direction, 'reduce')
    scaleFactor = 0.5;
    outputSize = ceil([M N]/scaleFactor);
    kernel = makePiecewiseConstantFunction( ...
        [3.5 2.5      1.5    0.5    -0.5   -1.5    -Inf], ...
        [0.0 0.0625   0.25   0.375   0.25   0.0625  0.0]);
    kernelWidth = 5;
    
else
    scaleFactor = 2;
    outputSize = scaleFactor*[M N];
    kernel = makePiecewiseConstantFunction(...
        [1.25   0.75    0.25   -0.25   -0.75   -1.25   -Inf],...
        [0.0    0.125   0.5     0.75    0.5    0.125    0.0]);
    kernelWidth = 3;
end

B = imresize(A, scaleFactor, {kernel, kernelWidth}, ...
    'OutputSize', outputSize, 'Antialiasing', false);
B = permute(B,[3 1 2]);

end

function fun = makePiecewiseConstantFunction(breakPoints, values)
% Constructs a piecewise constant function and returns a handle to it.
%
% breakPoints and values have to be vectors with the same number of
% elements.
%
% The elements in breakPoints have to be monotonically decreasing.
% 
% fun(x) returns values(1) if x >= breakPoints(1)
%
% else fun(x) returns values(2) if x >= breakPoints(2)
%
% else fun(x) returns values(3) if x >= breakPoints(3)
%
% etc.
%
% If x is less than breakPoint(end), then fun returns 0.
%
% If x is an array, then fun operates elementwise on x and returns an array
% of the same size.

%iptassert(all(diff(breakPoints) < 0), ...
%          'images:modified_impyramid:badBreakPointList')

fun = @piecewiseConstantFunction;

    function y = piecewiseConstantFunction(x)
        y = zeros(size(x));
        for k = 1:numel(x)
            yy = 0;
            xx = x(k);
            for p = 1:numel(breakPoints)
                if xx >= breakPoints(p)
                    yy = values(p);
                    break;
                end
            end
            y(k) = yy;
        end
    end
end