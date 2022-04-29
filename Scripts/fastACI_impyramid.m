function B = fastACI_impyramid(A, direction,keyvals)
%FASTACI_IMPYRAMID 
% My own modified version of MATLAB function impyramid. There are two small
% differences:
%   1. The variable outputSize = 2*[M N] (instead of 2*[M N] - 1) when direction == 'expand'
%   2. The input variable A is assumed to have dimensions: samples x frequency x time.
%      A is then temporarily set to frequency x time x samples, and after the
%      impyramid processing has been performed is set back to samples x frequency x time
%
% Author: Leo Varnet
% Old name: Script4_Calcul_ACI_modified_impyramid.m (before 28/04/2022)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A = permute(A,[2 3 1]);

validateattributes(A, {'numeric', 'logical'}, {}, mfilename, 'A', 1);
direction = validatestring(direction, {'reduce', 'expand'}, ...
    mfilename, 'DIRECTION', 2);
        
switch keyvals.pyramid_script
    case 'imgaussfilt'
        % Updated method used in the fastACI toolbox as default as of 28/04/2022
        
        i_level         = keyvals.i_level;
        pyramid_shape   = keyvals.pyramid_shape; % 0 for a 'standard' pyramid decomposition, -1 to start subsampling at level 3 only
        pyramid_padding = keyvals.pyramid_padding;
        sigma = i_level-1;                        
        switch direction
            case 'reduce'
                Ablur = imgaussfilt(A,sigma,'Padding',pyramid_padding); % blur
                if pyramid_shape == -1
                    warning('Need further validation (AO on 29/04/2022)');
                    idx_step = 2^(i_level-1+pyramid_shape);
                else
                    idx_step = 2; % comment added by AO
                end
                B = Ablur(1:idx_step:end,1:idx_step:end,:); % downsample, old name for B was Ablursub
            case 'expand'
                
                N_iterations = size(A,3);
                if isfield(keyvals,'Pyra_size')
                    Nf_and_Nt = keyvals.Pyra_size(2,:);
                else
                    Nf_and_Nt = [size(A,1) size(A,2)];
                end
                FinalSize = [(2+pyramid_shape)*Nf_and_Nt N_iterations]; % Temporary solution
                % zero-padding
                Azeroed = zeros(FinalSize);
                idx_step = size(Azeroed,1)/size(A,1); % idx_step = 2^(i_level-1+pyramid_shape);
                Azeroed(1:idx_step:end,1:idx_step:end,:) = A;
                % blur
                UR = imgaussfilt(1,sigma,'Padding',keyvals.pyramid_padding);
                Ablur = (1/UR)*imgaussfilt(Azeroed,sigma,'Padding',keyvals.pyramid_padding);
                
                B = Ablur;
        end
                    
    case 'imresize'
        % First method used in the fastACI toolbox (default until 27/04/2022):
        
        M = size(A,1);
        N = size(A,2);

        switch direction
            case 'reduce'
                scaleFactor = 0.5;
                outputSize = ceil([M N]*scaleFactor);
                kernel = il_makePiecewiseConstantFunction( ...
                    [3.5 2.5      1.5    0.5    -0.5   -1.5    -Inf], ...
                    [0.0 0.0625   0.25   0.375   0.25   0.0625  0.0]);
                kernelWidth = 5;

            case 'expand'
                scaleFactor = 2;
                outputSize = scaleFactor*[M N];
                kernel = il_makePiecewiseConstantFunction(...
                    [1.25   0.75    0.25   -0.25   -0.75   -1.25   -Inf],...
                    [0.0    0.125   0.5     0.75    0.5    0.125    0.0]);
                kernelWidth = 3;
        end

        B = imresize(A, scaleFactor, {kernel, kernelWidth}, ...
            'OutputSize', outputSize, 'Antialiasing', false);
end

B = permute(B,[3 1 2]);

end
% End of file

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions used in the script:
function fun = il_makePiecewiseConstantFunction(breakPoints, values)
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