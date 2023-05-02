function ax = nexttile(N)
% function ax = nexttile(N)
%
% 1. Description:
%      This function allows forward compatibility of the function nexttile
%      if the MATLAB version being used is older than R2019b.
%
% Author: Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bExist = exist('nexttile.p','file'); % nexttile
if bExist
    if nargin == 0
        nexttile;
    else
        nexttile(N);
    end
else
    if nargin ~=0
        if N == 1
            close;
        end
    end
    figure;
end

if nargout ~= 0
    ax = gcf;
end