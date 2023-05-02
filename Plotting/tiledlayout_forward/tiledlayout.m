function handle_tl = tiledlayout(N,M,TileSpacing,TileSpacing_option,varargin)
% function handle_tl = tiledlayout(N,M,TileSpacing,TileSpacing_option,varargin)
%
% 1. Description:
%      This function allows forward compatibility of the function tiledlayout
%      if the MATLAB version being used is older than R2019b.
%
% Author: Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

handle_tl = [];
if nargin < 4
    TileSpacing_option = 'Compact';
end
bExist = exist('tiledlayout.p','file'); % tiledlayout.p
if bExist
    if isempty(varargin)
        handle_tl = tiledlayout(N,M,TileSpacing,TileSpacing_option);
    else
        handle_tl = tiledlayout(N,M,TileSpacing,TileSpacing_option,varargin{:});
    end
else
    % for i = 1:N*M
    %     handle_tl(i) = subplot(N,M,i);
    % end
    warning('We programmed this script to use more recent graphic options from MATLAB. It might be that you won''t be able to visualise these results as we have foreseen...')
end