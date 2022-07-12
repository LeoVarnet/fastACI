function [h_image,h_colourbar] = my_image_plot(indata, opts)
% MY_PLOT_CORRMAT  % Calculate and/or Plot Correlation Matrix
%   
% Inputs:
%   indata: Input matrix. It represents the matrix that will be plotted.
%
%   opts:   Struct with optional parameters that will overwrite the defaults
%       in this function.
%
%   The following fields can be defined as input to the opts struct:
%   'title',            string denoting title (default: 'Correlations')
% 
%   'labels',           a cell array of ROI labels, such as 
%                       {'ROI_1','ROI_2',...etc.}
%   
%   'label_FontSize'
% 
%   'corrmat',          this input specifies an already calculated
%                       correlation matrix which is input solely for
%                       plotting purposes; in this case, "Plot" and
%                       "Preprocessing" menus will be disabled;
%                       'timeSeries' input should be specified as an empty
%                       matrix ('[]') in this case; Other valid inputs for
%                       corrmat mode include 'plot','title','labels',
%                       'colormap,'sort_ind', and 'print'
% 
%
%   The following fields are not yet dynamic and therefore defaults are used
%       for them:
%   'ROI',              -Character/String: specified as a filename for a 3D 
%                       image of ROIs, denoted by sequential integer voxel 
%                       values. -Cell Array: wherein each cell contains the
%                       subscript indices of an ROI (row: voxel #; column:
%                       [x,y,z] subscripts).
% 
%   'nuisance',         a cell array (one cell per scan) containing a 2D
%                       matrix (rows = time points; cols = nuisance 
%                       variables); this option can be useful, e.g., if one
%                       needs to regress out subject-specific motion
%                       parameters from the data, or if one wants to
%                       investigate the effect made on correlation matrices
%                       by different preprocessing options
% 
%   'detrend',          integer value specifying whether and what type of
%                       detrending is desired before computation of the
%                       power spectrum, wherein 0 = no detrending, 1 =
%                       linear detrend, 2 = quadratic detrend (default = 0)
%                       -note: quadratic detrending removes both linear and
%                       quadratic trends
% 
%   'bandpass',         vector specified as [HighPass,LowPass], wherein
%                       HighPass = the lower frequency limit (Hz), and
%                       LowPass = the upper frequency limit (Hz); bandpass
%                       filtering is implemented using Brainstorm's 
%                       bst_bandpass_filtfilt method (John Mosher, Francois
%                       Tadel, 2014)
% 
%   'Fs',               numeric: sampling frequency (Hz) = 1 / sample_time (s)
%                       -this only needs to be specified if bandpass
%                       filtering is requested
% 
%   'plot',             TRUE(1)/FALSE(0); (default: TRUE)
% 
%   'colormap',         char/string: 'jet','hot','pink',etc.; this can be
%                       edited interactively as well
% 
%   'which_scan',       integer specifying which scan number to display
%                       correlation matrix for, or 0 to display Grand Mean
%                       (default is 0 if multiple scans, 1 otherwise)
% 
%   'sort_ind',         numeric vector: indices for sorting rows & cols of 
%                       correlation matrices
% 
%   'print',            characer/string: 'path\filename'; if specified, the
%                       fully initialized figure will automatically print
%                       to the filename specified; possible extensions
%                       (.png, .tiff, .bmp) are automatically detected from
%                       from the filename; if not specified, defaults to
%                       .png (in 300 dpi)
% 
% 'insert_axes',        this specifies a pre-existing axes object in which
%                       to plot the corrmat, rather than generate a new
%                       axes and figure
% 
% 
% Preprocessing Types:
% 1 = raw; 2 = linear detrend; 3 = quadratic detrend; 4 = bandpass; 
% 5 = linear detrend & bandpass; 6 = quadratic detrend & bandpass; 
% 7 = nuisance only; 8 = nuisance, linear detrend; 
% 9 = nuisance, quadratic detrend; 10 = nuisance, linear detrend, bandpass;
% 11 = nuisance, quadratic detrend, bandpass; 12 = nuisance, bandpass
% 
% This script was originally coded by Elliot Layden (2016)
% I (Alejandro Osses) modified the script to simplify it and change the way
%    the optional structure is passed.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get name-value pair arguments:
if nargin < 2
    opts = [];
end

indata_orig = indata;
if isfield(opts,'cmin') && isfield(opts,'cmax')
    cmin = opts.cmin;
    cmax = opts.cmax;
else
    cmin = min( min(min(indata)), 0);
    cmax = max(max(indata));
end

if ~isfield(opts,'bText')
    opts.bText = 1;
end
bText = opts.bText;

opts = Ensure_field(opts,'bColourbar',1);
bColourbar = opts.bColourbar;

alpha_data = ones(size(indata));

idx = find(isnan(indata));
indata(idx) = 0;

crange = cmax-cmin;
m = 257;
indata = min(m,round((m-1)*(indata-cmin)/(cmax-cmin)));
h_image = image(indata,'AlphaData',alpha_data);

set(gca,'XTick',1:size(indata,1));
set(gca,'YTick',1:size(indata,1));

if isfield(opts,'Labels')
    set(gca,'XTickLabel',opts.Labels);
    set(gca,'YTickLabel',opts.Labels);
end

if bText
    for i = 1:size(indata,1)
        for j = 1:size(indata,2)
            if ~isnan(indata_orig(i,j))
                text(j,i,sprintf('%.1f',indata_orig(i,j)),'HorizontalAlignment','center');
            end
        end
    end
end

if bColourbar
    cmap = il_bluewhitered(m,cmin,cmax);
    colormap(cmap);
    h_colourbar = colorbar; % ('southoutside');%,'Position',colorbar_pos); 
    Ticks      = get(h_colourbar,'Ticks');
    TickLabels = round(10*crange*(Ticks-1)/m)/10;
    set(h_colourbar,'TickLabels',TickLabels);
else
    h_colourbar = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function newmap = il_bluewhitered(m,cmin,cmax)
%BLUEWHITERED   Blue, white, and red color map.
%   BLUEWHITERED(M) returns an M-by-3 matrix containing a blue to white
%   to red colormap, with white corresponding to the CAXIS value closest
%   to zero.  This colormap is most useful for images and surface plots
%   with positive and negative values.  BLUEWHITERED, by itself, is the
%   same length as the current colormap.
% 
% Adapted from:
% Author: Nathan Childress (2003) 
% https://www.mathworks.com/matlabcentral/fileexchange/4058-bluewhitered
if nargin < 1, m = size(get(gcf,'colormap'),1); end
bottom = [0 0 0.5]; botmiddle = [0 0.5 1]; middle = [1 1 1]; 
topmiddle = [1 0 0]; top = [0.5 0 0];
% Find middle
if nargin < 2
    lims = get(gca, 'CLim');
else
    lims = [cmin,cmax];
end
% Find ratio of negative to positive
if (lims(1) < 0) && (lims(2) > 0)
    % It has both negative and positive
    % Find ratio of negative to positive
    ratio = abs(lims(1)) / (abs(lims(1)) + lims(2));
    neglen = round(m*ratio);
    poslen = m - neglen;
    % Just negative
    new = [bottom; botmiddle; middle];
    len = length(new);
    oldsteps = linspace(0, 1, len);
    newsteps = linspace(0, 1, neglen);
    newmap1 = zeros(neglen, 3);
    for i=1:3
        % Interpolate over RGB spaces of colormap
        newmap1(:,i) = min(max(interp1(oldsteps,new(:,i),newsteps),0), 1);
    end
    % Just positive
    new = [middle; topmiddle; top];
    len = length(new);
    oldsteps = linspace(0, 1, len);
    newsteps = linspace(0, 1, poslen);
    newmap = zeros(poslen, 3);
    for i=1:3
        % Interpolate over RGB spaces of colormap
        newmap(:,i) = min(max(interp1(oldsteps,new(:,i),newsteps), 0), 1);
    end
    % combine
    newmap = [newmap1; newmap];
elseif lims(1) >= 0
    % Just positive
    new = [middle; topmiddle; top];
    len = length(new);
    oldsteps = linspace(0, 1, len);
    newsteps = linspace(0, 1, m);
    newmap = zeros(m, 3);
    for i=1:3
        % Interpolate over RGB spaces of colormap
        newmap(:,i) = min(max(interp1(oldsteps,new(:,i),newsteps), 0), 1);
    end   
else
    % Just negative
    new = [bottom; botmiddle; middle];
    len = length(new);
    oldsteps = linspace(0, 1, len);
    newsteps = linspace(0, 1, m);
    newmap = zeros(m, 3);
    for i=1:3
        % Interpolate over RGB spaces of colormap
        newmap(:,i) = min(max(interp1(oldsteps,new(:,i),newsteps),0), 1);
    end  
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function c = il_redblue(m)
% 
% if nargin < 1, m = size(get(gcf,'colormap'),1); end
% 
% if (mod(m,2) == 0) % if even
%     % From [0 0 1] to [1 1 1], then [1 1 1] to [1 0 0];
%     m1 = .5*m;
%     r = (0:m1-1)'/max(m1-1,1);
%     g = r;
%     r = [r; ones(m1,1)];
%     g = [g; flipud(g)];
%     b = flipud(r);
% else % if odd
%     % From [0 0 1] to [1 1 1] to [1 0 0];
%     m1 = floor(m*0.5);
%     r = (0:m1-1)'/max(m1,1);
%     g = r;
%     r = [r; ones(m1+1,1)];
%     g = [g; 1; flipud(g)];
%     b = flipud(r);
% end
% 
% c = [r g b];
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Nuisance Regression: perform_nuisance_regression(parsed_inputs.nuisance)
% % function nuisance_regression
% %     data_types(7) = 1;
% %     for ixxx = 1:nSeries
% %         X = parsed_inputs.nuisance{ixxx};
% %         % Check for constant term
% %         if ~any(all(X==1)); X = [ones(size(X,1),1),X]; end %#ok
% %         for jxxx = 1:numrois
% %             [~,~,alldata{7}{ixxx}(:,jxxx)] = regress(alldata{1}{ixxx}(:,jxxx),X);
% %         end
% %     end
% % end
% 
% % Detrending:  perform_detrend(parsed_inputs.detrend)
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function il_perform_detrend(type)
%     switch type
%         case 1 % parsed_inputs.detrend==1
%             if ~data_types(2)
%                 disp('Performing linear detrending.')
%                 data_types(2) = 1;
%                 for ixxx = 1:nSeries
%                     alldata{2}{ixxx} = detrend(alldata{1}{ixxx});
%                 end
%             end
%             if ~data_types(8)
%                 disp('Performing linear detrending on denoised data.')
%                 data_types(8) = 1;
%                 for ixxx = 1:nSeries
%                     alldata{8}{ixxx} = detrend(alldata{7}{ixxx});
%                 end
%             end
%         case 2 % parsed_inputs.detrend==2
%             if ~data_types(3)
%                 disp('Performing quadratic detrending.')
%                 x = (1:nTime)';
%                 for ixxx = 1:nSeries
%                     for jxxx = 1:numrois
%                         p = polyfit(x,alldata{1}{ixxx}(:,jxxx),2);
%                         predicted = polyval(p,x);
%                         alldata{3}{ixxx}(:,jxxx) = alldata{1}{ixxx}(:,jxxx)-predicted;
%                     end
%                 end
%                 data_types(3) = 1;
%             end
%             if ~data_types(9)
%                 disp('Performing quadratic detrending on denoised data.')
%                 x = (1:nTime)';
%                 for ixxx = 1:nSeries
%                     for jxxx = 1:numrois
%                         p = polyfit(x,alldata{7}{ixxx}(:,jxxx),2);
%                         predicted = polyval(p,x);
%                         alldata{9}{ixxx}(:,jxxx) = alldata{7}{ixxx}(:,jxxx)-predicted;
%                     end
%                 end
%                 data_types(9) = 1;
%             end
%     end
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Bandpass
% function il_perform_bandpass(initial)
%     if all(parsed_inputs.bandpass==0)
%         % Prompt User:
%         prompt = {'Min Frequency (Hz):','Max Frequency (Hz):'};
%         dlg_title = 'Band-Pass Settings'; num_lines = [1,40;1,40]; defaultans = {'.008','.1'};
%         answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
%         if isempty(answer); disp('User cancelled action.'); return; end
%         if ~any(isnan(str2double(answer)))
%             parsed_inputs.bandpass(1) = str2double(answer{1});
%             parsed_inputs.bandpass(2) = str2double(answer{2});
%         end
%     end
%     if parsed_inputs.bandpass(1) <= (.1*(fs/4))
%         parsed_inputs.bandpass(1) = (.1*(fs/4))+.0001;
%         warning('HighPass set too low for sampling frequency (Fs), adjusted to %.4g Hz',parsed_inputs.bandpass(1))
%     end
%     if parsed_inputs.bandpass(2) >= (fs/2)
%         parsed_inputs.bandpass(2) = (fs/2)-.001;
%         warning('LowPass set too high for sampling frequency (Fs), adjusted to %.4g Hz',parsed_inputs.bandpass(2))
%     end
%     for iter = find((data_types(1:3)-data_types(4:6))==1)
%         for ixxx = 1:nSeries
%             try
%             [bandpassed,~,~] = il_bst_bandpass_filtfilt(alldata{iter}{ixxx}',fs,...
%                 parsed_inputs.bandpass(1), parsed_inputs.bandpass(2), 0, 'iir');
%             catch
%                 error('Bandpass Filter Error: Check if ''Fs'' was specified accurately.')
%             end
%             alldata{iter+3}{ixxx} = bandpassed';
%         end
%         if ~initial
%             calc_corrs(iter+3);
%         end
%         data_types(iter+3) = 1;
%     end
%     % Now repeat for nuisance regression:
%     if nuisance_on
%         iter2 = [12, 10, 11];
%         ind_nuis = find((data_types([7, 8, 9])-data_types(iter2))==1);
%         for iter = ind_nuis
%             out_ind = iter2(iter);
%             for ixxx = 1:nSeries
%                 try
%                 [bandpassed,~,~] = il_bst_bandpass_filtfilt(alldata{iter+6}{ixxx}',fs,...
%                     parsed_inputs.bandpass(1), parsed_inputs.bandpass(2), 0, 'iir');
%                 catch
%                     error('Bandpass Filter Error: Check if ''Fs'' was specified accurately.')
%                 end
%                 alldata{out_ind}{ixxx} = bandpassed';
%             end
%             if ~initial
%                 calc_corrs(out_ind);
%             end
%             data_types(out_ind) = 1;
%         end
%     end
