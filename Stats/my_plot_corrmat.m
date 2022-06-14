function [h_corrmat, h_colorbar, corrmat] = my_plot_corrmat(indata, opts)
% MY_PLOT_CORRMAT  % Calculate and/or Plot Correlation Matrix
%   
% Inputs:
%   indata: Input matrix to be correlated. Each column is assumed to be 
%       related to a different treatment.
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

parsed_inputs = struct('ROI',[],'nuisance',[],'detrend',0,'bandpass',[0,0],...
    'Fs',[],'plot',1,'title','','labels',[],'colormap','jet','which_scan',0,...
    'sort_ind',[],'print',0,'corrmat',[],'insert_axes',[],'outline',0,...
    'label_FontSize',12,'colorbar_FontSize',12); % defaults

if isfield(opts,'corrmat')
    parsed_inputs.corrmat = opts.corrmat;
end

if isfield(opts,'label_FontSize')
    parsed_inputs.label_FontSize = opts.label_FontSize;
else; fprintf('The default label_FontSize was used\n'); end
   
if isfield(opts,'outline')
    parsed_inputs.outline = opts.outline;
else; fprintf('The default outline was used\n'); end

if isfield(opts,'labels')
    parsed_inputs.labels = opts.labels;
else; fprintf('The default labels were used\n'); end

if isfield(opts,'title')
    parsed_inputs.title = opts.title;
else; fprintf('The default labels are being used\n'); end

non_corrmat_mode = false;
numrois = size(parsed_inputs.corrmat,1);
if ~isempty(indata)
    
    if ~isfield(opts,'type')
        types = {'Spearman','Pearson'};
        Show_cell(types);
        bInput = input('Choose the type you want to use (numeric input): ');
        type = types{bInput};
    else
        type = opts.type;
    end
    corrmat = [];
    for i = 1:size(indata,2)
        ref_ma = indata(:,i);
        for j = 1:size(indata,2)
            corr_val   = corr(ref_ma,indata(:,j),'type',type);
            corrmat(i,j) = corr_val;
            corrmat(j,i) = corr_val;
        end
    end  
    numrois = size(corrmat,1);
    corrmat = {corrmat};
    
else
    corrmat = {parsed_inputs.corrmat};
end
which_scan = 1;
% multi_series = false;
curr_data_type = 1;

%% Initialise:
if non_corrmat_mode
    data_types = zeros(1,12); data_types(1) = 1;
    data_types_str = {'Raw: ','Linear Detrend: ','Quadratic Detrend: ',...
        'Bandpass Filter: ','Linear/Bandpass: ','Quadratic/Bandpass: ',...
        'Nuisance Regression: ','Nuisance/Linear: ','Nuisance/Quadratic: ',...
        'Nuisance/Linear/Bandpass: ','Nuisance/Quadratic/Bandpass:',...
        'Nuisance/Bandpass: '};
     
    % Initialise final time series cell array:
    alldata = cell(1,12);
end

N_cmap = '64';
h_image = 28.1873;
h_colorbar = 71.3747;
h_line = zeros(1,size(corrmat{1},1)); % outline handles
if isfield(opts,'cmin') && isfield(opts,'cmax')
    caxis_auto = false;
    cmin = opts.cmin;
    cmax = opts.cmax;
    error('Alejandro has to continue debugging...')
else
    caxis_auto = true;
    cmin = -.1; cmax = .1;
end
labels_on = true;
colorbar_on = true;
FontName = 'Arial'; FontSize = parsed_inputs.label_FontSize; FontWeight = 'normal';
title_FontName = 'Helvetica'; title_FontSize = 11; title_FontWeight = 'bold';
rb3_cmap = false;
main_colormap_selection = 'blue-red';
info_popup = 18.37272;
highlight_on = false;
highlight_color = [1,1,0]; % initialize to yellow

% ROI Indices:
sort_ind = 1:numrois;
% warning('''sort_ind'' should be of length equal to the # of ROIs')

alpha_data = triu(ones(numrois),1);
alpha_data2 = tril(ones(numrois),-1);

% Positioning:
if isfield(opts,'Position')
    figure_pos = [50 50 opts.Position];
    bCustom_pos = 1;
else
    screen_res = get(0,'MonitorPositions'); % get(0,'ScreenSize');
    figure_pos = [.236*screen_res(3),.063*screen_res(4),...
        .537*screen_res(3), .86*screen_res(4)];
    bCustom_pos = 0;
end

ax_pos = [.14,.1,.78,.88]; % leave x at .45 to leave room when large decimals on y-axis
colorbar_pos = [ax_pos(1),.034,ax_pos(3),.06];
ax_nocolorbar_pos = [.11,.01,.84,.97];
    
% Labels:
if isempty(parsed_inputs.labels) % default labels:
    parsed_inputs.labels = cell(1,numrois);
    for i = 1:numrois
        parsed_inputs.labels{i} = sprintf('ROI %02g',i);
    end
end

%% Calculate CorrMat

% Calculate Mean CorrMat if Multiple Series:
% 1 = raw, 2 = linear detrend, 3 = quadratic detrend, 4 = bandpass, 
% 5 = linear detrend & bandpass, 6 = quadratic detrend & bandpass
function calc_corrs(type) 
    corrmat{type} = zeros(numrois,numrois,nSeries+1);
    for iter2 = 1:nSeries
        r = corr(alldata{type}{iter2});
        r(logical(eye(numrois))) = nan; % null main diagonal
        corrmat{type}(:,:,iter2) = r;
    end 
    % Grand Mean:
    corrmat{type}(:,:,nSeries+1) = mean(corrmat{type},3);
end

%% Plot
if parsed_inputs.plot
    % Initialise Figure
    h_corrmat = figure('Position',figure_pos,'MenuBar','none',...
            'Name',parsed_inputs.title,'NumberTitle','on','Color',[1,1,1]); % [.8,.88,.98]
        
    file_menu = uimenu(h_corrmat,'Label','File');
    uimenu(file_menu,'Label','Save Figure','Callback',@il_save_figure_callback);
    uimenu(file_menu,'Label','Print','Callback',{@print_callback,0});
    
    % Plot Menus
    which_scan = 1;

        % Display Menu:
        display_menu = uimenu(h_corrmat,'Label','Display');
        colormap_menu = uimenu(display_menu,'Label','Colormap');
        uimenu(colormap_menu,'Label','blue-red'    ,'Callback',@colormap_callback);
        uimenu(colormap_menu,'Label','blue-red (2)','Callback',@colormap_callback);
        uimenu(colormap_menu,'Label','blue-red (3)','Callback',@colormap_callback);
        uimenu(colormap_menu,'Label','jet','Callback',@colormap_callback);
        uimenu(colormap_menu,'Label','hot','Callback',@colormap_callback);
        uimenu(colormap_menu,'Label','cool','Callback',@colormap_callback);
        uimenu(colormap_menu,'Label','gray','Callback',@colormap_callback);
        uimenu(colormap_menu,'Label','hsv','Callback',@colormap_callback);
        uimenu(colormap_menu,'Label','bone','Callback',@colormap_callback);
        uimenu(colormap_menu,'Label','copper','Callback',@colormap_callback);
        uimenu(colormap_menu,'Label','spring','Callback',@colormap_callback);
        uimenu(colormap_menu,'Label','summer','Callback',@colormap_callback);
        uimenu(colormap_menu,'Label','winter','Callback',@colormap_callback);
        uimenu(colormap_menu,'Label','pink','Callback',@colormap_callback);
        caxis_menu = uimenu(display_menu,'Label','Color Axis');
        uimenu(caxis_menu,'Label','Automatic','Callback',{@caxis_callback,1});
        uimenu(caxis_menu,'Label','Specify...','Callback',{@caxis_callback,0});
        title_menu = uimenu(display_menu,'Label','Title');
        uimenu(title_menu,'Label','Edit','Callback',@title_edit);
        uimenu(title_menu,'Label','FontName','Callback',{@title_font,1}); 
        uimenu(title_menu,'Label','FontSize','Callback',{@title_font,2}); 
        title_font_bold = uimenu(title_menu,'Label','FontWeight'); 
        uimenu(title_font_bold,'Label','Normal','Callback',{@title_font,3}); 
        uimenu(title_font_bold,'Label','Bold','Callback',{@title_font,4}); 
        labels_menu = uimenu(display_menu,'Label','Labels');
        uimenu(labels_menu,'Label','Enable/Disable','Callback',{@il_labels_callback,1});
        labels_font = uimenu(labels_menu,'Label','Font');
        uimenu(labels_font,'Label','FontName','Callback',{@il_labels_callback,2}); 
        uimenu(labels_font,'Label','FontSize','Callback',{@il_labels_callback,3}); 
        font_bold = uimenu(labels_font,'Label','FontWeight'); 
        uimenu(font_bold,'Label','Normal','Callback',{@il_labels_callback,4}); 
        uimenu(font_bold,'Label','Bold','Callback',{@il_labels_callback,5}); 
        colorbar_menu = uimenu(display_menu,'Label','ColorBar');
        uimenu(colorbar_menu,'Label','Enable/Disable','Callback',{@il_colorbar_callback,1});
        colorbar_font_menu = uimenu(colorbar_menu,'Label','Font');
        uimenu(colorbar_font_menu,'Label','FontName','Callback',{@il_colorbar_callback,2}); 
        uimenu(colorbar_font_menu,'Label','FontSize','Callback',{@il_colorbar_callback,3}); 
        colorbar_font_bold = uimenu(colorbar_font_menu,'Label','FontWeight'); 
        uimenu(colorbar_font_bold,'Label','Normal','Callback',{@il_colorbar_callback,4}); 
        uimenu(colorbar_font_bold,'Label','Bold','Callback',{@il_colorbar_callback,5}); 
        highlight_menu = uimenu(display_menu,'Label','Highlight');
        choose_color_menu = uimenu(highlight_menu,'Label','Choose Color');
        uimenu(choose_color_menu,'Label','Yellow','Callback',{@highlight_callback,1});
        uimenu(choose_color_menu,'Label','Green','Callback',{@highlight_callback,2});
        uimenu(choose_color_menu,'Label','Magenta','Callback',{@highlight_callback,3});
        uimenu(choose_color_menu,'Label','Cyan','Callback',{@highlight_callback,4});
        uimenu(choose_color_menu,'Label','Red','Callback',{@highlight_callback,5});
        uimenu(choose_color_menu,'Label','Blue','Callback',{@highlight_callback,6});
        uimenu(choose_color_menu,'Label','White','Callback',{@highlight_callback,7});
        uimenu(choose_color_menu,'Label','Black','Callback',{@highlight_callback,8});

        uimenu(highlight_menu,'Label','Disable','Callback',{@highlight_callback,0});
        outline_menu = uimenu(display_menu,'Label','Outline','Checked','off','Callback',@il_change_outline);
        
        % Initialise Axes
        ax = axes('Position',ax_pos,'XDir','reverse','YDir','reverse','XLim',[.5,numrois+.5],...
            'YLim',[.5,numrois-.5],'Box','off','XTick',[],'YTick',[],'Visible','off'); 
        hold(ax,'on'); view(-90,90);

        % Add Title:
        h_title = title(ax,parsed_inputs.title,'Visible','on','Units',...
            'normalized','PickableParts','all','FontName',title_FontName,...
            'FontSize',title_FontSize,'FontWeight',title_FontWeight); 
        centre_axes = (ax_pos(3)+ax_pos(1))-.5*ax_pos(3); % axes center in norm fig units
        adjust1 = (.5-centre_axes)/ax_pos(3); % adjustment needed in axes units
        title_pos = get(h_title,'Position');
        if bCustom_pos == 0
            title_pos(1) = .5 - .5*title_pos(3) + adjust1; 
            title_pos(2) = .99;
        else
            title_pos(2) = 0.95;
        end
        if isempty(parsed_inputs.insert_axes)
            set(h_title,'Position',title_pos);
        end
        % Update Plot:
        il_update_plot(1,curr_data_type);
end

% If 'outline' input is specified:
if parsed_inputs.outline
    hObject = struct('Checked','off');
    il_change_outline(hObject)
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function il_update_plot(initialise,curr_data_type,~)
    if ~initialise && isgraphics(h_image,'image')
        cla(ax); delete(h_colorbar);
    end
    % Determine CData Idx:
    imdata = corrmat{curr_data_type}(:,:,which_scan);
    if caxis_auto
        cmin = min(imdata(:)); cmax = max(imdata(:));  % min & max color value
    end
    if cmin==cmax; 
        cmax = cmin + 1; 
    end
    if nargin==3 %&& update_cmap % update cmap:
        colormap_callback([],[],[]);
    end
    if rb3_cmap
        m = 11;
    else
        m = str2double(N_cmap);
    end
    idx1 = min(m,round((m-1)*(imdata-cmin)/(cmax-cmin))+1);
    idx1(idx1<=0) = 1; % assure no negative or 0 indices
    main_colorbar_lim = [min(idx1(:)),max(idx1(:))];
    if main_colorbar_lim(1)>=main_colorbar_lim(2)
        main_colorbar_lim(2) = main_colorbar_lim(2)+1;
    end
    % Create Image:
    h_image = image('Parent',ax,'CData',idx1(sort_ind,sort_ind),'AlphaDataMapping','none',...
        'AlphaData',alpha_data,'ButtonDownFcn',@click_corrmat);
    if colorbar_on
        axes(ax); 
        h_colorbar = colorbar('southoutside','Position',colorbar_pos); 
        % Set Colorbar Ticks:
        h_colorbar.Limits = main_colorbar_lim;
        h_colorbar.LimitsMode = 'manual';
        h_colorbar.FontName = 'Arial';
        h_colorbar.FontSize = parsed_inputs.colorbar_FontSize;
        h_colorbar.FontWeight = 'bold';
        colorbar_main_cvec = linspace(cmin,cmax,m);
        if (cmax-cmin)>(.15*m) % if colorbar ticks should be integers
            colorbar_main_cvec = round(colorbar_main_cvec);
            h_colorbar.TickLabels = cellstr(sprintf('%1g\n',...
                colorbar_main_cvec(h_colorbar.Ticks)));                    
        else
            h_colorbar.TickLabels = cellstr(sprintf('%4.2g\n',...
                colorbar_main_cvec(h_colorbar.Ticks)));
        end
    end
    % Colormap:
    if initialise
        main_colormap = il_bluewhitered(str2double(N_cmap),cmin,cmax);
        colormap(ax,main_colormap);
    end
    % Add Labels:
    if labels_on
        for l = 2:numrois
            text(l,.42,parsed_inputs.labels{sort_ind(l)},'FontName',FontName,...
                'FontSize',FontSize,'FontWeight',FontWeight,...
                'HorizontalAlignment','right','Parent',ax);
        end
        for l = 1:numrois-1
            text(l+.2,l-.45,parsed_inputs.labels{sort_ind(l)},'FontName',FontName,... % .25, .4
                'FontSize',FontSize,'FontWeight',FontWeight,'Parent',ax);
        end
    end
    % Outline
    if strcmp(get(outline_menu,'Checked'),'on')
        outline_menu.Checked = 'off';
        il_change_outline(outline_menu)
    end
end

%% Callbacks:
rect_count = 0; h_rect = 12.3399;
function click_corrmat(~,event,~)
    click_ind = round(event.IntersectionPoint(1:2));
    if ~highlight_on && alpha_data2(click_ind(1),click_ind(2))
        if non_corrmat_mode
            % Get Data:
            stats = cell(1,sum(data_types)+1);
            stats{1} = [parsed_inputs.labels{sort_ind(click_ind(1))},...
                ' to ',parsed_inputs.labels{sort_ind(click_ind(2))}]; 
            count = 1;
            for ixxx = find(data_types)
                count = count+1;
                r1 = corrmat{ixxx}(sort_ind(click_ind(1)),sort_ind(click_ind(2)),which_scan);
                r1 = round(r1*100); r1 = r1*.01; % correct number of sig digits
                str = sprintf('%2g',r1);
                if r1>0
                    stats{count} = [data_types_str{ixxx},'r = ',str(2:end)];
                else
                    stats{count} = [data_types_str{ixxx},'r = -',str(3:end)];
                end
            end
        else
            stats = {};
            r1 = corrmat{1}(sort_ind(click_ind(1)),sort_ind(click_ind(2)));
            stats{1} = sprintf('%2.3g',r1);
        end
        if ishandle(info_popup); delete(info_popup); end % Only allow one at a time
        reversed = numrois:-1:1;
        info_pos = [click_ind(2)/numrois,reversed(click_ind(1))/numrois,.05,.05];
        if (click_ind(1)/numrois)>=.85; info_pos(2) = reversed(round(.85*numrois))/numrois; end
        if (click_ind(2)/numrois)>=.64; info_pos(1) = .64; end
        info_popup = annotation(h_corrmat,'textbox','Position',info_pos,...
            'String',stats,'BackgroundColor',[1,1,1],'EdgeColor',[0,0,0],...
            'FaceAlpha',1,'FitBoxToText','on','ButtonDownFcn',@remove_obj);
    elseif alpha_data2(click_ind(1),click_ind(2)) % highlight is on
        rect_count = rect_count+1;
        rect_pos = [click_ind(1)-.5,click_ind(2)-.5,1,1]; % [x,y,w,h]
        h_rect(rect_count) = rectangle('Parent',ax,'Position',rect_pos,...
            'Curvature',[.1,.1],'EdgeColor',highlight_color,'LineWidth',2,...
            'ButtonDownFcn',@remove_obj);
    end
end

function remove_obj(hObject,~,~)
    delete(hObject);
end

function highlight_callback(~,~,type)
    switch type
        case 0, highlight_on = false;
            for ixxx = 1:length(h_rect)
                if ishandle(h_rect(ixxx)); delete(h_rect(ixxx)); end
            end
            h_rect = 18.3727; rect_count = 0;
            for ixxx = 1:6
                if ishandle(h_oval(ixxx)); delete(h_oval(ixxx)); end
            end
            return;
        case 1, highlight_on = true;
            highlight_color = [1 1 0]; % yellow
        case 2, highlight_on = true;
            highlight_color = [0 1 0]; % green
        case 3, highlight_on = true;
            highlight_color = [1 0 1]; % magenta
        case 4, highlight_on = true;
            highlight_color = [0 1 1]; % cyan
        case 5, highlight_on = true;
            highlight_color = [1 0 0]; % red
        case 6, highlight_on = true;
            highlight_color = [0 0 1]; % blue
        case 7, highlight_on = true;
            highlight_color = [1 1 1]; % white
        case 8, highlight_on = true;      
            highlight_color = [0 0 0]; % black
    end
end

h_oval = rand(1,6);

function export_var(~,~,~)
    prompt = 'Specify variable name:';
    dlg_title = 'Export'; 
    answer = inputdlg(prompt,dlg_title);
    if ~isempty(answer) && ischar(answer{1})
        data2save = corrmat{curr_data_type}(sort_ind,sort_ind,which_scan);
        assignin('base',answer{1},data2save);
    end
end

function which_scan_callback(~,~,type)
    if which_scan~=type
        set(h_which_scan_menu(type),'Label',[get(h_which_scan_menu(type),'Label'),' *']);
        lab2 = get(h_which_scan_menu(which_scan),'Label');
        set(h_which_scan_menu(which_scan),'Label',lab2(1:end-2));
        which_scan = type;
        il_update_plot(0,curr_data_type);
    end
end

function nuisance_callback(~,~,~)
    if nuisance_curr
        nuisance_menu.Label = 'Nuisance Regression';
    else
        nuisance_menu.Label = 'Nuisance Regression *';
    end
    nuisance_curr = ~nuisance_curr;
    determine_curr_data_type;
    il_update_plot(0,curr_data_type);
end

function change_detrend_callback(~,~,type)
    if parsed_inputs.detrend~=type
        switch type
            case 0 
                detrend0_menu.Label = 'None *';
                detrend1_menu.Label = 'Linear';
                detrend2_menu.Label = 'Quadratic';   
            case 1
                detrend0_menu.Label = 'None';
                detrend1_menu.Label = 'Linear *';
                detrend2_menu.Label = 'Quadratic';
                if ~data_types(2) || ~data_types(8)
                    il_perform_detrend(1);
                    calc_corrs(2) ;
                    if nuisance_on; calc_corrs(8); end
                    if bandpass_on 
                        if ~data_types(5) || ~data_types(10)
                            il_perform_bandpass(0);
                        end
                    end
                end
            case 2
                detrend0_menu.Label = 'None';
                detrend1_menu.Label = 'Linear';
                detrend2_menu.Label = 'Quadratic *';
                if ~data_types(3) || ~data_types(9)
                    il_perform_detrend(2);
                    calc_corrs(3);
                    if nuisance_on; calc_corrs(9); end
                    if bandpass_on
                        if ~data_types(6) || ~data_types(11)
                            il_perform_bandpass(0);
                        end
                    end
                end
        end
        parsed_inputs.detrend = type;
        determine_curr_data_type;
        il_update_plot(0,curr_data_type);
    end
end

function change_bandpass_callback(~,~,type)
    if type~=2
        if bandpass_on~=type
            if type==1
                if any((data_types([1:3,7:9])-data_types([4:6,12,10,11]))==1)
                    il_perform_bandpass(0);
                end
                bandpass_on_menu.Label = 'On *';
                bandpass_off_menu.Label = 'Off';
            else
                bandpass_on_menu.Label = 'On';
                bandpass_off_menu.Label = 'Off *';
            end
            bandpass_on = type;
            determine_curr_data_type;
            il_update_plot(0,curr_data_type);
        end
    else
        parsed_inputs.bandpass = [0,0];
        data_types([4:6,10:12]) = 0;
        il_perform_bandpass(0);
        bandpass_on = 1;
        bandpass_on_menu.Label = 'On *';
        bandpass_off_menu.Label = 'Off';
        determine_curr_data_type;
        il_update_plot(0,curr_data_type);
    end
end

function determine_curr_data_type
    if ~nuisance_curr
        switch parsed_inputs.detrend
            case 0
                if bandpass_on
                    curr_data_type = 4;
                else
                    curr_data_type = 1;
                end
            case 1
                if bandpass_on
                    curr_data_type = 5;
                else
                    curr_data_type = 2;
                end
            case 2
                if bandpass_on
                    curr_data_type = 6;
                else
                    curr_data_type = 3;
                end  
        end
    else % nuisance regression is on
        switch parsed_inputs.detrend
            case 0
                if bandpass_on
                    curr_data_type = 12;
                else
                    curr_data_type = 7;
                end
            case 1
                if bandpass_on
                    curr_data_type = 10;
                else
                    curr_data_type = 8;
                end
            case 2
                if bandpass_on
                    curr_data_type = 11;
                else
                    curr_data_type = 9;
                end  
        end
    end
end

% Colormap Callback:
function colormap_callback(hObject,~,~)
    if ~isempty(hObject)
        main_colormap_selection = get(hObject,'Label');
    end
    if strcmp(main_colormap_selection,'blue-red') % default
        main_colormap = il_bluewhitered(str2double(N_cmap),cmin,cmax);
        rb3_cmap = false;
    elseif strcmp(main_colormap_selection,'blue-red (2)') % default
        main_colormap = il_redblue(str2double(N_cmap));
        rb3_cmap = false;
    elseif strcmp(main_colormap_selection,'blue-red (3)') % default
        main_colormap = redbluecmap(11);
        rb3_cmap = true;
    else
        rb3_cmap = false;
        main_colormap = eval([main_colormap_selection,'(',N_cmap,')']);
    end
    colormap(ax,main_colormap);
    % only update if called from menu item 
    % (otherwise being called from within update_plot itself)
    if ~isempty(hObject) 
        il_update_plot(0,curr_data_type);
    end
end

% Color Axis Callback:
function caxis_callback(~,~,type)
    if type==1 % Automatic
        caxis_auto = true;
    else
        prompt = {'Min Correlation:','Max Correlation:'};
        dlg_title = 'CAxis'; num_lines = [1,20;1,20]; 
        answer1 = inputdlg(prompt,dlg_title,num_lines);
        if isempty(answer1); disp('User cancelled action.'); return; end
        cmin = str2double(answer1{1});
        cmax = str2double(answer1{2});
        if cmin~=cmax && ~isnan(cmin) && ~isnan(cmax)
            caxis_auto = false;
        else
            warning('Color axis values should be numeric increasing. Using automatic caxis')
            caxis_auto = true;
        end
    end
    il_update_plot(0,curr_data_type,1);
end

function title_edit(~,~,~)
    set(get(ax,'Title'),'Editing','on'); 
end

function title_font(~,~,type)
    switch type
        case 1 % title_FontName
            dlg_title = 'Font';
            prompt = 'Title Font Name:'; num_lines = [1,20]; 
            answer = char(inputdlg(prompt,dlg_title,num_lines));
            if isempty(answer)
                disp('User cancelled action.'); 
                return; 
            else
                fonts_list = listfonts;
                if sum(strcmp(answer,fonts_list))
                    title_FontName = answer;
                else
                    error('Error: You must enter a valid FontName.')
                end
            end
            h_title.FontName = title_FontName;
        case 2 % title_FontSize
            dlg_title = 'Size';
            prompt = 'Title Font Size:'; num_lines = [1,20]; 
            answer = str2double(inputdlg(prompt,dlg_title,num_lines));
            if (rem(answer,1)==0) && (sign(answer)==1) % check if positive integer
                title_FontSize = answer;
            else
                error('Error: Please enter a positive integer value.')
            end
            h_title.FontSize = title_FontSize;
        case 3, h_title.FontWeight = 'normal';
        case 4, h_title.FontWeight = 'bold';
    end
    il_update_plot(0,curr_data_type);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function il_labels_callback(~,~,type)
    switch type
        case 1 % Labels On/Off
            if labels_on
                labels_on = false;
            else
                labels_on = true;
            end
        case 2 % FontName
            dlg_title = 'Font';
            prompt = 'Font Name:'; num_lines = [1,15]; 
            answer = char(inputdlg(prompt,dlg_title,num_lines));
            if isempty(answer)
                disp('User cancelled action.'); 
                return; 
            else
                fonts_list = listfonts;
                if sum(strcmp(answer,fonts_list))
                    FontName = answer;
                else
                    error('Error: You must enter a valid FontName.')
                end
            end
        case 3 % FontSize
            dlg_title = 'Size';
            prompt = 'Font Size:'; num_lines = [1,15]; 
            answer = str2double(inputdlg(prompt,dlg_title,num_lines));
            if (rem(answer,1)==0) && (sign(answer)==1) % check if positive integer
                FontSize = answer;
            else
                error('Error: Please enter a positive integer value.')
            end
        case 4, FontWeight = 'normal';
        case 5, FontWeight = 'bold';
    end
    il_update_plot(0,curr_data_type);
end

function il_colorbar_callback(~,~,type)
    switch type
        case 1 % Colorbar On/Off
            if colorbar_on
                colorbar_on = false;
                set(ax,'Position',ax_nocolorbar_pos);
            else
                colorbar_on = true;
                set(ax,'Position',ax_pos);
            end
            il_update_plot(0,curr_data_type);
        case 2 % FontName
            if colorbar_on
                dlg_title = 'Font';
                prompt = 'Font Name:'; num_lines = [1,15]; 
                answer = char(inputdlg(prompt,dlg_title,num_lines));
                if isempty(answer)
                    disp('User cancelled action.'); 
                    return; 
                else
                    fonts_list = listfonts;
                    if sum(strcmp(answer,fonts_list))
                        h_colorbar.FontName = answer;
                    else
                        error('Error: You must enter a valid FontName.')
                    end
                end
            else
                warning('Turn on colorbar before editing this property')
            end
        case 3 % FontSize
            if colorbar_on
                dlg_title = 'Size';
                prompt = 'Font Size:'; num_lines = [1,15]; 
                answer = str2double(inputdlg(prompt,dlg_title,num_lines));
                if (rem(answer,1)==0) && (sign(answer)==1) % check if positive integer
                    h_colorbar.FontSize = answer;
                else
                    error('Error: Please enter a positive integer value.')
                end
            else
                warning('Turn on colorbar before editing this property')
            end
        case 4 
            if colorbar_on 
                h_colorbar.FontWeight = 'normal'; 
            else
                warning('Turn on colorbar before editing this property')
            end
        case 5
            if colorbar_on
                h_colorbar.FontWeight = 'bold';
            else
                warning('Turn on colorbar before editing this property')
            end
    end
end

% Save Figure Function:
function il_save_figure_callback(~,~,~)
    % Identify File Extension:
    [title1,path1] = uiputfile('*.fig','Specify filename:');
    if isnumeric(path1) && path1==0
        disp('User cancelled printing.'); return;
    end
    figure_title = fullfile(path1,title1);
    savefig(h_corrmat,figure_title)
    disp(['Saved figure: ', figure_title])
end

function il_change_outline(hObject,~,~)
    if strcmp(hObject.Checked,'on')
        outline_menu.Checked = 'off';
        for jx = 1:length(h_line); delete(h_line(jx)); end
    else
        outline_menu.Checked = 'on';
        dim = size(corrmat{1},1);
        h_line(1) = plot([1.5,dim+.5],[0+.5;0+.5],'LineWidth',.05,'Color',zeros(1,3),'AlignVertexCenters','on','Parent',ax);
        for jx = 1:dim
            h_line(jx+1) = plot([jx+.5,dim+.5],[jx+.5;jx+.5],'LineWidth',.05,'Color',zeros(1,3),'AlignVertexCenters','on','Parent',ax);
            h_line(dim+jx+1) = plot([jx+.5;jx+.5],[-.5,jx+.5],'LineWidth',.05,'Color',zeros(1,3),'AlignVertexCenters','on','Parent',ax);
        end    
    end
end

end % end funct_view_spectrum

%% Bandpass Function
% The following function, used for applying band-pass filtration, is 
% provided from Brainstorm (http://neuroimage.usc.edu/brainstorm),
% distributed under the GNU license(http://www.gnu.org/copyleft/gpl.html).
% Author Credits: John Mosher, Francois Tadel, 2014

function [x,b,a] = il_bst_bandpass_filtfilt(x, Fs, HighPass, LowPass, isStopBand, FilterType)
% BST_BANDPASS_FILTFILT: Bandpass filter for the signal x, using the filtfilt function (used by default after Nov 2014)
%
% USAGE:  [x,b,a] = bst_bandpass_filtfilt(x, Fs, HighPass, LowPass, isStopBand=0, FilterType='fir')
% 
% INPUT: 
%    - x          : [nChannels,nTime] signal to process
%    - Fs         : Sampling frequency
%    - HighPass   : Frequency below this value are filtered (set to 0 for low-pass filter only)
%    - LowPass    : Frequency above this value are filtered (set to 0 for high-pass filter only)
%    - isStopBand : If 1, create a stop-band filter instead of a pass-band filter
%    - FilterType : 'fir' or 'iir'
%
% OUTPUT:
%    - x   : Filtered signals
%    - b,a : Filter coefficients, as defined in all the Matlab functions
%
% DESCRIPTION: 
%    - A linear phase FIR filter is created and applied both forward and backward using filtfilt.
%    - Function "filtfilt" is used to employ the filtering "mirror" trick, and to reset automatically the group delay.
%    - Function "kaiserord" and "kaiser" are used to set the necessary order for fir1. 
%    - The transition band is hard-coded. 
%    - Requires Signal Processing Toolbox for the following functions: kaiserord, kaiser, ellipord, 
%      If not, using Octave-based alternatives

% @=============================================================================
% This function is part of the Brainstorm software:
% http://neuroimage.usc.edu/brainstorm
% 
% Copyright (c)2000-2016 University of Southern California & McGill University
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPLv3
% license can be found at http://www.gnu.org/copyleft/gpl.html.
% 
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% For more information type "brainstorm license" at command prompt.
% =============================================================================@
%
% Authors: John Mosher, Francois Tadel, 2014

% ===== PARSE INPUTS =====
if (nargin < 6) || isempty(FilterType)
    FilterType = 'fir';
end
if (nargin < 5) || isempty(isStopBand)
    isStopBand = 0;
end
if isempty(HighPass)
    HighPass = 0;
end
if isempty(LowPass)
    LowPass = 0;
end
% If both high-pass and low-pass are zero: return signals unaltered
if (HighPass == 0) && (LowPass == 0)
    disp('BST_BANDPASS> Error: No frequency band in input');
    return;
end


% ===== FILTER PARAMETERS =====
PASSBAND_RIPPLE = 5;    % percent pass band ripple
PASSBAND_DB     = 1;    % dB of pass band ripple
STOP_ATTEN_DB   = 40;   % dB of attenuation in the stop band
TRANSITION_BAND = 0.05; % normalized to Nyquist, the allowed transition band
% We use filtfilt, which doubles the effect (attenuation and ripple)
PASSBAND_RIPPLE = PASSBAND_RIPPLE/2;
STOP_ATTEN_DB   = STOP_ATTEN_DB/2;
% Conversion from percent
Ripple = PASSBAND_RIPPLE/100;    
% Stop band attenuation
Atten  = 10^(-STOP_ATTEN_DB/20);

% Convert frequencies to normalized form
Nyquist = Fs/2;
f_highpass = HighPass / Nyquist;
f_lowpass  = LowPass  / Nyquist;
% Reasonable digital transition band
f_highstop = f_highpass - TRANSITION_BAND; 
f_lowstop  = f_lowpass  + TRANSITION_BAND;


% ===== CREATE FILTER =====
switch FilterType
    % ===== FIR =====
    case 'fir'
        % Build the general case first
        fcuts = [f_highstop, f_highpass, f_lowpass, f_lowstop];
        % Stop-band
        if isStopBand
            mags = [1 0 1];               % filter magnitudes
            devs = [Ripple Atten Ripple]; % deviations
        % Pass-band
        else
            mags = [0 1 0];               % filter magnitudes
            devs = [Atten Ripple Atten];  % deviations
        end
        % Now adjust for desired properties
        fcuts = max(0,fcuts);     % Can't go below zero
        fcuts = min(1-eps,fcuts); % Can't go above or equal to 1
        % We have implicitly created a bandpass, but now adjust for desired filter
        if (f_lowpass == 0)  % User didn't want a lowpass
            fcuts(3:4) = [];
            mags(3) = [];
            devs(3) = [];
        end
        if (f_highpass == 0)  % User didn't want a highpass
            fcuts(1:2) = [];
            mags(1) = [];
            devs(1) = [];
        end
        % Generate FIR filter
        [n,Wn,beta,ftype] = kaiserord(fcuts,mags,devs,2);
        n = n + rem(n,2);  % ensure even order
        b = fir1(n,Wn,ftype,kaiser(n+1,beta),'noscale');
        a = 1;
        
    % ===== IIR =====
    case 'iir'
        % Stop-band
        if isStopBand
            ftype = 'stop';
            Ws = [f_highpass f_lowpass]; % the range of stopped
            Wp = [f_highstop f_lowstop]; % the transition band
        % Pass-band
        else
            ftype = 'bandpass';
            Ws = [f_highstop f_lowstop]; % the transition band
            Wp = [f_highpass f_lowpass]; % the passband
        end
        % Now handle extremes
        Ws = max(eps,Ws);   % Can't be zero or less
        Wp = max(eps,Wp);
        Ws = min(1-eps,Ws); % Can't be 1
        Wp = min(1-eps,Wp);
        % Now handle highpass or lowpass only
        if (f_lowpass == 0)  % User didn't want a lowpass
            ftype = 'high';
            Ws(2) = [];
            Wp(2) = [];
        end
        if (f_highpass == 0)  % User didn't want a highpass
            ftype = 'low';
            Ws(1) = [];
            Wp(1) = [];
        end
        % Generate IIR filter
        [n,WP] = ellipord(Wp,Ws,PASSBAND_DB,STOP_ATTEN_DB);
        [b,a]  = ellip(n,PASSBAND_DB,STOP_ATTEN_DB,WP,ftype);
end


% ===== FILTER THE DATA =====
Ntime = size(x,2);
% Remove the mean of the data before filtering, which wrecks most filters
xmean = mean(x,2);
x = bsxfun(@minus, x, xmean)';    % Transposed output (time is now down the columns)
% Using filtfilt to use the mirroring trick and remove group delay
try
    x = filtfilt(b,a,x)';     % Transposed output
    
    % OCTAVE IMPLEMENTATION
    % http://octave-signal.sourcearchive.com/documentation/1.0.8/filtfilt_8m-source.html
catch 
    fprintf('Sequence too short for filtfilt, using alternate approach.\n')  
    xmirror = [flipud(x); x; flipud(x)];      % Mirror either end
    xmirror = filter(b,1,xmirror);            % Filter
    xmirror = flipud(xmirror);                % Reverse in time
    xmirror = flipud(filter(b,1,xmirror));    % Filter and flip again
    x = xmirror(Ntime + (1:Ntime),:)';        % Transposed output
end
% Restore the mean of the signal (only if there is no high-pass filter)
if (f_highpass == 0) 
    x = bsxfun(@plus, x, xmean);
end
end

%% REDBLUE colormap generation:

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
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function c = il_redblue(m)

if nargin < 1, m = size(get(gcf,'colormap'),1); end

if (mod(m,2) == 0) % if even
    % From [0 0 1] to [1 1 1], then [1 1 1] to [1 0 0];
    m1 = .5*m;
    r = (0:m1-1)'/max(m1-1,1);
    g = r;
    r = [r; ones(m1,1)];
    g = [g; flipud(g)];
    b = flipud(r);
else % if odd
    % From [0 0 1] to [1 1 1] to [1 0 0];
    m1 = floor(m*0.5);
    r = (0:m1-1)'/max(m1,1);
    g = r;
    r = [r; ones(m1+1,1)];
    g = [g; 1; flipud(g)];
    b = flipud(r);
end

c = [r g b];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nuisance Regression: perform_nuisance_regression(parsed_inputs.nuisance)
% function nuisance_regression
%     data_types(7) = 1;
%     for ixxx = 1:nSeries
%         X = parsed_inputs.nuisance{ixxx};
%         % Check for constant term
%         if ~any(all(X==1)); X = [ones(size(X,1),1),X]; end %#ok
%         for jxxx = 1:numrois
%             [~,~,alldata{7}{ixxx}(:,jxxx)] = regress(alldata{1}{ixxx}(:,jxxx),X);
%         end
%     end
% end

% Detrending:  perform_detrend(parsed_inputs.detrend)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function il_perform_detrend(type)
    switch type
        case 1 % parsed_inputs.detrend==1
            if ~data_types(2)
                disp('Performing linear detrending.')
                data_types(2) = 1;
                for ixxx = 1:nSeries
                    alldata{2}{ixxx} = detrend(alldata{1}{ixxx});
                end
            end
            if ~data_types(8)
                disp('Performing linear detrending on denoised data.')
                data_types(8) = 1;
                for ixxx = 1:nSeries
                    alldata{8}{ixxx} = detrend(alldata{7}{ixxx});
                end
            end
        case 2 % parsed_inputs.detrend==2
            if ~data_types(3)
                disp('Performing quadratic detrending.')
                x = (1:nTime)';
                for ixxx = 1:nSeries
                    for jxxx = 1:numrois
                        p = polyfit(x,alldata{1}{ixxx}(:,jxxx),2);
                        predicted = polyval(p,x);
                        alldata{3}{ixxx}(:,jxxx) = alldata{1}{ixxx}(:,jxxx)-predicted;
                    end
                end
                data_types(3) = 1;
            end
            if ~data_types(9)
                disp('Performing quadratic detrending on denoised data.')
                x = (1:nTime)';
                for ixxx = 1:nSeries
                    for jxxx = 1:numrois
                        p = polyfit(x,alldata{7}{ixxx}(:,jxxx),2);
                        predicted = polyval(p,x);
                        alldata{9}{ixxx}(:,jxxx) = alldata{7}{ixxx}(:,jxxx)-predicted;
                    end
                end
                data_types(9) = 1;
            end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bandpass
function il_perform_bandpass(initial)
    if all(parsed_inputs.bandpass==0)
        % Prompt User:
        prompt = {'Min Frequency (Hz):','Max Frequency (Hz):'};
        dlg_title = 'Band-Pass Settings'; num_lines = [1,40;1,40]; defaultans = {'.008','.1'};
        answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
        if isempty(answer); disp('User cancelled action.'); return; end
        if ~any(isnan(str2double(answer)))
            parsed_inputs.bandpass(1) = str2double(answer{1});
            parsed_inputs.bandpass(2) = str2double(answer{2});
        end
    end
    if parsed_inputs.bandpass(1) <= (.1*(fs/4))
        parsed_inputs.bandpass(1) = (.1*(fs/4))+.0001;
        warning('HighPass set too low for sampling frequency (Fs), adjusted to %.4g Hz',parsed_inputs.bandpass(1))
    end
    if parsed_inputs.bandpass(2) >= (fs/2)
        parsed_inputs.bandpass(2) = (fs/2)-.001;
        warning('LowPass set too high for sampling frequency (Fs), adjusted to %.4g Hz',parsed_inputs.bandpass(2))
    end
    for iter = find((data_types(1:3)-data_types(4:6))==1)
        for ixxx = 1:nSeries
            try
            [bandpassed,~,~] = il_bst_bandpass_filtfilt(alldata{iter}{ixxx}',fs,...
                parsed_inputs.bandpass(1), parsed_inputs.bandpass(2), 0, 'iir');
            catch
                error('Bandpass Filter Error: Check if ''Fs'' was specified accurately.')
            end
            alldata{iter+3}{ixxx} = bandpassed';
        end
        if ~initial
            calc_corrs(iter+3);
        end
        data_types(iter+3) = 1;
    end
    % Now repeat for nuisance regression:
    if nuisance_on
        iter2 = [12, 10, 11];
        ind_nuis = find((data_types([7, 8, 9])-data_types(iter2))==1);
        for iter = ind_nuis
            out_ind = iter2(iter);
            for ixxx = 1:nSeries
                try
                [bandpassed,~,~] = il_bst_bandpass_filtfilt(alldata{iter+6}{ixxx}',fs,...
                    parsed_inputs.bandpass(1), parsed_inputs.bandpass(2), 0, 'iir');
                catch
                    error('Bandpass Filter Error: Check if ''Fs'' was specified accurately.')
                end
                alldata{out_ind}{ixxx} = bandpassed';
            end
            if ~initial
                calc_corrs(out_ind);
            end
            data_types(out_ind) = 1;
        end
    end
end
