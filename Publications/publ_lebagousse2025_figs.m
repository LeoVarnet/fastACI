function publ_lebagousse2025_figs(varargin)
% function publ_labagousse2025_figs(varargin)
%
% % To display Fig. 1 of Lebagousse & Varnet (2025) use :::
%     publ_labagousse2025_figs('fig1');
%
% fig1: see the original publication by Ahumada Marken and Sandusky (AMS75)
% fig2: Results obtained for the 4 participants in the 11.8-dB condition (SNR matched with AMS75).
% fig3: Results obtained for the 9 participants in the 18-dB condition (performance matched with AMS75).
% fig4A: mean ACI in the in the 18-dB condition
% fig4B: mean ACI in the in the 11.8-dB condition
% fig5A: Marginalized target-absent ACIs along the temporal and spectral dimensions (18-dB condition)
% fig5B: Marginalized target-absent ACIs along the temporal and spectral dimensions (11.8-dB condition)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    help publ_carranante2024_figs;
    return
end

h = [];
hname = [];

definput.flags.type={'missingflag','fig1','fig3','fig2','fig4A','fig4B','fig5A','fig5B'};
% definput.keyvals.models=[];
definput.keyvals.dir_out=[];
definput.flags.local={'local','zenodo'};
definput.keyvals.dir_zenodo=[];

[flags,keyvals]  = ltfatarghelper({},definput,varargin);

bZenodo = flags.do_zenodo;
if bZenodo == 1
    fprintf('%s: The data from Zenodo will be loaded... \n',mfilename);
    if ~isempty(keyvals.dir_zenodo)
        fprintf('\t The data from Zenodo will be loaded from %s \n',keyvals.dir_zenodo);
        if exist(keyvals.dir_zenodo,'dir')
            bIs_dir = 1;
        else
            bIs_dir = 0;
        end
    else
        bIs_dir = 0;
    end
    if bIs_dir == 0
        msg_here = sprintf('\t You indicated that you want to proccess the data obtained from Zenodo. Please indicate the directory where it is...\n');
        fprintf(msg_here)
        keyvals.dir_zenodo = uigetdir([pwd filesep],msg_here);
        keyvals.dir_zenodo = [keyvals.dir_zenodo filesep];
    end
end

if flags.do_fig1
    error('Figure adapted from Ahumada, Marken & Sandusky (1975). Please refer to the original publication.')
end

if flags.do_fig3 || flags.do_fig4A || flags.do_fig5A
    Subjects = {'SQ01','SQ02','SQ03','SQ04','SQ05','SQ06','SQ07','SQ08','SQ09'};%
elseif flags.do_fig2 || flags.do_fig4B || flags.do_fig5B
    Subjects = {'S01','S03','S05','S09'};
else 
    Subjects = {};
end

if flags.do_fig3 || flags.do_fig2
    script_dataload = 'replication_ahumada1975_dataload';
    dir_out = ['ResultsAhumada' filesep];
elseif flags.do_fig4A || flags.do_fig4B || flags.do_fig5A || flags.do_fig5B
    script_dataload = 'fastACI_getACI_dataload';
    dir_out = ['ResultsACI/' filesep];
else 
    Subjects = {};
end

glmfct = 'classic_revcorr';%'l1glm';%

masker = 'white';

if flags.do_fig3 || flags.do_fig2 || flags.do_fig4A || flags.do_fig4B
    trialtype_analysis = 'total';
else
    trialtype_analysis = 't2';
end

TF_type = 'gammatone';
N_subjects = length(Subjects);

resume_trial = [1,401,801,1201,1601,2001,2401,2801,3201];%resume_trial = data_passation.resume_trial; resume_trial(1)=1;resume_trial(end+1)=data_passation.i_current;
N_sessions  = length(resume_trial)-1;

flags_for_input = {TF_type,...
    'trialtype_analysis', trialtype_analysis,...
    'N_folds', 10, ...
    'bwmul',0.25,...
    'no_plot', ...
    'add_signal',0, ...
    'apply_SNR',0, ...
    'script_dataload',script_dataload,...
    'skip_if_on_disk',1, ...
    'no_permutation', ...,
    'dir_out',dir_out,...
    };

switch glmfct
    case 'l1glm'
        N_lambda = 30;
        Lambdas = logspace(-4, -1, N_lambda);
        idx = find(Lambdas >= 10^-3);
        Lambdas = Lambdas(idx);
        flags_for_input = [flags_for_input, ...
            'lambda',Lambdas, ...
            'pyramid_script','imgaussfilt', ...
            'pyramid_shape',-1, ...
            ];
    case 'classic_revcorr'
        
    otherwise
        error('For this analysis, only the ''l1glm'' and ''classic_revcorr'' were used.\n')
end

%dir_out = '';
%dir_out_figs = [pwd filesep];
fullpath = fastACI_paths('dir_data');

%% Generate targets

if flags.do_fig1_suppl
        
    cfg_inout = replication_ahumada1975_set([]);
    
    basef = 500;
    flags_gamma = {'basef',basef,'flow',1,'fhigh',5000,'bwmul',0.05,'dboffset',100,'no_adt','binwidth',0.01, ...
        'no_outerear','no_middleear'};
    
    [G_tone, fc, t_spec, outs] = Gammatone_proc(cfg_inout.signal, cfg_inout.fs, flags_gamma{:});
    cfg_ACI.t = t_spec;
    cfg_ACI.f = fc;
    % 
    % [~,color_aba,color_apa] = selectcolormap('speechACI_Logatome-abpa-S43M');
    % [~,color_ada,color_ata] = selectcolormap('speechACI_Logatome-adta-S43M');
    % [~,color_aga,color_aka] = selectcolormap('speechACI_Logatome-agka-S43M');
end

%% suppl Fig 1. Stimuli

if flags.do_fig1_suppl
    h_fig1 = figure('Position', [100,100,700,450]); 
    affichage_tf(log(G_tone)', 'pow', t_spec, fc); %title('/aka/', 'Color', color_aka,'FontSize', fontSize)
    %hold on; affichage_tf_add_Praat_metrics_one_sound(path_aka, cfg_ACI, [], '-', 'w')
    %set(gca,'YTickLabels',{}); colorbar off
    xlabel('time (s)'); ylabel('frequency (Hz)');
    clim([-7 -4.8])
colorbar off
end


%% Load ACI data

if flags.do_fig3 || flags.do_fig2 || flags.do_fig4A || flags.do_fig4B || flags.do_fig5A || flags.do_fig5B

        for i_subject = 1:N_subjects
            
            subject = Subjects{i_subject};
            
            dir_subject = [fastACI_dir_data 'replication_ahumada1975' filesep subject filesep];
            dir_results = [dir_subject 'Results' filesep];
            
            fname_results = Get_filenames(dir_results,['savegame_*' masker '.mat']);
            
            if ~isempty(fname_results)
                fname_results = fname_results{end};
                fprintf(['Processing ' fname_results '\n'])
                load([dir_results fname_results]);
                cd(dir_results)
                
                cfg_games{i_subject} = Check_cfg_crea_dirs(cfg_game); % Updates dir_target/dir_noise to local folders
                
                Data_matrix = [];
                
                flags_here = {'dir_noise',cfg_games{i_subject}.dir_noise, ...
                    'skip_if_on_disk',1};
                
                switch glmfct
                    case {'lassoslow','lassoglmslow','l1glm'}
                        % calculate ACI lassoslow
                        [ACI,cfg_ACI{i_subject},results{i_subject}, Data_matrix] = fastACI_getACI(fname_results, glmfct, flags_for_input{:}, flags_here{:}, 'Data_matrix', Data_matrix);
                        
                        lambda = results{i_subject}.FitInfo.Lambda;
                        Dev_test(i_subject,:,:) = results{i_subject}.FitInfo.Dev_test./results{i_subject}.FitInfo.CV.TestSize;
                        PC_test(i_subject,:,:) = results{i_subject}.FitInfo.PC_test;
                        idxlambda(i_subject) = results{i_subject}.idxlambda;
                        Nfolds = results{i_subject}.FitInfo.CV.NumTestSets;
                        
                        if flags.do_fig2 || flags.do_tab1
                        % statistical thresholding 
                        % The thresholds are dependent on idxlambda. They
                        % were obtained by randomization (based on 1000
                        % permutations). 
                            B = results{i_subject}.B(:,results{i_subject}.idxlambda);
                            
                            switch idxlambda(i_subject)
                                case 12
                                    B(B<0.0240 & B>-0.0243) = 0;
                                case 13
                                    B(B<0.0142 & B>-0.0138) = 0;
                                case 14
                                    B(B<0.0043 & B>-0.0041) = 0;
                            end
                            [ACI] = Convert_lasso_B2ACI(B, cfg_ACI{i_subject}, 1, cfg_ACI{i_subject}.keyvals);
                        end
                    otherwise
                        % calculate ACI classic revcorr
                        [ACI,cfg_ACI{i_subject},results{i_subject}] = fastACI_getACI(fname_results, glmfct, flags_for_input{:}, 'Data_matrix', Data_matrix);
                end
                ACI_matrix(:,:,i_subject) = ACI;
                Subject_matrix{i_subject} = subject;
                i_subject = i_subject+1;
            end
        end

end

%% fig3. and Fig 3. Ahumada's like ACIs for the two groups

if flags.do_fig3 || flags.do_fig2
            
    if flags.do_fig3 
        hfig2 = figure('Position', [10,10,680,500]);
        tiledlayout(2,ceil(N_subjects/2),'TileSpacing', 'compact');
        displayfactor = 400;
    else
        hfig2 = figure('Position', [10,10,550,250]);
        tiledlayout(1,N_subjects,'TileSpacing', 'compact');
        displayfactor = 200;
    end
        
        [cmap, color1, color2] = selectcolormap('speechACI_Logatome-abda-S43M');
        
        for i_subject = 1:N_subjects
            
            % plot ACI lassoslow, classic revcorr, lasso
            %idx_tile = i_condition+(i_subject-1)*N_conditions;
            figure(hfig2); nexttile;%(idx_tile); %(i_subject,i_masker)

            % ACI-like display
            % affichage_tf(ACI_matrix(:,:,i_subject),'CI', 'cfg', cfg_ACI, 'colorbar_map', cmap, 'NfrequencyTicks', 4); hold on

            %Ahumada-like display
            plot([0.05 0.15 0.25 0.35 0.45],ACI_matrix(:,:,i_subject)'*displayfactor+[400 450 500 550 600],'k','LineWidth',2);hold on
            plot([0.05 0.15 0.25 0.35 0.45],zeros([5,5])+[400 450 500 550 600],'k-+','LineWidth',0.5,'MarkerSize',5)
            xlim([0 0.5]); set(gca,'XTick',[0 0.1 0.2 0.3 0.4 0.5]);
            ylim([350 650]);
            set(gca,'YTick',[400 450 500 550 600]);
            xlabel('time (s)'); ylabel('frequency (Hz)'); 

            %text(0.95,0.95,Subject_matrix(i_subject),'Units','normalized','HorizontalAlignment', 'right','VerticalAlignment', 'top')
            title(Subject_matrix(i_subject))
            colorbar off
            %ylabel('');set(gca,'YTickLabels',[]);
            if ~(i_subject == 1 || i_subject == 6)
                ylabel('');
                set(gca,'YTickLabels',[]);
            end
            if flags.do_fig3 && i_subject < 6
                xlabel('');
                set(gca,'XTickLabels',[]);
            end
            box off
        end
            
end

%% Fig.4 High-resolution ACIs for the two groups

if flags.do_fig4A || flags.do_fig4B
            
        hfig4 = figure('Position', [10,10,300,300]);

        [cmap, color1, color2] = selectcolormap('speechACI_Logatome-abda-S43M');

        figure(hfig4); nexttile;%(idx_tile); %(i_subject,i_masker)

        % ACI-like display
        %affichage_tf(ACI_matrix(:,:,i_subject),'CI', 'cfg', cfg_ACI{i_subject}, 'colorbar_map', cmap, 'NfrequencyTicks', 4); hold on
        affichage_tf(mean(ACI_matrix,3),'CI', 'cfg', cfg_ACI{1}, 'colorbar_map', cmap, 'NfrequencyTicks', 7); hold on
        set(gca,'YTickLabels',round(str2num(get(gca,'YTickLabels'))/10,0)*10)
        %text(0.95,0.95,Subject_matrix(i_subject),'Units','normalized','HorizontalAlignment', 'right','VerticalAlignment', 'top')
        if flags.do_fig4A
            title('18-dB condition')
        else
            title('11.8-dB condition')
        end
        colorbar off
end

%% Fig.5 High-resolution ACIs, target-absent

if flags.do_fig5A || flags.do_fig5B
            
        hfig4 = figure('Position', [10,10,600,300]);

        %[cmap, color1, color2] = selectcolormap('speechACI_Logatome-abda-S43M');

        figure(hfig4); tiledlayout(1,2);

        % ACI-like display
        %affichage_tf(ACI_matrix(:,:,i_subject),'CI', 'cfg', cfg_ACI{i_subject}, 'colorbar_map', cmap, 'NfrequencyTicks', 4); hold on
        % affichage_tf(mean(ACI_matrix,3),'CI', 'cfg', cfg_ACI{1}, 'colorbar_map', cmap, 'NfrequencyTicks', 7); hold on
        % set(gca,'YTickLabels',round(str2num(get(gca,'YTickLabels'))/10,0)*10)
        %text(0.95,0.95,Subject_matrix(i_subject),'Units','normalized','HorizontalAlignment', 'right','VerticalAlignment', 'top')

        idx_f = cfg_ACI{1, 1}.f>=450 & cfg_ACI{1, 1}.f<=550;
        idx_t = cfg_ACI{1, 1}.t>=0.2 & cfg_ACI{1, 1}.t<=0.3;

        nexttile(2); 
        plot(cfg_ACI{1, 1}.t,squeeze(mean(mean(ACI_matrix(idx_f,:,:),1),3)))
        xlabel('time (s)')
        ylabel('correlation')
        hold on; plot([0 0.5],[0 0],'k--')
        ylim([-0.1 0.1])
        xlim([0 0.5])
        title({'average correlation' 'in the 450-550-Hz band'})

        nexttile(1); 
        semilogx(cfg_ACI{1, 1}.f,squeeze(mean(mean(ACI_matrix(:,idx_t,:),2),3)))
        xlabel('frequency (Hz)')
        ylabel('correlation')
        hold on; plot([50 5000],[0 0],'k--')
        ylim([-0.01 0.06])
        xlim([50 5000])
        title({'average correlation', 'in the 0.2-0.3-s segment'})
end

%% Supplementary Fig 2. Diference of spectrograms

if flags.do_fig3_suppl

    for i_condition = 1:N_conditions
        switch i_condition
            case 1
                cmap = 'DG_red_blue';
                spec_diff = G_aba'-G_ada';
            case 2
                cmap = 'DG_blue_yellow';
                spec_diff = G_ada'-G_aga';
            case 3
                cmap = 'DG_red_green';
                spec_diff = G_apa'-G_ata';
            case 4
                cmap = 'DG_green_purple';
                spec_diff = G_aka'-G_ata';
            case 5
                cmap = 'DG_blue_purple';
                spec_diff = G_aba'-G_apa';
            case 6
                cmap = 'DG_cyan_purple';
                spec_diff = G_ada'-G_ata';
            case 7
                cmap = 'DG_yellow_cyan';
                spec_diff = G_aga'-G_aka';
        end

        hfig3_suppl{i_condition} = figure('Position', [100,100,400,300]);
        experiment = Conditions{i_condition};
        [cmap, color1, color2] = selectcolormap(experiment);

        affichage_tf(spec_diff,'CI', 'cfg', cfg_ACI, 'colorbar_map', cmap); hold on

        %outs_aff = affichage_tf_add_Praat_metrics(cfg_games{i_condition}.dir_target,cfg_ACI);%,[], {'-','-'},{color1,color2},1.5);

        title(upper(experiment(20:23)), 'interpreter','none');
        %cfg_ACI{i_condition}.keyvals.trialtype_analysis

    end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handle_tl = il_tiledlayout(N,M,TileSpacing,TileSpacing_option,varargin)

handle_tl = [];
if nargin < 4
    TileSpacing_option = 'Compact';
end
bExist = exist('tiledlayout','file'); % tiledlayout.p
if bExist
    if isempty(varargin)
        handle_tl = tiledlayout(N,M,TileSpacing,TileSpacing_option);
    else
        handle_tl = tiledlayout(N,M,TileSpacing,TileSpacing_option,varargin{:});
    end
else
    warning('We programmed this script to use more recent graphic options from MATLAB. It might be that you won''t be able to visualise these results as we have foreseen...')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ax = il_nexttile(N)

bExist = exist('nexttile','file'); % nexttile
if bExist
    if nargin == 0
        nexttile;
    else
        nexttile(N);
    end
else
    if N == 1
        close;
    end
    figure;
end

if nargout ~= 0
    ax = gcf;
end

function [cmap,color1,color2] = selectcolormap(experiment)
switch experiment
    case 'speechACI_Logatome-abda-S43M'
        cmap = 'DG_red_blue';
        color1 = [0.6,0,0];
        color2 = [0,0,0.6];
    case 'speechACI_Logatome-adga-S43M'
        cmap = 'DG_blue_yellow';
        color1 = [0,0,0.6];
        color2 = [0.6,0.6,0];
    case 'speechACI_Logatome-abpa-S43M'
        cmap = 'DG_red_green';
        color1 = [0.6,0,0];
        color2 = [0,0.6,0];
    case 'speechACI_Logatome-apta-S43M'
        cmap = 'DG_green_purple';
        color1 = [0,0.6,0];
        color2 = [0.6,0,0.6];
    case 'speechACI_Logatome-adta-S43M'
        cmap = 'DG_blue_purple';
        color1 = [0,0,0.6];
        color2 = [0.6,0,0.6];
    case 'speechACI_Logatome-akta-S43M'
        cmap = 'DG_cyan_purple';
        color1 = [0,0.6,0.6];
        color2 = [0.6,0,0.6];
    case 'speechACI_Logatome-agka-S43M'
        cmap = 'DG_yellow_cyan';
        color1 = [0.6,0.6,0];
        color2 = [0,0.6,0.6];
end

function colors = distinguishable_colors(n_colors,bg,func)
% DISTINGUISHABLE_COLORS: pick colors that are maximally perceptually distinct
%
% When plotting a set of lines, you may want to distinguish them by color.
% By default, Matlab chooses a small set of colors and cycles among them,
% and so if you have more than a few lines there will be confusion about
% which line is which. To fix this problem, one would want to be able to
% pick a much larger set of distinct colors, where the number of colors
% equals or exceeds the number of lines you want to plot. Because our
% ability to distinguish among colors has limits, one should choose these
% colors to be "maximally perceptually distinguishable."
%
% This function generates a set of colors which are distinguishable
% by reference to the "Lab" color space, which more closely matches
% human color perception than RGB. Given an initial large list of possible
% colors, it iteratively chooses the entry in the list that is farthest (in
% Lab space) from all previously-chosen entries. While this "greedy"
% algorithm does not yield a global maximum, it is simple and efficient.
% Moreover, the sequence of colors is consistent no matter how many you
% request, which facilitates the users' ability to learn the color order
% and avoids major changes in the appearance of plots when adding or
% removing lines.
%
% Syntax:
%   colors = distinguishable_colors(n_colors)
% Specify the number of colors you want as a scalar, n_colors. This will
% generate an n_colors-by-3 matrix, each row representing an RGB
% color triple. If you don't precisely know how many you will need in
% advance, there is no harm (other than execution time) in specifying
% slightly more than you think you will need.
%
%   colors = distinguishable_colors(n_colors,bg)
% This syntax allows you to specify the background color, to make sure that
% your colors are also distinguishable from the background. Default value
% is white. bg may be specified as an RGB triple or as one of the standard
% "ColorSpec" strings. You can even specify multiple colors:
%     bg = {'w','k'}
% or
%     bg = [1 1 1; 0 0 0]
% will only produce colors that are distinguishable from both white and
% black.
%
%   colors = distinguishable_colors(n_colors,bg,rgb2labfunc)
% By default, distinguishable_colors uses the image processing toolbox's
% color conversion functions makecform and applycform. Alternatively, you
% can supply your own color conversion function.
%
% Example:
%   c = distinguishable_colors(25);
%   figure
%   image(reshape(c,[1 size(c)]))
%
% Example using the file exchange's 'colorspace':
%   func = @(x) colorspace('RGB->Lab',x);
%   c = distinguishable_colors(25,'w',func);
% Copyright 2010-2011 by Timothy E. Holy
  % Parse the inputs
  if (nargin < 2)
    bg = [1 1 1];  % default white background
  else
    if iscell(bg)
      % User specified a list of colors as a cell aray
      bgc = bg;
      for i = 1:length(bgc)
	bgc{i} = parsecolor(bgc{i});
      end
      bg = cat(1,bgc{:});
    else
      % User specified a numeric array of colors (n-by-3)
      bg = parsecolor(bg);
    end
  end
  
  % Generate a sizable number of RGB triples. This represents our space of
  % possible choices. By starting in RGB space, we ensure that all of the
  % colors can be generated by the monitor.
  n_grid = 30;  % number of grid divisions along each axis in RGB space
  x = linspace(0,1,n_grid);
  [R,G,B] = ndgrid(x,x,x);
  rgb = [R(:) G(:) B(:)];
  if (n_colors > size(rgb,1)/3)
    error('You can''t readily distinguish that many colors');
  end
  
  % Convert to Lab color space, which more closely represents human
  % perception
  if (nargin > 2)
    lab = func(rgb);
    bglab = func(bg);
  else
    C = makecform('srgb2lab');
    lab = applycform(rgb,C);
    bglab = applycform(bg,C);
  end
  % If the user specified multiple background colors, compute distances
  % from the candidate colors to the background colors
  mindist2 = inf(size(rgb,1),1);
  for i = 1:size(bglab,1)-1
    dX = bsxfun(@minus,lab,bglab(i,:)); % displacement all colors from bg
    dist2 = sum(dX.^2,2);  % square distance
    mindist2 = min(dist2,mindist2);  % dist2 to closest previously-chosen color
  end
  
  % Iteratively pick the color that maximizes the distance to the nearest
  % already-picked color
  colors = zeros(n_colors,3);
  lastlab = bglab(end,:);   % initialize by making the "previous" color equal to background
  for i = 1:n_colors
    dX = bsxfun(@minus,lab,lastlab); % displacement of last from all colors on list
    dist2 = sum(dX.^2,2);  % square distance
    mindist2 = min(dist2,mindist2);  % dist2 to closest previously-chosen color
    [~,index] = max(mindist2);  % find the entry farthest from all previously-chosen colors
    colors(i,:) = rgb(index,:);  % save for output
    lastlab = lab(index,:);  % prepare for next iteration
  end
%end

function c = parsecolor(s)
  if ischar(s)
    c = colorstr2rgb(s);
  elseif isnumeric(s) && size(s,2) == 3
    c = s;
  else
    error('MATLAB:InvalidColorSpec','Color specification cannot be parsed.');
  end
%end

function c = colorstr2rgb(c)
  % Convert a color string to an RGB value.
  % This is cribbed from Matlab's whitebg function.
  % Why don't they make this a stand-alone function?
  rgbspec = [1 0 0;0 1 0;0 0 1;1 1 1;0 1 1;1 0 1;1 1 0;0 0 0];
  cspec = 'rgbwcmyk';
  k = find(cspec==c(1));
  if isempty(k)
    error('MATLAB:InvalidColorString','Unknown color string.');
  end
  if k~=3 || length(c)==1,
    c = rgbspec(k,:);
  elseif length(c)>2,
    if strcmpi(c(1:3),'bla')
      c = [0 0 0];
    elseif strcmpi(c(1:3),'blu')
      c = [0 0 1];
    else
      error('MATLAB:UnknownColorString', 'Unknown color string.');
    end
  end
%end
