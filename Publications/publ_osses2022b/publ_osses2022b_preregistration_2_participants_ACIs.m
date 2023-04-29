function publ_osses2022b_preregistration_2_participants_ACIs(varargin)
% function publ_osses2022b_preregistration_2_participants_ACIs(varargin)
%
% 1. Description:
%
% 2. Examples:
%       publ_osses2022b_preregistration_2_participants_ACIs('subjectname','S01');
%
% Other functions:
%   l20220325_analysis_pipeline_newversion.m
%   l20211201_analysis_pipeline_versie_Alejandro.m on 13/01/2022
%
% Original name: g20220401_publ_osses2022b_ACIs (from fastACI_sim)
% Author: Leo Varnet, Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc, close all

bUse_idx_PA = 0;

if nargin == 0
    help publ_osses2022b_preregistration_2_participants_ACIs;
    return
end

h = [];
hname = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
definput.flags.plot={'plot','no_plot'};

definput.keyvals.subjectname=[];
definput.keyvals.dir_out=[];

[flags,keyvals]  = ltfatarghelper({},definput,varargin);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flags.do_plot
    bPlot = 1;
else
    bPlot = 0;
end

experiment = 'speechACI_Logatome-abda-S43M';
dir_where = fastACI_paths('dir_data');

if bPlot
    h = []; hname = [];
    dir_out_figs = fastACI_paths('dir_output');
end

if ~isempty(keyvals.subjectname)
    if ischar(keyvals.subjectname)
        Subjects = {keyvals.subjectname};
    elseif iscell(keyvals.subjectname)
        Subjects = keyvals.subjectname;
    end
end

if ~isempty(definput.keyvals.dir_out)
    bDir_out_default = 0;
    dir_out = '/home/alejandro/Desktop/fastACI_today/fastACI_dataproc/';
else
    bDir_out_default = 1;
    dir_out = fastACI_dir_datapost;
end
N_subjects = length(Subjects);

Maskers = {'white','bumpv1p2_10dB', 'sMPSv1p3'}; 
N_maskers  = length(Maskers); %     N_maskers  = length(Maskers);

maskercolors = {'b','r','g'};

%%% Analysis configuration:
trialtype_analysis = 'total';   % prereg, Sec. 6.1, after Eq. (1)
TF_type = 'gammatone';          % prereg, Sec. 6.2.1
expvar_after_reversal = 4;      % prereg, Sec. 6.1, after Eq. (1)
%%%

glmfct = 'l1glm';
switch glmfct
    case {'l1lm', 'l1glm'}
        N_lambda = 30;
        Lambdas = logspace(-4, -1, N_lambda);
        idx = find(Lambdas >= 10^-3);
        Lambdas = Lambdas(idx);
end

if bPlot
    hACI = figure; tiledlayout(N_subjects,N_maskers); % 1 x 3
    il_set_position(hACI,N_subjects,N_maskers)
    h(end+1) = hACI;
    hname{end+1} = ['ACI-' glmfct '-' trialtype_analysis];
 
    hMSEcurve = figure; hold on; % 1 x 1
    h(end+1) = hMSEcurve;
    hname{end+1} = ['MSE-' trialtype_analysis];

    hDEVcurve = figure; hold on; % 1 x 1
    h(end+1) = hDEVcurve;
    hname{end+1} = ['DEV-' trialtype_analysis];

    hPAcurve = figure; hold on; % 1 x 1
    h(end+1) = hPAcurve;
    hname{end+1} = ['PA-' trialtype_analysis];
end

lMSE = [];
lDEV = [];
lPA = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i_subject = 1:N_subjects
    subject = Subjects{i_subject};
    dir_subject = [fastACI_dir_data experiment filesep subject filesep];
    if bDir_out_default
        dir_exp = [dir_out experiment filesep];
        if ~exist(dir_exp,'dir')
            mkdir(dir_exp);
        end
        dir_exp = [dir_exp subject filesep];
        if ~exist(dir_exp,'dir')
            mkdir(dir_exp);
        end
        dir_out = dir_exp;
    else
        %%% Nothing to do
    end
    
    if bPlot
        switch subject
            case 'osses2021'
                YL_fig1_2 = [-30 -10];
            otherwise
                YL_fig1_2 = [-20 0];
        end
    end
    
    for i_masker = 1:N_maskers
        masker      = Maskers{i_masker};
        colour_here = maskercolors{i_masker};
        
        % loading data
        switch subject
            case {'osses2021','osses2022a','dau1997','king2019','maxwell2020','relanoiborra2019'}
                % Case for simulations:
                dirs = Get_filenames(dir_subject,'Run*');
                Show_cell(dirs);
                % bInput = input('Enter the number of the run you want to analyse: ');
                disp('the last folder will be used...')
                disp('pausing for 5 seconds, cancel (ctrl+c) if you want to abort');
                pause(5);
                bInput = length(dirs);
                
                dir_results = [dir_subject dirs{bInput} filesep 'Results' filesep];
                
            otherwise
                % dir_results = [dir_subject 'Results' filesep];
                dir_results = [fastACI_basepath 'Publications' filesep 'publ_osses2022b' filesep ...
                    'data_' subject filesep '1-experimental_results' filesep];
        end

        fname_results = Get_filenames(dir_results,['savegame_*' masker '.mat']);
        fname_results = fname_results{end};
        
        %%% Loading the data:
        [cfg_game,data_passation] = Convert_ACI_data_type([dir_results fname_results]);
        %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Prepare ACI analysis
        if ~isempty(dir_where)
            cfg_game = Check_cfg_crea_dirs(cfg_game,dir_where); 
        else
            cfg_game = Check_cfg_crea_dirs(cfg_game); % Updates dir_target/dir_noise to local folders
        end
        
        Data_matrix = [];
            
        flags_for_input = {TF_type, ...
            'trialtype_analysis',trialtype_analysis,... % prereg, Sec. 6.1, after Eq. (1)
            'N_folds', 10, ...                % prereg, Sec. 6.1, after Eq. (1)
            'dir_noise', cfg_game.dir_noise, ...
            'dir_target',cfg_game.dir_target, ...
            'add_signal',0, ...               % prereg, Sec. 6.1, after Eq. (1)
            'apply_SNR',0, ...                % prereg, Sec. 6.1, after Eq. (1)
            'skip_if_on_disk',1, ...          
            'no_permutation', ...
            'no_bias', ...                    % prereg, Sec. 6.1, after Eq. (1)
            'no_plot', ...                                                 % 'idx_trialselect',idx_trialselect,
            'expvar_after_reversal',expvar_after_reversal, ... % prereg, Sec. 6.1, after Eq. (1)
            'lambda',Lambdas, ...
            'dir_out',dir_out ...
            'pyramid_script','imgaussfilt', ... % new parameter as of 29/04
            'pyramid_shape', -1 ...             % new parameter as of 29/04
            };

        [ACI,cfg_ACI,results, Data_matrix] = fastACI_getACI([dir_results fname_results], glmfct, flags_for_input{:}, 'Data_matrix', Data_matrix);
        if bPlot
            handle_here = hACI;
        end
        
        Me = mean(results.FitInfo.PC_test,2); Me = 100*(Me-Me(end));
        [xx,idx2use] = max(Me);
        if idx2use ~= results.idxlambda
            fprintf('Using PA the ACI with idx=%.0f should be chosen instead of idx=%.0f\n',idx2use,results.idxlambda);                
        end    
        if bUse_idx_PA            
            ACI = squeeze(results.ACIs(idx2use,:,:));    
        end
        
        if bPlot
            % plot ACI lassoslow, classic revcorr, lasso
            figure(handle_here); nexttile; %(i_subject,i_masker)
            affichage_tf(ACI,'CI', 'cfg', cfg_ACI); hold on
            outs_aff = affichage_tf_add_Praat_metrics(cfg_game.dir_target,cfg_ACI,[], {'-','-.'},{[0.6,0,0],[0,0,0.6]},1.5);
            if i_masker ~= N_maskers
                set(outs_aff.hl,'Visible','off'); % Removes the legend
            end
            title(['ACI ' cfg_ACI.Condition ', ' trialtype_analysis], 'interpreter','none')
            set(handle_here,'Name',glmfct);
        end

        FI = results.FitInfo;
        %MSE curve
        x_var = FI.Lambda;
        y_var = FI.MSE_test - FI.MSE_test(end,:);
        idx_here = results.idxlambda;

        Me = mean(y_var,2);
        Std = std(FI.MSE_test,[],2);

        if bPlot
            figure(hMSEcurve)
            h_here = errorbar(x_var, Me, Std, 'Color',colour_here,'LineWidth',1); hold on, grid on

            plot(x_var(idx_here), Me(idx_here),'*','Color',colour_here,'LineWidth',2);
            lMSE(end+1) = h_here;

            if i_masker == N_maskers
                set(gca,'XScale','log'); 
                %ylim([0.2 0.35])
                title('MSE'); 
                xlabel('\lambda'); 
                ylabel('MSE');
                legend(lMSE, Maskers, 'interpreter', 'none')
            end
        end

        %DEV curve
        x_var = FI.Lambda;
        y_var = FI.Dev_test./FI.CV.TestSize;
        y_var_bias = y_var(end,:); % y_var = FI.Dev_test./FI.CV.TestSize - FI.Dev_test(end,:)./FI.CV.TestSize;
        y_var = y_var-y_var_bias;
        idx_here = results.idxlambda;

        Me = mean(y_var,2);
        Std = std(y_var,[],2);

        if bPlot
            figure(hDEVcurve)
            h_here = errorbar(x_var, Me, Std, 'Color',colour_here,'LineWidth',1); hold on, grid on

            plot(x_var(idx_here), Me(idx_here),'*','Color',colour_here,'LineWidth',2);
            lDEV(end+1) = h_here;

            if i_masker == N_maskers
                set(gca,'XScale','log'); 
                %ylim([0.2 0.35])
                title('Deviance/trial'); 
                xlabel('\lambda'); 
                ylabel('Dev/trial');
                legend(lDEV, Maskers, 'interpreter', 'none')
            end
        end
        
        % Percent accuracy curve
        % Same 'x_var', but new 'y_var':
        y_var_bias = FI.PC_test(end,:);
        y_var      = FI.PC_test-y_var_bias;
        Me     = 100*mean(y_var,2);
        y_bars = 100*std(y_var,[],2); %100*sem(y_var,[],2);

        if bPlot
            figure(hPAcurve)
            h_here = errorbar(x_var, Me, y_bars, 'Color', colour_here,'LineWidth',1); hold on, grid on
            plot(x_var(idx_here), Me(idx_here),'*','Color',colour_here,'LineWidth',2);

            lPA(end+1) = h_here;
            if i_masker == N_maskers
                set(gca,'XScale','log'); 
                ylim([-10 30])
                title('Percentage of responses explained by ACI'); 
                xlabel('\lambda'); 
                ylabel('% accuracy difference');
                legend(lPA, Maskers, 'interpreter', 'none')
            end
        end

    end
    
    if bPlot
        bSave = 1;

        if bSave == 1

            for i = 1:length(h)
                opts = [];
                opts.format = 'epsc';
                fname = [dir_out_figs 'Subj-' subject '-' hname{i}];
                Saveas(h(i),fname,opts);

                opts.format = 'png';
                Saveas(h(i),fname,opts);

                opts.format = 'fig';
                Saveas(h(i),fname,opts);
            end
        end
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function il_set_position(handle_fig,N_rows,N_cols)
% % Regular position:
% Pos34 = [ 570 420]; % for 1 x 1
% Pos34 = [ 800 420]; % for 1 x 2
% Pos34 = [1100 420]; % for 1 x 3
% Pos34 = [ 570 700]; % for 3 x 1

Pos = get(handle_fig,'Position');
switch N_rows
    case 1
        Pos(4) = 420;
    case 2
        % No figure with this config
    case 3
        Pos(4) = 700;
end
switch N_cols
    case 1
        Pos(3) = 570;
    case 2
        Pos(3) = 800;
    case 3
        Pos(3) = 1100;
end
set(handle_fig,'Position',Pos);
