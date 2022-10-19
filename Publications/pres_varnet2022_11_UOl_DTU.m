function data = pres_varnet2022_11_UOl_DTU(varargin)
% function data = pres_varnet2022_11_UOl_DTU(varargin)
%
% Generates the figures
%
% % To display Fig. 1 of Osses and Varnet, (2022, JASA) use :::
%     pres_varnet2022_11_UOl_DTU('fig8'); % Cross predictions, experimental
%
% % To display Fig. 2 of Osses and Varnet, (2022, JASA) use :::
%     pres_varnet2022_11_UOl_DTU('fig8b'); % Cross predictions, simulations
%
% Author: Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% close all, clc

if nargin == 0
    help pres_varnet2022_11_UOl_DTU;
    return
end

h = [];
hname = [];

definput.flags.type={'missingflag', ...
    'fig8', ... % Cross predictions, L1608
    'fig8b', ... % Cross predictions 'between participant' for simulations
    'fig9', ... % Cross predictions between noise
    'fig12' ... % Cross predictions between noise for simulations
}; % different lambdas
definput.flags.plot={'plot','no_plot'};
% definput.keyvals.models=[];
definput.keyvals.dir_out=[];

[flags,keyvals]  = ltfatarghelper({},definput,varargin);

% dir_fastACI_results = fastACI_paths('dir_data'); % '/home/alejandro/Documents/Databases/data/fastACI/';
% 
% experiment = 'speechACI_varnet2013';
% dir_exp = [dir_fastACI_results experiment filesep];

% Common variables:
noise_str         = {'white','bumpv1p2_10dB','sMPSv1p3'};
noise_types_label = {'white','bump','MPS'};
masker_colours    = {[0,0,1],[1,0,0],rgb('Green')}; % Leo's colours

Show_cell(noise_types_label);
bIdx = input('Choose the noises to be processed ([1 2 3], or, 1, 2, or 3): ');
noise_str = noise_str(bIdx);
noise_types_label = noise_types_label(bIdx);
masker_colours = masker_colours(bIdx);

experiment        = 'speechACI_Logatome-abda-S43M';
N_maskers = length(noise_str);

% End common variables

dir_savegame_sim = [fastACI_dir_data 'speechACI_Logatome-abda-S43M' filesep 'osses2022a_debug_0715' filesep]; % '/home/alejandro/Documents/Databases/data/fastACI_data/speechACI_Logatome-abda-S43M/osses2022a_debug_0715/';
dir_savegame_exp = [fastACI_basepath 'Publications' filesep 'publ_osses2022b' filesep]; 

dir_ACI_exp = '/home/alejandro/Desktop/fastACI_today/fastACI_dataproc_Leo/'; % My run: dir_where = '/home/alejandro/Desktop/fastACI_today/fastACI_dataproc/';
    % Run-S01/ACI-osses2022a-speechACI_Logatome-white-nob-gt-l1glm-rev4.mat
    % Run-S02/ACI-osses2022a-speechACI_Logatome-white-nob-gt-l1glm-rev4.mat
    % ...
    % Run-S12/ACI-osses2022a-speechACI_Logatome-white-nob-gt-l1glm-rev4.mat
dir_ACI_sim = '/home/alejandro/Documents/Databases/data/fastACI_data_z_tmp/20220715_sim_Q1_osses2022a/';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

if flags.do_fig8 || flags.do_fig9

    %  g20220720_behavstats_from_l20220325
    %%% Analysis for 12 participants:
    Subjects = {'S01','S02','S03','S04','S05','S06','S07','S08','S09','S10','S11','S12'}; 
    N_subjects = length(Subjects);
end

if flags.do_fig12
    Subjects = {'S01','S02','S03','S04','S05','S06','S07','S08','S09','S10','S11','S12'}; 
    N_subjects = length(Subjects);
end

if flags.do_fig8 || flags.do_fig9
    glmfct = 'l1glm';
    bIs_Alejandro = isunix;
    
    flags_base = {'gammatone', ... % gt
                  glmfct, ... % l1glm
                  'no_bias', ... % nob
                  'expvar_after_reversal',4, ... %'idx_trialselect', 1:4000 ... % tr4000
                  'pyramid_script','imgaussfilt', ...
                  'pyramid_shape',-1 ...
                  }; % the extra fields

    if bIs_Alejandro
        dir_where = dir_ACI_exp;
        % Location of the savegames:
                   % /home/alejandro/Desktop/fastACI_today/fastACI_data/speechACI_Logatome-abda-S43M/S01/Results
        % dir_res = {'/home/alejandro/Documents/MATLAB/MATLAB_ENS/fastACI/Publications/publ_osses2022b/',['1-experimental_results' filesep]};
        % dir_wav = '/home/alejandro/Documents/Databases/data/fastACI_data/speechACI_Logatome-abda-S43M/';
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flags.do_fig8 || flags.do_fig9
    % Loading cross predictions
    
    dir_res = {[fastACI_paths('dir_fastACI') 'Publications' filesep 'publ_osses2022b' filesep],['1-experimental_results' filesep]};
           
    trialtype_analysis = 'total';

    
    for i_masker = 1:N_maskers
        masker  = noise_str{i_masker};

        for i_subject = 1:N_subjects
         
            subject = Subjects{i_subject};
            
            %%%
            if flags.do_fig8 || flags.do_fig9
                % 2.1. Checking whether you previously stored already the cross-prediction data:
                if flags.do_fig8 || flags.do_fig3_suppl || flags.do_fig8b_old
                    if i_subject == 1 
                        bInclude = zeros(1,N_subjects);
                        cross_suffix = '';
                        fname_crosspred = [];
                        for ii = 1:N_subjects
                            fname_crosspred{ii} = [dir_ACI_exp 'ACI-' Subjects{ii} '-speechACI_Logatome-' masker '-nob-gt-' glmfct '+pyrga-rev4.mat'];
                            if exist(fname_crosspred{ii},'file')
                                bInclude(ii) = 1;
                            end    
                        end
                    end
                else
                    bInclude = zeros(1,N_subjects);
                    
                    cross_suffix = 'noise';
                    for ii = 1:N_maskers
                        % fname_crosspred{ii} = [dir_subject 'Results' filesep 'Results_ACI' filesep 'ACI-' subject '-speechACI_Logatome-' noise_str{ii} '-nob-gt-l1glm-rev4.mat'];
                        fname_crosspred{ii} = [dir_ACI_exp 'ACI-' subject '-speechACI_Logatome-' noise_str{ii} '-nob-gt-' glmfct '+pyrga-rev4.mat'];
                        if ii == i_masker
                            if exist(fname_crosspred{ii},'file')
                                bInclude(i_subject) = 1;
                            end
                        end
                    end
                end
            end
            
            if bInclude(i_subject)
                % loading data
                dir_results = [dir_res{1} 'data_' subject filesep dir_res{2}];
                files = Get_filenames(dir_results,['savegame_*' masker '.mat']);
                fname_results = [dir_results files{1}];

                [cfg_game,data_passation] = Convert_ACI_data_type(fname_results);

                dir_out_here = dir_where; %[dir_where subject filesep];

                cfg_game = Check_cfg_crea_dirs(cfg_game); % probably gives an error if the folders aren't visible

                flags_for_input = {...
                    'trialtype_analysis', trialtype_analysis, ...
                    'N_folds', 10, ...
                    'dir_noise',cfg_game.dir_noise, ...
                    'dir_target',cfg_game.dir_target, ...
                    'add_signal',0, ...
                    'apply_SNR',0, ...
                    'skip_if_on_disk',1, ...%'force_dataload', ...
                    'no_permutation', ...
                    'no_bias', ...
                    'no_plot'...
                    'dir_out',dir_out_here ...
                    };

                %%% 1. Obtaining the ACI of the current subject:
                [col1,cfg_ACI_here,res,Data_matrix_here] = fastACI_getACI(fname_results, ...
                    flags_base{:},flags_for_input{:}); % ,'force_dataload');
                % crosspred = [];
                
                % 2.1. Checking whether you previously stored already the cross-prediction data:
                [crosspred, fname_cross] = Check_if_crosspred(cfg_ACI_here.fnameACI,fname_crosspred,cross_suffix);

                bForceComplete = 0;

                if isempty(crosspred) || bForceComplete
                    % 2.2. If there is no cross-prediction in your computer, then
                    %      it obtains it
                    flags = {'ACI_crosspred',fname_crosspred,'Data_matrix',Data_matrix_here}; % the extra fields
                    [col1,cfg_ACI_here,res,col4] = fastACI_getACI(fname_results, flags_base{:},flags{:},flags_for_input{:});
                else
                    % (Part of 2.1): places 'crosspred' as a field of res. (local results)
                    res.crosspred = crosspred;
                end

                %%% 3. It places the cross-prediction reults in a cell array:
                ACI{i_subject,i_masker}     = col1;
                cfg_ACI{i_subject,i_masker} = cfg_ACI_here; 
                results{i_subject,i_masker} = res;
                % idxlambda(i_subject) = res.idxlambda;

                [cross,outs_cross] = Read_crosspred(fname_cross,res.idxlambda);
                PC_matrix(i_subject,:,i_masker)  = outs_cross.PA_mean_re_chance;
                Dev_matrix(i_subject,:,i_masker) = outs_cross.Dev_mean;
                Dev_matrix_plus_SEM(i_subject,:,i_masker) = outs_cross.Dev_mean+outs_cross.Dev_SEM;
                Dev_matrix_min_SEM(i_subject,:,i_masker)  = outs_cross.Dev_mean-outs_cross.Dev_SEM;
                
            end % end if bInclude
        end % end if i_subject
    end % end i_masker
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if flags.do_fig8
    % Code by Leo, taken from l20220325_crossprediction_newversion.m
 
    XTL = [];
    XT = 1:N_subjects;
    for ii = 1:N_subjects
        txt2use = 'S';
        if ii < 10 
            txt2use = [txt2use '0'];
        end
        txt2use = [txt2use num2str(ii)];
        YTL{ii} = txt2use;
        if mod(ii,2)==0
            XTL{ii} = txt2use;
        else
            if N_maskers > 1
                XTL{ii} = '';
            else
                XTL{ii} = txt2use;
            end
        end
    end
    
    % figure('Position',[100 100 800 350]);
    figure('Position',[100 100 720 415]);
    il_tiledlayout(1,N_maskers,'TileSpacing','tight'); % 'tight');

    idx_sig = [0 0 0]; % idx of significant ACIs
    for i_column = 1:N_maskers
        if N_maskers == 1
            Colour_here = masker_colours{1};
        else
            switch mod(i_column,N_maskers)
                case 1 % white
                    Colour_here = masker_colours{1};
                case 2 % bump
                    Colour_here = masker_colours{2};
                case 0 % MPS
                    Colour_here = masker_colours{3};
            end
        end
        Dev_matrix_here = Dev_matrix(:,:,i_column);
        Dev_matrix_plus = Dev_matrix_plus_SEM(:,:,i_column);
        % Dev_matrix_min  = Dev_matrix_min_SEM(:,:,i_column);
        % Dev_diag(:,i_column) = diag(Dev_matrix_here);
        Dev_diag(:,i_column) = diag(Dev_matrix_plus);
        Dev_matrix_nodiag = Dev_matrix_here-repmat(Dev_diag(:,i_column)',N_subjects,1); % matrix with no diagonal

        PC_matrix_here = PC_matrix(:,:,i_column);
        PC_diag(:,i_column) = diag(PC_matrix_here);
        PC_matrix_nodiag = PC_matrix_here-repmat(PC_diag(:,i_column)',N_subjects,1); % matrix with no diagonal

        [Me_diag(i_column),Me_nodiag(i_column), opts] = il_avg_diag_nodiag(PC_matrix_here);
        
        if flags.do_fig8
            il_nexttile(i_column)
            imagesc(1:N_subjects,1:N_subjects,PC_matrix_here); hold on
            colormap('gray')
            caxis([0 10])
            set(gca,'XTick',XT);
            set(gca,'YTick',XT);
            if i_column ~= 1
                set(gca,'YTickLabel',[]);
            end
            if i_column == N_maskers
                c = colorbar;
            end
            for ii = 1:N_subjects
                offxy = 0.5;
                offxy = offxy*.9;
                if Dev_diag(ii,i_column) < 0 % significance using 'Dev'
                    % idx = find(PC_matrix_nodiag(:,ii)>0); % Criterion using PA
                    idx = find(Dev_matrix_plus(:,ii)<0);
                    if ~isempty(idx)
                        for jj = 1:length(idx)
                            il_plot_square(ii,idx(jj),'m',offxy,'--'); % Magenta colour
                            if ii ~= jj % only counting cross predictions
                                idx_sig(i_column) = idx_sig(i_column)+1;
                            end
                        end
                    end
                else
                    %%% Add arrow:
                    Y_start = 2;
                    [xaf,yaf] = ds2nfu(ii*[1 1],[Y_start Y_start-1]); % Convert to normalized figure units
                    annotation('arrow',xaf,yaf,'color','r');
                    %%%
                end
                % Plots a box:
                il_plot_square(ii,ii,Colour_here,offxy);
            end
        end
    end
    fprintf('Delta PAexp values ranged between %.1f and %.1f\n', ...
        min(min(min(PC_matrix))),max(max(max(PC_matrix))));
    fprintf('ACIs leading to significant predictions using the CVD criterion: %.1f, %.1f, %.1f (white, bump, MPS)\n', ...
        idx_sig);
    FS_here = 10;

    opts_txt = {'Units','Normalized','FontWeight','bold','FontSize',14};
    for i_masker = 1:N_maskers
        nexttile(i_masker);
        ylabel('Data from');
        switch i_masker
            case 1
                set(gca,'YTickLabel',YTL);
                if N_maskers ~= 1
                    text(0,1.05,'A.',opts_txt{:});
                    set(gca,'XTickLabel',[]);
                end
            case 2
                set(gca,'YTickLabel',[]);
                text(0,1.05,'B.',opts_txt{:});
            case 3
                set(gca,'YTickLabel',[]);
                text(0,1.05,'C.',opts_txt{:});
                set(gca,'XTickLabel',[]);
        end
        set(gca,'XTickLabel',XTL);
        set(gca,'FontSize',FS_here);
    end
    if flags.do_fig8
        % c.Label.String = 'Percentage accuracy benefit (%)';
        set(c,'Limits',[0 10.5]);
        fig_name = 'fig8a-crosspred_participant';
        
        txt_title = 'Experimental data';
        if N_maskers == 1
            txt_title = [noise_types_label{1} ': ' txt_title];
        end
        ht = title(txt_title);
        set(ht,'FontSize',FS_here);
    end
    
    for i_masker = 1:N_maskers
        xlabel('ACI from');
    end
    h(end+1) = gcf;
    hname{end+1} = fig_name;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flags.do_fig9
    for j = 1:N_subjects
        for k = 1:N_maskers
            for i = 1:N_maskers 
                idxlambda = results{j,i}.idxlambda;

                PC_test  = results{j,i}.crosspred(k).PC_test;
                yvar = PC_test(idxlambda,:)-PC_test(end,:);
                PC_matrix(j,k,i) = 100*mean(yvar); % mean converted to a percentage

                Dev_test = results{j,i}.crosspred(k).Dev_test; 
                yvar = Dev_test(idxlambda,:)-Dev_test(end,:);
                Dev_matrix(j,k,i) = mean(yvar); % mean 
                Dev_matrix_plus_SEM(j,k,i) = mean(yvar)+1.64*sem(yvar);
                Dev_matrix_min_SEM(j,k,i)  = mean(yvar)-1.64*sem(yvar);
            end
        end
    end
        
    % Build the exact matrix to plot
    for j = 1:N_subjects
        for k = 1:N_maskers
            PC_pool_here{j}(k,:) = transpose(squeeze(PC_matrix(j,k,:)));
            Dev_pool_here{j}(k,:) = transpose(squeeze(Dev_matrix(j,k,:)));
        end
    end
    PC= [];
    PC_rem_tmp = [];
    Dev = [];
    for j = 1:4
        PC = [PC PC_pool_here{j}];
        PC_rem_tmp = [PC_rem_tmp PC_pool_here{j}-repmat(diag(PC_pool_here{j})',3,1)];
        Dev = [Dev Dev_pool_here{j}];
    end
    PC_pool = PC;
    PC_pool_diag_removed = PC_rem_tmp;
    Dev_pool = Dev;

    PC= [];
    PC_rem_tmp = [];
    Dev = [];
    for j = 5:8
        PC = [PC PC_pool_here{j}];
        PC_rem_tmp = [PC_rem_tmp PC_pool_here{j}-repmat(diag(PC_pool_here{j})',3,1)];
        Dev = [Dev Dev_pool_here{j}];
    end
    PC_pool  = [PC_pool; PC];
    PC_pool_diag_removed = [PC_pool_diag_removed; PC_rem_tmp];
    Dev_pool = [Dev_pool; Dev];

    PC= [];
    PC_rem_tmp = [];
    Dev = [];
    for j = 9:12
        PC = [PC PC_pool_here{j}];
        PC_rem_tmp = [PC_rem_tmp PC_pool_here{j}-repmat(diag(PC_pool_here{j})',3,1)];
        Dev = [Dev Dev_pool_here{j}];
    end
    PC_pool = [PC_pool; PC];
    PC_pool_diag_removed = [PC_pool_diag_removed; PC_rem_tmp];
    Dev_pool = [Dev_pool; Dev];
    
    % figure;
    figure('Position',[100 100 1000 420]); % 560 420]);
    il_tiledlayout(1,2,'TileSpacing','compact'); % 'tight');

    il_nexttile(1);
    x_var = 1:size(Dev_pool,2);
    y_var = 1:size(Dev_pool,1);
    imagesc(x_var,y_var,Dev_pool); hold on
    colormap('gray')
    caxis([-10 0])

    c = colorbar;
    c.Label.String = 'Deviance / trial benefit (adim)';

    for ii = 1:4
        text((ii-1)*3+.7,0.7,Subjects{ii},'FontSize',10,'Color','b');
    end
    for ii = 5:8
        text(mod((ii-1),4)*3+.7,4-.3,Subjects{ii},'FontSize',10,'Color','b');
    end
    for ii = 9:12
        text(mod((ii-1),4)*3+.7,7-.3,Subjects{ii},'FontSize',10,'Color','b');
    end
    %%% Drawing a separator between participants:
    plot(3.5*[1 1],[0 9.5],'Color','b');
    plot(3.5*[1 1],[0 9.5],'Color','b');

    plot(6.5*[1 1],[0 9.5],'Color','b');
    plot(6.5*[1 1],[0 9.5],'Color','b');

    plot(9.5*[1 1],[0 9.5],'Color','b');
    plot(9.5*[1 1],[0 9.5],'Color','b');

    plot([0 12.5],3.5*[1 1],'Color','b');
    plot([0 12.5],6.5*[1 1],'Color','b');
    plot([0 12.5],9.5*[1 1],'Color','b');
    set(gca,'FontSize',10);
    %%%

    disp('')
    set(gca,'XTick',x_var)
    set(gca,'YTick',y_var)

    XTL = repmat({'W','BP','MPS'},1,4);
    YTL = repmat({'W','BP','MPS'},1,3);
    set(gca,'XTickLabel',XTL);
    set(gca,'YTickLabel',YTL);

    xlabel('ACI from condition');
    ylabel('Data from condition');

    text(0,1.05,'D.','Units','Normalized','FontSize',14,'FontWeight','Bold');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % figure;
    il_nexttile(2);
    x_var = 1:size(PC_pool,2);
    y_var = 1:size(PC_pool,1);
    imagesc(x_var,y_var,PC_pool); hold on;
    colormap('gray')
    caxis([0 10])

    %%% Drawing a separator between participants:
    plot(3.5*[1 1],[0 9.5],'Color','b');
    plot(3.5*[1 1],[0 9.5],'Color','b');

    plot(6.5*[1 1],[0 9.5],'Color','b');
    plot(6.5*[1 1],[0 9.5],'Color','b');

    plot(9.5*[1 1],[0 9.5],'Color','b');
    plot(9.5*[1 1],[0 9.5],'Color','b');

    plot([0 12.5],3.5*[1 1],'Color','b');
    plot([0 12.5],6.5*[1 1],'Color','b');
    plot([0 12.5],9.5*[1 1],'Color','b');
    set(gca,'FontSize',10);
    %%%

    PC_pool_diag_removed = transpose(PC_pool_diag_removed);
    offxy = 0.5;
    offxy = offxy*.9;
    for jj = 1:size(PC_pool_diag_removed,1)
        idx_here = find(PC_pool_diag_removed(jj,:)>0);
        for ii = 1:length(idx_here)
            il_plot_square(jj,idx_here(ii),'m',offxy,'--');
        end
    end

    c = colorbar;
    c.Label.String = 'Percentage accuracy benefit (%)';

    for ii = 1:4
        text((ii-1)*3+.7,0.7,Subjects{ii},'FontSize',10,'Color','b');
    end
    for ii = 5:8
        text(mod((ii-1),4)*3+.7,4-.3,Subjects{ii},'FontSize',10,'Color','b');
    end
    for ii = 9:12
        text(mod((ii-1),4)*3+.7,7-.3,Subjects{ii},'FontSize',10,'Color','b');
    end

    disp('')
    set(gca,'XTick',x_var)
    set(gca,'YTick',y_var)

    XTL = repmat({'W','BP','MPS'},1,4);
    YTL = repmat({'W','BP','MPS'},1,3);
    set(gca,'XTickLabel',XTL);
    set(gca,'YTickLabel',YTL);

    xlabel('ACI from condition');
    ylabel('Data from condition');

    text(0,1.05,'A.','Units','Normalized','FontSize',14,'FontWeight','Bold');

    h(end+1) = gcf;
    hname{end+1} = 'fig9a-crosspred_masker';

    % c.Label.String = 'Deviance / trial benefit (adim)';

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % figure;
    figure('Position',[100 100 560 250]); % 560 420]);
    % il_tiledlayout(1,2,'TileSpacing','compact'); % 'tight');

    PC_all = (PC_pool(1:3,:)+PC_pool(4:6,:)+PC_pool(7:9,:))/3;
    PC_all = (PC_all(:,1:3)+PC_all(:,4:6)+PC_all(:,7:9)+PC_all(:,10:12))/4;

    x_var = 1:size(PC_all,2);
    y_var = 1:size(PC_all,1);
    imagesc(x_var,y_var,PC_all); hold on
    colormap('gray')
    caxis([0 10])
    set(gca,'FontSize',10);

    c = colorbar;
    c.Label.String = 'Percentage accuracy benefit (%)';

    for ii = 1:3
        for jj = 1:3
            text(ii-0.1,jj,sprintf('%.1f',PC_all(ii,jj)),'FontSize',10,'Color','b');
        end
    end
    % %%%
    set(gca,'XTick',x_var)
    set(gca,'YTick',y_var)

    XTL = repmat({'W','BP','MPS'},1,4);
    YTL = repmat({'W','BP','MPS'},1,3);
    set(gca,'XTickLabel',XTL);
    set(gca,'YTickLabel',YTL);

    xlabel('ACI from condition');
    ylabel('Data from condition');

    text(0,1.05,'B.','Units','Normalized','FontSize',14,'FontWeight','Bold');
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    xlabel('ACI from condition');
    ylabel('Data from condition');

    % text(0,1.05,'A.','Units','Normalized','FontSize',14,'FontWeight','Bold');

    h(end+1) = gcf;
    hname{end+1} = 'fig9b-crosspred_masker';

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figs 10 and 11
if flags.do_fig8b || flags.do_fig12
    
    % All simulation results were stored in one output folder
    if isunix
        dir_results = dir_ACI_sim;
    else
        error('Leo put your folder here...')
        dir_results = [];
    end
    Subjects_ID = {'Results-S01-v1','Results-S02-v1','Results-S03-v1','Results-S04-v1', ...
                   'Results-S05-v1','Results-S06-v1','Results-S07-v1','Results-S08-v1', ...
                   'Results-S09-v1','Results-S10-v1','Results-S11-v1','Results-S12-v1'};
    model_str = 'osses2022a';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flags.do_fig8b || flags.do_fig12
    N_noises = length(noise_str);
    N_subjects = length(Subjects_ID);
    
    for i_noise = 1:N_noises
        PA_mean_chance = []; % emptied after having plotted every noise
        for i_subj = 1:N_subjects
            if i_noise == 1
                Subj_ID = strsplit(Subjects_ID{i_subj},'-');
                labels2use{i_subj} = Subj_ID{2};
            end
        
            dir_subj = [dir_results Subjects_ID{i_subj} filesep];
        
            ACI_fname_folder = [dir_subj 'ACI-' model_str '-speechACI_Logatome-' noise_str{i_noise} '-nob-gt-l1glm+pyrga-rev4' filesep];
            Crossfile = [ACI_fname_folder 'Crosspred.mat'];
            
            cfg_ACI_here = load([ACI_fname_folder(1:end-1) '.mat']);
            res_here = cfg_ACI_here.results;
            %%% Save some data for group-level analysis
            idxlambda = res_here.idxlambda;
            
            % [cross,outs_cross] = Read_crosspred(Crossfile, idxlambda);
            % PA_mean_chance(i_subj,:) = outs_cross.PA_mean_re_chance;
            [cross,outs_cross] = Read_crosspred(Crossfile, idxlambda);
            PC_matrix(i_subj,:,i_noise)  = outs_cross.PA_mean_re_chance;
            Dev_matrix(i_subj,:,i_noise) = outs_cross.Dev_mean;
            Dev_matrix_plus_SEM(i_subj,:,i_noise) = outs_cross.Dev_mean+outs_cross.Dev_SEM;
            Dev_matrix_min_SEM(i_subj,:,i_noise)  = outs_cross.Dev_mean-outs_cross.Dev_SEM;
                
            % num1 = res_here.FitInfo.Dev_test./res_here.FitInfo.CV.TestSize;
            % num2 = num1(end,:); % referential null ACI 
            % G_Devtrial_opt(:,i_subject,i_masker) = num1(idxlambda,:) - num2;
             
            num1 = res_here.FitInfo.PC_test;
            num2 = num1(end,:); % referential null ACI 
            G_PA_sim(:,i_subj,i_noise) = 100*(num1(idxlambda,:) - num2);
            
            %%%
            res = [];
            res.idxlambda = idxlambda;
            res.crosspred = cross;
            results{i_subj,i_noise} = res;
            %%%
            
            if flags.do_fig12
                % i_subj: ACI from
                % k: Data from
                for k = 1:N_subjects
                    PC_test  = cross(k).PC_test;
                    yvar = PC_test(idxlambda,:)-PC_test(end,:);
                    PC_matrix(i_subj,k,i_noise) = 100*mean(yvar);
                end
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flags.do_fig8b
    %%% 4. Plotting the results:
	% 3 matrices of cross-prediction accuracy
    % xvar = results{1,1}.crosspred(1).lambdas;
    XTL = [];
    XT = 1:N_subjects;
    for ii = 1:N_subjects
        txt2use = 'S';
        if ii < 10 
            txt2use = [txt2use '0'];
        end
        txt2use = [txt2use num2str(ii)];
        YTL{ii} = txt2use;
        if mod(ii,2)==0
            XTL{ii} = txt2use;
        else
            if N_maskers > 1
                XTL{ii} = '';
            else
                XTL{ii} = txt2use;
            end
        end
    end

    figure('Position',[100 100 720 415]);
    il_tiledlayout(1,N_maskers,'TileSpacing','tight');

    for i_column = 1:N_maskers
        if N_maskers == 1
            Colour_here = masker_colours{1};
        else
            switch mod(i_column,N_maskers)
                case 1 % white
                    Colour_here = masker_colours{1};
                case 2 % bump
                    Colour_here = masker_colours{2};
                case 0 % MPS
                    Colour_here = masker_colours{3};
            end
        end
        Dev_matrix_here = Dev_matrix_plus_SEM(:,:,i_column);
        Dev_diag(:,i_column) = diag(Dev_matrix_here);

        il_nexttile(i_column)
        
        PC_matrix_here = PC_matrix(:,:,i_column);
        [Me_diag(i_column),Me_nodiag(i_column), opts] = il_avg_diag_nodiag(PC_matrix_here);
        
        if flags.do_fig8b
            imagesc(1:N_subjects,1:N_subjects,PC_matrix_here); hold on
            colormap('gray')
            caxis([14 24.5])
            set(gca,'XTick',XT);
            set(gca,'YTick',XT);
            idx_total = 0;
            for ii = 1:N_subjects
                offxy = 0.5;
                offxy = offxy*.9;
                if Dev_diag(ii,i_column) < 0 % significance using 'Dev'
                    % idx = find(PC_matrix_nodiag(:,ii)>0); % Criterion using PA
                    idx = find(Dev_matrix_plus_SEM(:,ii,i_column)<0);
                    if ~isempty(idx)
                        for jj = 1:length(idx)
                            il_plot_square(ii,idx(jj),'m',offxy,'--'); % Magenta colour
                        end
                        idx_total = idx_total + length(idx);
                    end
                else
                    %%% Add arrow:
                    Y_start = 2;
                    [xaf,yaf] = ds2nfu(ii*[1 1],[Y_start Y_start-1]); % Convert to normalized figure units
                    annotation('arrow',xaf,yaf,'color','r');
                    %%%
                end
                % Plots a box:
                il_plot_square(ii,ii,Colour_here,offxy);
            end
        end

        fprintf('PA - Noise cond=%s: idx=%.0f of %.0f with significant predictions\n',noise_types_label{i_noise},idx_total,N_subjects*N_subjects);
    end

    FS_here = 10;

    opts_txt = {'Units','Normalized','FontWeight','bold','FontSize',14};
    for i_masker = 1:N_maskers
        il_nexttile(i_masker);
        if i_masker == 1
            ylabel('Simulation data from');
            set(gca,'YTickLabel',YTL);
            if N_maskers > 1
                text(0,1.05,'D.',opts_txt{:});
            end
        else
            set(gca,'YTickLabel',[]);
            switch i_masker
                case 2
                    text(0,1.05,'E.',opts_txt{:});
                case 3
                    text(0,1.05,'F.',opts_txt{:});
            end
        end
    
    end
    set(gca,'FontSize',FS_here);
    set(gca,'XTickLabel',XTL);
    
    c = colorbar;

    for i_masker = 1:N_maskers
        nexttile(i_masker);
        ht = title([noise_types_label{i_masker} ': Simulated data']);
        set(ht,'FontSize',FS_here);
        xlabel('ACI_sim from','interpreter','none');
    end
    h(end+1) = gcf;
    if flags.do_fig8b
        hname{end+1} = 'fig8b-crosspred_model_participant'; 
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flags.do_fig12
    between_noise = 1;
    if between_noise == 1
        % Build the exact matrix to plot
        for j = 1:N_subjects
            for k = 1:N_maskers
                PC_pool_here{j}(k,:) = transpose(squeeze(PC_matrix(j,k,:)));
                % Dev_pool_here{j}(k,:) = transpose(squeeze(Dev_matrix(j,k,:)));
            end
        end
        PC= [];
        PC_rem_tmp = [];
        % Dev = [];
        for j = 1:4
            PC = [PC PC_pool_here{j}];
            PC_rem_tmp = [PC_rem_tmp PC_pool_here{j}-repmat(diag(PC_pool_here{j})',3,1)];
            % Dev = [Dev Dev_pool_here{j}];
        end
        PC_pool = PC;
        PC_pool_diag_removed = PC_rem_tmp;
        % Dev_pool = Dev;
         
        PC= [];
        PC_rem_tmp = [];
        % Dev = [];
        for j = 5:8
            PC = [PC PC_pool_here{j}];
            PC_rem_tmp = [PC_rem_tmp PC_pool_here{j}-repmat(diag(PC_pool_here{j})',3,1)];
            % Dev = [Dev Dev_pool_here{j}];
        end
        PC_pool  = [PC_pool; PC];
        PC_pool_diag_removed = [PC_pool_diag_removed; PC_rem_tmp];
        % Dev_pool = [Dev_pool; Dev];
        
        PC= [];
        PC_rem_tmp = [];
        % Dev = [];
        for j = 9:12
            PC = [PC PC_pool_here{j}];
            PC_rem_tmp = [PC_rem_tmp PC_pool_here{j}-repmat(diag(PC_pool_here{j})',3,1)];
            % Dev = [Dev Dev_pool_here{j}];
        end
        PC_pool = [PC_pool; PC];
        PC_pool_diag_removed = [PC_pool_diag_removed; PC_rem_tmp];
        % Dev_pool = [Dev_pool; Dev];
    end
      
    if between_noise == 0
        XTL = [];
        XT = 1:N_subjects;
        for ii = 1:N_subjects
            txt2use = 'S';
            if ii < 10 
                txt2use = [txt2use '0'];
            end
            txt2use = [txt2use num2str(ii)];
            YTL{ii} = txt2use;
            if mod(ii,2)==0
                XTL{ii} = txt2use;
            else
                XTL{ii} = '';
            end
        end
   
        figure('Position',[100 100 800 600]);
        il_tiledlayout(2,3,'TileSpacing','tight');
         
        for i_column = 1:3
            switch mod(i_column,3)
                case 1 % white
                    Colour_here = masker_colours{1};
                case 2 % bump
                    Colour_here = masker_colours{2};
                case 0 % MPS
                    Colour_here = masker_colours{3};
            end
            il_nexttile(i_column)
            PC_matrix_here = PC_matrix(:,:,i_column);

            il_nexttile(i_column)
            imagesc(1:N_subjects,1:N_subjects,PC_matrix_here); hold on
            colormap('gray')
            caxis([14 24.5])
            set(gca,'XTick',XT);
            set(gca,'YTick',XT);
            
            for ii = 1:N_subjects
                
                offxy = 0.5;
                offxy = offxy*.9;
                il_plot_square(ii,ii,Colour_here,offxy);
                % if Dev_diag(ii,i_column) ~= 0
                %     idx = find(Dev_matrix_nodiag(:,ii)<0);
                %     if ~isempty(idx)
                %         for jj = 1:length(idx)
                %             il_plot_square(ii,idx(jj),'m',offxy,'--'); % Magenta colour
                %         end
                %     end
                % end
            end
        end
        
        FS_here = 10;
        
        il_nexttile(1);
        set(gca,'XTickLabel',[]);
        ylabel('Data from');
        set(gca,'YTickLabel',YTL);
        opts_txt = {'Units','Normalized','FontWeight','bold','FontSize',14};
        text(0,1.05,'A.',opts_txt{:});
        set(gca,'FontSize',FS_here);
        
        il_nexttile(2);
        set(gca,'XTickLabel',[]);
        set(gca,'YTickLabel',[]);
        text(0,1.05,'B.',opts_txt{:});
        set(gca,'FontSize',FS_here);
        
        il_nexttile(3);
        text(0,1.05,'C.',opts_txt{:});
        c = colorbar;
        c.Label.String = 'Deviance / trial benefit (adim)';
        set(gca,'XTickLabel',[]);
        set(gca,'YTickLabel',[]);
        % set(c,'Limits',[-10.5 0])
        set(gca,'FontSize',FS_here);
        
        nexttile(4);
        ylabel('Simulation data from');
        set(gca,'YTickLabel',YTL);
        set(gca,'XTickLabel',XTL);
        set(gca,'FontSize',FS_here);
        
        il_nexttile(5);
        set(gca,'YTickLabel',[]);
        set(gca,'XTickLabel',XTL);
        set(gca,'FontSize',FS_here);
        
        il_nexttile(6);
        c = colorbar;
        c.Label.String = 'Percentage accuracy benefit (%)';
        % set(c,'Limits',[0 10.5])
        
        set(gca,'YTickLabel',[]);
        set(gca,'XTickLabel',XTL);
        set(gca,'FontSize',FS_here);
        
        for i_masker = 1:N_maskers
            nexttile(i_masker);
            ht = title(noise_types_label{i_masker});
            set(ht,'FontSize',FS_here);
            
            nexttile(i_masker+3);
            xlabel('ACI_s_i_m from');
        end
        h(end+1) = gcf;
        hname{end+1} = 'fig11-crosspred_model_participant';
        
        %%%
        PC_rem = [];
        for idx_noise = 1:3
            PC_rem{idx_noise} = PC_matrix(:,:,idx_noise)-repmat(diag(PC_matrix(:,:,idx_noise))',N_subjects,1);
            idx = find(PC_rem{idx_noise}>0); 
            length(idx)
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Replicating the figure, but showing averages of cross predictions (between subjects)
        figure('Position',[100 100 700*2-200 550]);
        il_tiledlayout(2,2,'TileSpacing','tight');

        XL = [0.5 3.5]; 

        il_nexttile(1);
        plot(XL,0*[1 1],'k--'); hold on; grid on

        il_nexttile(2);
        plot(XL,0*[1 1],'k--'); hold on; grid on

        il_nexttile(3);
        plot(XL,0*[1 1],'k--'); hold on; grid on

        il_nexttile(4);
        plot(XL,0*[1 1],'k--'); hold on; grid on

        for i_masker = 1:N_maskers
            for i_subject = 1:N_subjects

                offx = 0.05*(i_subject - ((N_subjects+2)/2-.5)); %  0.05*i_subject;

                idx_here = 1:N_subjects;
                idx_here(i_subject) = []; % removing the diagonal
                y2plot  = G_PA_sim(:,i_subject,i_masker); % 10 folds x 12 subjects x 
                y2plot2 = PC_rem{i_masker}(idx_here,i_subject);

                MCol = 'w';
                il_nexttile(1);
                errorbar(i_masker+offx, mean(y2plot), 1.64*sem(y2plot),'d', 'Color', masker_colours{i_masker},'MarkerFaceColor',MCol); hold on

                il_nexttile(3);
                Me = median(y2plot2);
                errU = max(y2plot2)-Me;
                errL = Me-min(y2plot2);
                errorbar(i_masker+offx, Me, errL, errU, 'Color', masker_colours{i_masker},'MarkerFaceColor',MCol); hold on

                Benefit_re_diag_MI(i_subject,i_masker) = min(y2plot2);
                Benefit_re_diag_MA(i_subject,i_masker) = max(y2plot2);
                
                %%% Incorrect trials:
                y2plot  = nan([N_subjects 1]); 
                y2plot2 = nan([N_subjects 1]); 

                il_nexttile(2);
                Me = mean(y2plot);
                plot(i_masker+offx, Me, 'd', 'Color', masker_colours{i_masker},'MarkerFaceColor',MCol); hold on

                il_nexttile(4);
                Me = mean(y2plot2);
                plot(i_masker+offx,Me, 'd', 'Color', masker_colours{i_masker},'MarkerFaceColor',MCol); hold on            
            end
            y2plot  = mean(G_PA_sim(:,:,i_masker));
            il_nexttile(1);
            offx = 0.05*(N_subjects+2 - (N_subjects+2)/2);
            errorbar(i_masker+offx, mean(y2plot) , 1.64*sem(y2plot) , 'o', 'Color', masker_colours{i_masker},'LineWidth',2,'MarkerFaceColor',masker_colours{i_masker}); hold on
        end

        XLabel = 'Masker';
        % subplot(2,1,1)
        il_nexttile(1);
        title('All trials')
        grid on
        xlim(XL);
        % xlabel(XLabel); 
        ylim([4 26])
        set(gca, 'XTick', 1:N_maskers); 
        set(gca, 'XTickLabels', []); 
        ylabel('PA benefit (%)')
        xcoor = 0.01;
        ycoor = 0.94;
        txt_opts = {'Units','Normalized','FontWeight','Bold','FontSize',14};
        text(xcoor,ycoor,'D',txt_opts{:});
        set(gca,'FontSize',FS_here);
        
        il_nexttile(2);
        title('Incorrect trials only')
        grid on
        xlim(XL);
        set(gca, 'XTick', 1:N_maskers); 
        set(gca, 'XTickLabels', []); 
        set(gca, 'YTickLabels',[]);
        text(xcoor,ycoor,'NaN',txt_opts{:});
        set(gca,'FontSize',FS_here);

        il_nexttile(3);
        grid on
        xlim(XL)
        ylim([-10 5])
        xlabel(XLabel); 
        set(gca, 'XTick', 1:N_maskers); 
        set(gca, 'XTickLabels', noise_types_label); 
        ylabel('PA benefit re. auto prediction (%)')
        text(xcoor,ycoor,'E',txt_opts{:});
        set(gca,'FontSize',FS_here);
        
        il_nexttile(4);

        xlim(XL)
        ylim([-10 5])
        xlabel(XLabel); 
        set(gca, 'XTick', 1:N_maskers); 
        set(gca, 'XTickLabels', noise_types_label); 
        set(gca, 'YTickLabels',[]);
        % ylabel('Percent accuracy benefit (%)')
        text(xcoor,ycoor,'NaN',txt_opts{:});
        set(gca,'FontSize',FS_here);
        
        h(end+1) = gcf;
        hname{end+1} = ['fig11b-metrics-benefit-model-cross'];
        
        % fprintf('Delta PAsim values ranged between %.2f and %.2f\n', ...
        %     min(min(min(G_PA_sim))),max(max(max(G_PA_sim))));
        fprintf('Delta PAsim values ranged between %.1f and %.1f\n', ...
            min(min(PA_mean_chance)),max(max(PA_mean_chance)));
        fprintf('Benefit re. main diagonal between %.1f and %.1f (median=%.1f)\n', ...
            min(min(Benefit_re_diag_MI)), ...
            max(max(Benefit_re_diag_MA)), ...
            prctile([Benefit_re_diag_MI(:) Benefit_re_diag_MA(:)],50));
        
    else
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % figure;
        % figure('Position',[100 100 1000 420]); 
        figure('Position',[100 100 560 420]); 
        il_tiledlayout(1,1,'TileSpacing','compact'); % 'tight');
        
        % figure;
        il_nexttile(1);
        x_var = 1:size(PC_pool,2);
        y_var = 1:size(PC_pool,1);
        imagesc(x_var,y_var,PC_pool); hold on;
    	colormap('gray')
        % caxis([0 10])
        
        %%% Drawing a separator between participants:
        plot(3.5*[1 1],[0 9.5],'Color','b');
        plot(3.5*[1 1],[0 9.5],'Color','b');
        
        plot(6.5*[1 1],[0 9.5],'Color','b');
        plot(6.5*[1 1],[0 9.5],'Color','b');
        
        plot(9.5*[1 1],[0 9.5],'Color','b');
        plot(9.5*[1 1],[0 9.5],'Color','b');
        
        plot([0 12.5],3.5*[1 1],'Color','b');
        plot([0 12.5],6.5*[1 1],'Color','b');
        plot([0 12.5],9.5*[1 1],'Color','b');
        set(gca,'FontSize',10);
        %%%
        
        PC_pool_diag_removed = transpose(PC_pool_diag_removed);
        offxy = 0.5;
        offxy = offxy*.9;
        for jj = 1:size(PC_pool_diag_removed,1)
            idx_here = find(PC_pool_diag_removed(jj,:)>0);
            for ii = 1:length(idx_here)
                il_plot_square(jj,idx_here(ii),'m',offxy,'--');
            end
        end
        
        c = colorbar;
        c.Label.String = 'Percentage accuracy benefit (%)';
        
        for ii = 1:4
            text((ii-1)*3+.7,0.7,Subjects{ii},'FontSize',10,'Color','b');
        end
        for ii = 5:8
            text(mod((ii-1),4)*3+.7,4-.3,Subjects{ii},'FontSize',10,'Color','b');
        end
        for ii = 9:12
            text(mod((ii-1),4)*3+.7,7-.3,Subjects{ii},'FontSize',10,'Color','b');
        end
        
        disp('')
        set(gca,'XTick',x_var)
        set(gca,'YTick',y_var)
        
        XTL = repmat({'W','BP','MPS'},1,4);
        YTL = repmat({'W','BP','MPS'},1,3);
        set(gca,'XTickLabel',XTL);
        set(gca,'YTickLabel',YTL);
        
        xlabel('ACI_s_i_m from condition');
        ylabel('Data from condition');
        
        text(0,1.05,'A.','Units','Normalized','FontSize',14,'FontWeight','Bold');
        
        h(end+1) = gcf;
        hname{end+1} = 'fig12a-crosspred_masker-sim';
        
        % c.Label.String = 'Deviance / trial benefit (adim)';
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % figure;
        figure('Position',[100 100 560 250]); % 560 420]);
        % il_tiledlayout(1,2,'TileSpacing','compact'); % 'tight');
         
        PC_all = (PC_pool(1:3,:)+PC_pool(4:6,:)+PC_pool(7:9,:))/3;
        PC_all = (PC_all(:,1:3)+PC_all(:,4:6)+PC_all(:,7:9)+PC_all(:,10:12))/4;
        
        x_var = 1:size(PC_all,2);
        y_var = 1:size(PC_all,1);
        imagesc(x_var,y_var,PC_all); hold on
        colormap('gray')
        caxis([14.5 23.5])
        set(gca,'FontSize',10);
        
        c = colorbar;
        c.Label.String = 'PA benefit (%)';
         
        for ii = 1:3
            for jj = 1:3
                text(ii-0.1,jj,sprintf('%.1f',PC_all(ii,jj)),'FontSize',10,'Color','b');
            end
        end
        % %%%
        set(gca,'XTick',x_var)
        set(gca,'YTick',y_var)
         
        XTL = repmat({'W','BP','MPS'},1,4);
        YTL = repmat({'W','BP','MPS'},1,3);
        set(gca,'XTickLabel',XTL);
        set(gca,'YTickLabel',YTL);
        
        xlabel('ACI_s_i_m from condition');
        ylabel('Data from condition');
         
        text(0,1.05,'B.','Units','Normalized','FontSize',14,'FontWeight','Bold');
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        xlabel('ACI_s_i_m from condition');
        ylabel('Data from condition');
         
        % text(0,1.05,'A.','Units','Normalized','FontSize',14,'FontWeight','Bold');

        h(end+1) = gcf;
        hname{end+1} = 'fig12b-crosspred_masker-sim';
        
        fprintf('Between-noise: Delta PAsim values ranged between %.1f and %.1f\n', ...
            min(min(PC_pool)),max(max(PC_pool)));
        % fprintf('Benefit re. main diagonal between %.1f and %.1f (median=%.1f)\n', ...
        %     min(min(Benefit_re_diag_MI)), ...
        %     max(max(Benefit_re_diag_MA)), ...
        %     prctile([Benefit_re_diag_MI(:) Benefit_re_diag_MA(:)],50));
        
        idx_here = find(PC_pool_diag_removed(:,:)>0); 
        fprintf('Number of cases where PA ''between'' was higher than the autoprediction: %.0f (out of %.0f)\n', ...
            length(idx_here),length(PC_pool_diag_removed(:))-3*N_subjects);
        
        idx_here_above2   = find(PC_pool_diag_removed(:,:)>=2); 
        idx_here_above1p5 = find(PC_pool_diag_removed(:,:)>=1.5); 
        
        idx_here_below1   = find(PC_pool_diag_removed(idx_here)<1);
        
        prctile(PC_pool_diag_removed(idx_here),90)
        prctile(PC_pool_diag_removed(idx_here),75)
        prctile(PC_pool_diag_removed(idx_here),50)
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data.h = h;
data.hname = hname;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function il_tiledlayout(N,M,TileSpacing,TileSpacing_option)

if nargin < 4
    TileSpacing_option = 'Compact';
end
bExist = exist('tiledlayout','file'); % tiledlayout.p
if bExist
    tiledlayout(N,M,TileSpacing,TileSpacing_option);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function il_nexttile(N)

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [PC_fold, Dev_fold] = il_tbtpred_inc(cfg_ACI,results,data_passation)
% Adapted version by Alejandro

N_lambdas = length(results.lambdas);
N_folds = cfg_ACI.N_folds;

N = length(results.FitInfo.PC_test_t);
PC_all = nan([N_lambdas, N_folds, N]);
Dev_all = nan(size(PC_all));

if nargin >= 3
    idxs_inc  = find(data_passation.is_correct(cfg_ACI.idx_analysis)==0);
    idxs_corr = find(data_passation.is_correct(cfg_ACI.idx_analysis)==1); % correct indexes
else
    fprintf('All trials will be returned...')
end
                
for i_lambda = 1:N_lambdas
    for i_fold = 1:N_folds
        PC_test_t  = squeeze(results.FitInfo.PC_test_t(i_lambda,i_fold,:));
        Dev_test_t = squeeze(results.FitInfo.Dev_test_t(i_lambda,i_fold,:));
        CVtest = find(results.FitInfo.CV.test(i_fold));
        
        if nargin >= 3
            [idx_corr_this_fold,idx_corr_sort] = intersect(CVtest,idxs_corr);
        end
        % idx_inc_this_fold  = intersect(CVtest,idxs_inc);
        N_here = length(CVtest);
        
        if length(CVtest)~=length(Dev_test_t)
            %fprintf(['unequal CVtest and MSEtest_t at fold # ' num2str(i_fold) '\n'])
            if (length(CVtest)==length(Dev_test_t)-1) && (Dev_test_t(end)==0)
                PC_test_t = PC_test_t(1:end-1);
                Dev_test_t = Dev_test_t(1:end-1);
            else
                error('Problem with the length of vectors Dev_test_t and PC_test_t')
            end
        end
        
        PC_all(i_lambda,i_fold,1:N_here)  = PC_test_t;
        Dev_all(i_lambda,i_fold,1:N_here) = Dev_test_t;
        
        if nargin >= 3
            % setting back to NaN those trials within the fold that were correct:
            PC_all(i_lambda,i_fold,idx_corr_sort)  = nan; 
            Dev_all(i_lambda,i_fold,idx_corr_sort) = nan;
        end
    end
end

PC_fold = nanmean(PC_all,3); % across third dimension
Dev_fold = nanmean(Dev_all,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function il_plot_square(xcoor,ycoor,Colour_here,offxy,LineStyle)

if nargin < 5
    LineStyle = '-';
end
plot_opts = {'Color',Colour_here,'LineStyle',LineStyle};
% Plots a box:
plot([xcoor-offxy xcoor+offxy],(ycoor-offxy)*[1 1]      ,plot_opts{:}); % bottom side
plot((xcoor-offxy)*[1 1]      ,[ycoor-offxy ycoor+offxy],plot_opts{:}); % left side
plot((xcoor+offxy)*[1 1]      ,[ycoor-offxy ycoor+offxy],plot_opts{:}); % right side
plot([xcoor-offxy xcoor+offxy],(ycoor+offxy)*[1 1]      ,plot_opts{:}); % top side

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function il_plot_square_XY(xcoor,ycoor,Colour_here,offx,offy,LineStyle,LW)

if nargin < 5
    LineStyle = '-';
end
plot_opts = {'Color',Colour_here,'LineStyle',LineStyle,'LineWidth',LW};
% Plots a box:
plot([xcoor-offx xcoor+offx],(ycoor-offy)*[1 1]      ,plot_opts{:}); % bottom side
plot((xcoor-offx)*[1 1]      ,[ycoor-offy ycoor+offy],plot_opts{:}); % left side
plot((xcoor+offx)*[1 1]      ,[ycoor-offy ycoor+offy],plot_opts{:}); % right side
plot([xcoor-offx xcoor+offx],(ycoor+offy)*[1 1]      ,plot_opts{:}); % top side

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Me_diag,Me_nodiag,opts] = il_avg_diag_nodiag(insig)

Me_diag = mean(diag(insig));
opts.Me_diag_sem = 1.64*sem(diag(insig));

for i = 1:size(insig,1)
    insig(i,i) = NaN;
end
Me_nodiag = nanmean(insig(:));
insig = insig( find(~isnan(insig)) );
opts.Me_nodiag_sem = 1.64*sem(insig(:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function outs = il_get_SNR(file_savegame)

% data_passation = [];
% cfg_game = [];
load(file_savegame,'data_passation');

n_responses = data_passation.n_responses;
N_trials    = length(n_responses); % completed number of trials
% n_signal    = data_passation.n_targets(1:N_trials);
SNR         = data_passation.expvar;

resume_trial = data_passation.resume_trial; resume_trial(1)=1;resume_trial(end+1)=data_passation.i_current;
N_sessions  = length(resume_trial)-1;

% Fig. 2: Reversal analysis ---------------------------------------

if N_sessions > 10
    % In this case, the participants might have taken one or more pauses...
    idx_odd = find(mod(resume_trial(1:end-1)-1,400)~=0);
    fprintf('Participant %s has %.0f shorter sessions\n',subject,length(idx_odd));
    fprintf('\tI am going to merge that session with the previous one...\n',subject,length(idx_odd));
    for i_odd = 1:length(idx_odd)
        fprintf('\t(session %.0f started at trial %.0f)\n',idx_odd(i_odd),resume_trial(idx_odd));
    end
    resume_trial(idx_odd) = [];
    N_sessions = length(resume_trial)-1;
end
            
idxs2use_all = [];
for i_session = 1:N_sessions

    %%% Leo's way: Analysis only using the reversals: -------------
    % sessions are redefined as blocks of 400 trials
    idxi = 1+400*(i_session-1);
    idxf = 400*(i_session);
    [rev, idx] = Get_mAFC_reversals(SNR(idxi:idxf));
                
    % % Then uses median:
    % medianSNR_old(i_subject, i_masker, i_session) = median(rev);
    % percSNR_L_old(i_subject, i_masker, i_session) = prctile(rev,5);
    % percSNR_U_old(i_subject, i_masker, i_session) = prctile(rev,95);
    iscorr = data_passation.is_correct(idx);
    % perc_corr_old(i_subject, i_masker, i_session) = sum(iscorr)/length(iscorr);

    %%% Alejandro's way:
    idxi = resume_trial(i_session);
    idxf = resume_trial(i_session+1)-1;

    idxi_rev = idx(4); % reversal number 4 (including it) and later
    idxi100 = idxf-100; % reversal number 4 (including it) and later

    idxs2use = idxi+idxi_rev:idxf;
    idxs2use_all = [idxs2use_all idxs2use]; % collating the idxs of the kept trials
    SNR_here = SNR(idxs2use);
    SNR_me(i_session)    = median(SNR_here);
    SNR_percL(i_session) = prctile(SNR_here,5);
    SNR_percU(i_session) = prctile(SNR_here,95);

    % iscorr = data_passation.is_correct(idxs2use);
    % perc_corr(i_subject, i_masker, i_session) = sum(iscorr)/length(iscorr);
    % 
    % iscorr = data_passation.is_correct(idxi:idxf);
    % perc_corr_all(i_subject, i_masker, i_session) = sum(iscorr)/length(iscorr);
    % 
    % iscorr = data_passation.is_correct(idxi_rev:idxf);
    % perc_corr_rev(i_subject, i_masker, i_session) = sum(iscorr)/length(iscorr);
    % 
    % iscorr = data_passation.is_correct(idxi100:idxf);
    % perc_corr_100(i_subject, i_masker, i_session) = sum(iscorr)/length(iscorr);
    % disp('')
end
outs = [];
outs.SNR_me = SNR_me;
outs.SNR_percL = SNR_percL;
outs.SNR_percU = SNR_percU;

outs.idxs2use_all = idxs2use_all;