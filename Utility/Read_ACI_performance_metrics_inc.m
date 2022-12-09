function outs = Read_ACI_performance_metrics_inc(res, data_passation, cfg_ACI)
% function outs = Read_ACI_performance_metrics_inc(res, data_passation, cfg_ACI)
%
% Always calculated on test trials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
idxlambda = res.idxlambda;
 
% switch type
%     case 'test'
%         Dev_raw = res.FitInfo.Dev_test;
%         CV_size = res.FitInfo.CV.TestSize;
%         PC_raw = res.FitInfo.PC_test;
%     case 'train'
%         Dev_raw = res.FitInfo.Dev_train;
%         CV_size = res.FitInfo.CV.TrainSize;
%         PC_raw = res.FitInfo.PC_train;
% end

% function [PC_fold, Dev_fold] = il_tbtpred_inc(cfg_ACI,results,data_passation)
% % Adapted version by Alejandro

N_lambdas = length(res.lambdas);
N_folds = cfg_ACI.N_folds;

N = length(res.FitInfo.PC_test_t);
PC_all = nan([N_lambdas, N_folds, N]);
Dev_all = nan(size(PC_all));

idxs_inc  = find(data_passation.is_correct(cfg_ACI.idx_analysis)==0);
idxs_corr = find(data_passation.is_correct(cfg_ACI.idx_analysis)==1); % correct indexes
                
for i_lambda = 1:N_lambdas
    for i_fold = 1:N_folds
        PC_test_t  = squeeze(res.FitInfo.PC_test_t(i_lambda,i_fold,:));
        Dev_test_t = squeeze(res.FitInfo.Dev_test_t(i_lambda,i_fold,:));
        CVtest = find(res.FitInfo.CV.test(i_fold));
        
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
        
        % setting back to NaN those trials within the fold that were correct:
        PC_all(i_lambda,i_fold,idx_corr_sort)  = nan; 
        Dev_all(i_lambda,i_fold,idx_corr_sort) = nan;
    end
end

PA_fold  = 100*nanmean(PC_all,3); % across third dimension (the trial dimension)
Dev_fold = nanmean(Dev_all,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num1 = Dev_fold; % averaged per trial already
num2 = num1(end,:); % referential null ACI 
yvar = num1(idxlambda,:) - num2;
 
Dev      = yvar;
Dev_mean = mean(yvar);
Dev_SEM  = 1.64*sem(yvar);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num1 = PA_fold; % PC_raw
num2 = num1(end,:); % referential null ACI 
factor_guess = repmat(100./(100-PA_fold(end,:)),N_lambdas,1); % correction for guessing, added on 9/10/2022
yvar = factor_guess.*( num1(idxlambda,:)-num2 );
PA_re_chance = mean(yvar);
PA_mean_re_chance = mean(PA_re_chance);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
outs.PA_fold  = PA_fold;
outs.Dev_fold = Dev_fold;

% outs.PA_mean           = PA_mean;
outs.PA_re_chance      = PA_re_chance;
outs.PA_mean_re_chance = PA_mean_re_chance;
outs.Dev               = Dev;
outs.Dev_mean          = Dev_mean;
outs.Dev_SEM           = Dev_SEM;
