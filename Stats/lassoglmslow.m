function [B,FitInfo] = lassoglmslow(X,y,N_folds,lambdas)
% function [B,FitInfo] = lassoglmslow(X,y,N_folds,lambdas)
%
% Author: Leo Varnet
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 4
    N_lambda = 30;
    FitInfo.Lambda = logspace(-4, -1, N_lambda); % lambdas;
else
    N_lambda = length(lambdas);
    FitInfo.Lambda = lambdas;
end

CV = cvpartition(size(X,1),'KFold',N_folds);
fprintf('\n')

% Memory allocation:
FitInfo.yhat_train = nan(N_lambda,N_folds,max(CV.TrainSize));
FitInfo.yhat_test  = nan(N_lambda,N_folds,max(CV.TestSize));
FitInfo.PCtest_t   = nan(N_lambda,N_folds,max(CV.TestSize));
FitInfo.MSEtest_t  = nan(N_lambda,N_folds,max(CV.TestSize));

for i_lambda = 1:N_lambda
    fprintf(['Slow lassoglm, computing lambda ' num2str(i_lambda) ' of ' num2str(N_lambda) '\n'])
    for i_fold = 1:N_folds
        
        idx_train = CV.training(i_fold); % idxs of the training set in this fold
        idx_test  = CV.test(i_fold); % idxs for the test (validation) in this fold
        
        [B_temp,FitInfo_temp] = lassoglm(X(idx_train,:),y(idx_train),'binomial','Lambda',FitInfo.Lambda(i_lambda));
        FitInfo.Intercept(i_lambda,i_fold) = FitInfo_temp.Intercept;
        FitInfo.DF(i_lambda,i_fold) = FitInfo_temp.DF;
        FitInfo.Devtrain(i_lambda,i_fold) = FitInfo_temp.Deviance;
        B(:,i_lambda,i_fold ) = B_temp;
        
        coef = [FitInfo_temp.Intercept; B_temp];
        
        %%% Training:
        yhat_train = glmval(coef,X(idx_train,:),'logit');
        [PC,MSE,Dev,MSE_rounded] = Get_prediction_metrics(yhat_train,y,idx_train);
        % Devtrain = Dev;
        
        FitInfo.MSEtrain(i_lambda,i_fold) = MSE;
        FitInfo.PCtrain(i_lambda,i_fold)  = PC;
        FitInfo.yhat_train(i_lambda,i_fold,1:length(yhat_train)) = yhat_train;
        
        %%% Test (or validation):
        yhat_test = glmval(coef,X(idx_test,:),'logit'); % X(CV.test(i_fold),:)*B_temp + FitInfo_temp.Intercept;
        [PC,MSE,Dev,MSE_rounded, yhat_test_rounded, PC_t, MSE_t, Dev_t] = Get_prediction_metrics(yhat_test,y,idx_test);
        
        FitInfo.Devtest(i_lambda,i_fold) = Dev;
        FitInfo.MSEtest(i_lambda,i_fold) = MSE;
        FitInfo.PCtest(i_lambda,i_fold)  = PC;
        FitInfo.yhat_test(i_lambda,i_fold,1:length(yhat_test)) = yhat_test;
        FitInfo.PCtest_t(i_lambda,i_fold,1:length(yhat_test))  = PC_t;
        FitInfo.Devtest_t(i_lambda,i_fold,1:length(yhat_test)) = Dev_t;
    end
end

% Temp solution
FitInfo.y = y;
FitInfo.B = B;
%B = mean(mean(B,3),2);
B = mean(B,3);

FitInfo.CV = CV;
end