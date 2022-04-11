function [B,FitInfo] = lassoslow(X,y,N_folds,lambdas)
% function [B,FitInfo] = lassoslow(X,y,N_folds,lambdas)
%
% This function requires the lasso.m implementation.
%
% Author: Leo Varnet
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 4
    N_lambda = 30;
    FitInfo.Lambda = logspace(-4, -1, N_lambda);%lambdas;
else
    N_lambda = length(lambdas);
    FitInfo.Lambda = lambdas;
end

CV = cvpartition(size(X,1),'KFold',N_folds);%cvpartition(y,'KFold',N_folds);%
fprintf('\n')

FitInfo.yhat_train = nan(N_lambda,N_folds,max(CV.TrainSize)); % memory allocation
FitInfo.yhat_test  = nan(N_lambda,N_folds,max(CV.TestSize)); % memory allocation
FitInfo.PCtest_t   = nan(N_lambda,N_folds,max(CV.TestSize)); % memory allocation
FitInfo.MSEtest_t  = nan(N_lambda,N_folds,max(CV.TestSize)); % memory allocation

for i_lambda = 1:N_lambda
    fprintf(['Slow lasso, computing lambda ' num2str(i_lambda) ' of ' num2str(N_lambda) '\n'])
    for i_fold = 1:N_folds
        
        idx_train = CV.training(i_fold); % idxs of the training set in this fold
        idx_test  = CV.test(i_fold); % idxs for the test (validation) in this fold
        [B_temp,FitInfo_temp] = lasso(X(idx_train,:),y(idx_train),'Lambda',FitInfo.Lambda(i_lambda));
        FitInfo.Intercept(i_lambda,i_fold) = FitInfo_temp.Intercept;
        FitInfo.DF(i_lambda,i_fold)        = FitInfo_temp.DF;
        FitInfo.MSE_train(i_lambda,i_fold)  = FitInfo_temp.MSE;
        B(:,i_lambda,i_fold ) = B_temp;
        
        %%% Training:
        yhat_train = X(idx_train,:)*B_temp + FitInfo_temp.Intercept;
        [PC,MSE,Dev,MSE_rounded] = Get_prediction_metrics(yhat_train,y,idx_train);
        
        FitInfo.MSE_train(i_lambda,i_fold) = MSE;
        FitInfo.PC_train(i_lambda,i_fold)  = PC;
        FitInfo.yhat_train(i_lambda,i_fold,1:length(yhat_train)) = yhat_train;
        
        %%% Test (or validation):
        yhat_test = X(idx_test,:)*B_temp + FitInfo_temp.Intercept;
        [PC,MSE,Dev,MSE_rounded,yhat_test_rounded,PC_t,MSE_t,Dev_t] = Get_prediction_metrics(yhat_test,y,idx_test);
        
        FitInfo.MSE_test(i_lambda,i_fold) = MSE;
        FitInfo.PC_test(i_lambda,i_fold)  = PC;
        FitInfo.yhat_test(i_lambda,i_fold,1:length(yhat_test)) = yhat_test;
        
        FitInfo.PC_test_t(i_lambda,i_fold,1:length(yhat_test))  = PC_t; % yhat_test_rounded==y(idx_test);
        FitInfo.MSE_test_t(i_lambda,i_fold,1:length(yhat_test)) = MSE_t; % (y(idx_test) - yhat_test).^2;

    end
end

% Temp solution
FitInfo.y = y;
FitInfo.B = B;
%B = mean(mean(B,3),2);
B = mean(B,3);

FitInfo.CV = CV;

end