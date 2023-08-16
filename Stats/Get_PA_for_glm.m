function [FitInfo, ACI_sub] = Get_PA_for_glm(Data_matrix,y)
% function [FitInfo, ACI_sub] = Get_PA_for_glm(Data_matrix,y)
%
% This script calculates the mean square error (MSE) or the percentage of
%  Each row of Data_matrix is a different trial
%  y - expected output
%  ACI_sub - ACI fitted using only test data (but not reshaped yet).
%
% Author: Alejandro Osses
% Date: 9/08/2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X = Data_matrix;

idx90 = round(.9*size(X,1));
idx_random = randperm(size(X,1));
idx_test = idx_random(idx90+1:end);
idx_random = idx_random(1:idx90);

X_train = X(idx_random,:);
y_train = y(idx_random);
[ACI_sub,Dev_sub,Stat_sub] = glmfit(X_train,y_train,'binomial','link','logit');
% ACI_sub = reshape(ACI_sub(2:end),ACI_size); % ACI_size is the size of the freq and time dimensions

X_test = X(idx_test,:);

y_hat = glmval(ACI_sub,X_test,'logit'); 

Decision = round(y_hat); % i.e., above 0.5 approximated to 1, below to 0
y_ref = y(idx_test);
PC = mean(Decision(:)==y_ref(:));

FitInfo.PC = PC;
FitInfo.idx_test = idx_test;
FitInfo.idx_train = idx_random;
