function [PC, MSE, Dev, MSE_rounded, yhat_rounded, PC_t, MSE_t, Dev_t] = Get_prediction_metrics(yhat,y,idx_train_or_test)
% function [PC, MSE, Dev, MSE_rounded, yhat_rounded, MSE_t, Dev_t] = Get_prediction_metrics(yhat,y,idx_train_or_test)
%
% This script calculates the mean square error (MSE) or the percentage of
%     coincidence between the variable 'yhat' and the variable 'y'. These
%     metrics are actually calculated between yhat and y(idx_train_or_test).
% This script is used in lassoslow.m and lassoglmslow.m.
%
% PC
% MSE
% Dev  : Deviance
% Dev_t: Deviance as a function of time sample
%
% See also: lassoslow.m, lassoglmslow.m
%
% Authors: Leo Varnet and Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

yhat_rounded = (yhat>=0.5); % equivalent to round(yhat). It assumes a hard
                            % threshold at 0.5
MSE_t = (y(idx_train_or_test) - yhat).^2;
MSE = mean(MSE_t);
MSE_rounded = mean((y(idx_train_or_test) - yhat_rounded).^2);

PC_t = yhat_rounded==y(idx_train_or_test);
PC = mean(PC_t);

%%%
% yref = y(idx_train_or_test); % reference values
% n_of_trials = 1; % 
% pdf_train_or_test = binopdf(yref,n_of_trials,yhat)
% pdf_saturated     = binopdf(yref,n_of_trials,yref) % saturated model?

% lassoglmslow
Dev_t = -2*(log(binopdf(y(idx_train_or_test),1,yhat))) - sum(log(binopdf(y(idx_train_or_test),1,y(idx_train_or_test))));
Dev   = sum(Dev_t);
