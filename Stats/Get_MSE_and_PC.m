function [PC, MSE, MSE_rounded, yhat_rounded] = Get_MSE_and_PC(yhat,y,idx_train_or_test)
% function [PC, MSE, MSE_rounded, yhat_rounded] = Get_MSE_and_PC(yhat,y,idx_train_or_test)
%
% This script calculates the mean square error (MSE) or the percentage of
%     coincidence between the variable 'yhat' and the variable 'y'. These
%     metrics are actually calculated between yhat and y(idx_train_or_test).
% This script is used in lassoslow.m and lassoglmslow.m.
%
% See also: lassoslow.m, lassoglmslow.m
%
% Authors: Leo Varnet and Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

yhat_rounded = (yhat>=0.5); % equivalent to round(yhat). It assumes a hard
                            % threshold at 0.5
MSE = mean((y(idx_train_or_test) - yhat).^2);
MSE_rounded = mean((y(idx_train_or_test) - yhat_rounded).^2);

PC = mean(yhat_rounded==y(idx_train_or_test));
