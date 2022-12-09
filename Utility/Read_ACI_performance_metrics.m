function outs = Read_ACI_performance_metrics(res, type)
% function outs = Read_ACI_performance_metrics(res, type)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 2
    type = 'test';
end
idxlambda = res.idxlambda;

switch type
    case 'test'
        Dev_raw = res.FitInfo.Dev_test;
        CV_size = res.FitInfo.CV.TestSize;
        PC_raw = 100*res.FitInfo.PC_test;
    case 'train'
        Dev_raw = res.FitInfo.Dev_train;
        CV_size = res.FitInfo.CV.TrainSize;
        PC_raw = 100*res.FitInfo.PC_train;
end
num1 = Dev_raw./CV_size; % averaged per trial
num2 = num1(end,:); % referential null ACI 
yvar = num1(idxlambda,:) - num2;

Dev      = yvar;
Dev_mean = mean(yvar);
Dev_SEM  = 1.64*sem(yvar);
%%%

N_lambdas = size(PC_raw,1);
num1 = PC_raw;
num2 = num1(end,:); % referential null ACI 
factor_guess = repmat(100./(100-PC_raw(end,:)),N_lambdas,1); % correction for guessing, added on 9/10/2022
yvar = factor_guess.*( num1(idxlambda,:)-num2 );
PA_re_chance = mean(yvar);
PA_mean_re_chance = mean(PA_re_chance);

% outs.PA_mean           = PA_mean;
outs.PA_re_chance      = PA_re_chance;
outs.PA_mean_re_chance = PA_mean_re_chance;
outs.Dev               = Dev;
outs.Dev_mean          = Dev_mean;
outs.Dev_SEM           = Dev_SEM;