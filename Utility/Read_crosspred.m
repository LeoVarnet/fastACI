function [crosspred,outs] = Read_crosspred(file, idxlambda_opt)
% function [crosspred,outs] = Read_crosspred(file, idxlambda)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    bLook_for_lambda = 1;
    fprintf('%s: Looking for the best idxlambda value for each cross prediction\n',mfilename);
else
    bLook_for_lambda = 0;
    fprintf('%s: Using a fixed idxlambda value\n',mfilename);
end
crosspred = [];
info_toolbox = [];
load(file);

outs.info_toolbox = info_toolbox;

N_subjects = length(crosspred);
for i = 1:N_subjects
    if crosspred(i).bProcessed
        dev_here = mean(crosspred(i).Dev_test,2);
        
        if bLook_for_lambda
            % Getting the lambda that minimises the deviance:
            [~,idxlambda(i)] = min(dev_here); 
        else
            idxlambda(i) = idxlambda_opt;
        end
        xvar = crosspred(i).lambdas;
        PC_test  = crosspred(i).PC_test;
        Dev_test = crosspred(i).Dev_test;
        Dev_test_t = crosspred(i).Dev_test_t;
        
        yvar = mean(PC_test-PC_test(end,:),2); % avoiding to write down 'results{1}.crosspred.PC_test' twice
        evar = std(PC_test-PC_test(end,:),[],2);%/size(PC_test,2);%SEM

        PA_mean_re_chance(i) = 100*yvar(idxlambda(i));
        
        yvar = mean(PC_test,2); % avoiding to write down 'results{1}.crosspred.PC_test' twice
        PA_mean(i) = 100*yvar(idxlambda(i));
        
        %%%
        yvar = Dev_test(idxlambda(i),:)-Dev_test(end,:);
        
        %%%
        Dev_trial_fold = squeeze(Dev_test_t(idxlambda(i),:,:));
        N_folds = size(Dev_trial_fold,1);
        Nt = nan([1 N_folds]);
        for kk = 1:N_folds
            Nt(kk) = sum(~isnan(Dev_trial_fold(kk,:)));
        end
        yvar = yvar./Nt; % average per trial
        Dev_mean(i) = mean(yvar); % mean 
        Dev_SEM(i)  = 1.64*sem(yvar);
    else
        idxlambda(i) = nan;
        PA_mean(i) = nan;
        PA_mean_re_chance(i) = nan;
        
        Dev_mean(i) = nan;
        Dev_SEM(i) = nan;
    end
end

outs.PA_mean           = PA_mean;
outs.PA_mean_re_chance = PA_mean_re_chance;
outs.Dev_mean          = Dev_mean;
outs.Dev_SEM           = Dev_SEM;