function [crosspred,outs] = Read_crosspred(file)
% function [crosspred,outs] = Read_crosspred(file)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

crosspred = [];
info_toolbox = [];
load(file);

outs.info_toolbox = info_toolbox;

N_subjects = length(crosspred);
for i = 1:N_subjects
    if crosspred(i).bProcessed
        dev_here = mean(crosspred(i).Dev_test,2);
        
        % Getting the lambda that minimises the deviance:
        [~,idxlambda(i)] = min(dev_here); 

        xvar = crosspred(i).lambdas;
        PC_test = crosspred(i).PC_test;
        Dev_test = crosspred(i).Dev_test;
        
        yvar = mean(PC_test-PC_test(end,:),2); % avoiding to write down 'results{1}.crosspred.PC_test' twice
        evar = std(PC_test-PC_test(end,:),[],2);%/size(PC_test,2);%SEM

        PA_mean_re_chance(i) = 100*yvar(idxlambda(i));
        
        yvar = mean(PC_test,2); % avoiding to write down 'results{1}.crosspred.PC_test' twice
        PA_mean(i) = 100*yvar(idxlambda(i));
        disp('')

    else
        idxlambda(i) = nan;
        PA_mean(i) = nan;
        PA_mean_re_chance(i) = nan;
    end
end

outs.PA_mean           = PA_mean;
outs.PA_mean_re_chance = PA_mean_re_chance;