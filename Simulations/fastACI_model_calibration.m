function keyvals = fastACI_model_calibration(experiment, model, Condition, keyvals)
% function fastACI_model_calibration(experiment,model)
%
% Calibration of a model: This script assesses the criterion value K.
%
% Original name: g20211207_calibrating_the_model.m (snapshot on 14/04/2022)
% Author: Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch experiment
    case 'speechACI_Logatome-abda-S43M'
        suff_exp = '';
        
    case 'speechACI_varnet2013'
        
        validated_conds = {'white','SSN'};
        Show_cell(validated_conds);
        
        bInput = input('Enter the number of the conditions you want to test (one or more numeric values are accepted, or enter 0 for all): ');
        if bInput == 0
            Conds = validated_conds;
        else
            Conds = validated_conds(bInput); 
        end
        
        suff_exp = ['-' experiment];
        
        if ~strcmp(Conds{1},'white')
            suff_exp = [suff_exp '-' Conds{1}];
        end
end

expvar = -40; % This is a very low experimental variable
switch model
    case 'king2019'
        bias_step = (5e-5)/100; % based on Osses2022_preprint Fig. 10B
    case 'maxwell2020'
        % bias_step = 4; % for normal configuraton
        bias_step = 0.1;
    otherwise
        bias_step = 0.1; % 1 MU over 10
end

bRe_init = 1;

% end
bContinue = 1;
if bRe_init
    thres_for_bias = -1.4813; % -1.2173; % 0;
    in_std = 0;
    % in_std = 3.14;
    % save(outfile,'thres_for_bias','in_std');
end
%%% 
Ni = 1;
N = 400;

Iterations = 1;
while bContinue == 1 && Iterations <= 20
    Nf = Ni + N;
    flags_here = {'Ni',Ni,'Nf',Nf,'thres_for_bias',thres_for_bias}; % ,'file_model_decision_config',file_model_decision_config};

    try
        [cfg_game, data_passation] = fastACI_experiment_constant(experiment,model,Condition,expvar,flags_here{:});
    catch
        %%% Run the following first, to create a create file
        % [cfg_game, data_passation] = fastACI_experiment(experiment,model,noise_cond);

        % It means that the experiment has to be initialised first
        fastACI_experiment_init(experiment,model, Condition);
        [cfg_game, data_passation] = fastACI_experiment_constant(experiment,model,Condition,expvar,flags_here{:});
    end

    if Iterations == 1
        % First adjustment, based on one session:
        mues = data_passation.decision_var_mue2choose;
        diffe = mues(:,2)-mues(:,1);
        thres_for_bias = median(diffe);
    end

    idxs = 1:N;
    bias_r1(Iterations) = 100*sum(data_passation.n_responses(idxs)==1)/N;
    bias_r2(Iterations) = 100*sum(data_passation.n_responses(idxs)==2)/N;
    score(Iterations)   = 100*sum(data_passation.is_correct(idxs)==1)/N;

    thres_tested(Iterations) = thres_for_bias;
    if bias_r1(Iterations) > 54
        if Iterations ~= 1
            thres_for_bias = thres_for_bias-bias_step;
        end
    elseif bias_r1(Iterations) < 46
        if Iterations ~= 1
            thres_for_bias = thres_for_bias+bias_step;
        end
    else
        bContinue = 0;
    end

    Iterations = Iterations + 1;

    Ni = Ni + N;
    if Ni > cfg_game.N
        Ni = 1;
    end
end

var2show = [thres_tested; bias_r1; bias_r2; score];
% var2latex(var2show);

res = var2show;

% thres_bias(j) = thres_for_bias;
date = Get_date;
date = date.date2print;

%%%
fname = ['optimal_detector-' model '-' Condition '-' date];

pars = [];
pars.thres_for_bias = thres_for_bias;
pars.in_std = in_std;
pars.description = ['Calibration using ' mfilename];

Save_model_calibration(fname,pars);

%%%

keyvals.thres_for_bias = thres_for_bias;
keyvals.in_std         = in_std;
