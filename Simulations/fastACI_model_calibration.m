function keyvals = fastACI_model_calibration(experiment, model, Condition, flags,keyvals)
% function fastACI_model_calibration(experiment,model)
%
% Calibration of a model: This script assesses the criterion value K.
%
% Original name: g20211207_calibrating_the_model.m (snapshot on 14/04/2022)
% Author: Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% definput.import={'fastACI_experiment','fastACI_simulations','fastACI_simulation_detect'};
% [flags,keyvals]  = ltfatarghelper({},definput,varargin);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isfield(keyvals,'expvar_cal')
    expvar_cal = keyvals.expvar_cal;
else
    expvar_cal = -40; % This is a very low experimental variable
    fprintf('%s: Using default expvar value of %.1f dB\n',upper(mfilename),expvar_cal);
end

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
    thres_for_bias = 0; % -1.2173; % 0;
    in_std = 0;
    % in_std = 3.14;
    % save(outfile,'thres_for_bias','in_std');
end
%%% 
if isempty(keyvals.Ni)
    Ni = 1; 
else
    Ni = keyvals.Ni;
end

if isempty(keyvals.Nf)
    N = 400;
else
    N = keyvals.Nf;
end

Iterations = 1;
while bContinue == 1 && Iterations <= 20 
    Nf = Ni + N;
    flags_here = Get_idle_flag;
        
    kv_here = keyvals;
    kv_here.Ni = Ni;
    kv_here.Nf = Nf;
    kv_here.thres_for_bias = thres_for_bias;
    
    try
        [cfg_game, data_passation] = fastACI_experiment_constant(experiment,model,Condition,expvar_cal,'argimport',flags_here,kv_here);
    catch
        %%% Run the following first, to create a create file
        % [cfg_game, data_passation] = fastACI_experiment(experiment,model,noise_cond);

        % It means that the experiment has to be initialised first
        fastACI_experiment_init(experiment,model, Condition);
        [cfg_game, data_passation] = fastACI_experiment_constant(experiment,model,Condition,expvar_cal,flags_here{:});
    end

    if Iterations == 1
        % First adjustment, based on one session:
        mues = data_passation.decision_var_mue2choose;
        diffe = mues(:,2)-mues(:,1);
        thres_for_bias = median(diffe);
    end

    Play_ready;
    
    idxs = Ni:Nf-1;
    bias_r1(Iterations) = 100*sum(data_passation.n_responses(idxs)==1)/N;
    bias_r2(Iterations) = 100*sum(data_passation.n_responses(idxs)==2)/N;
    score(Iterations)   = 100*sum(data_passation.is_correct(idxs)==1)/N;

    thres_tested(Iterations) = thres_for_bias;
    if mod(Iterations,2) == 1
        factor = 0.8; % arbitrary to avoid to fall in a loop
    else
        factor = 1;
    end
    if bias_r1(Iterations) > 54
        if Iterations ~= 1
            thres_for_bias = thres_for_bias-factor*bias_step;
        end
    elseif bias_r1(Iterations) < 46
        if Iterations ~= 1
            thres_for_bias = thres_for_bias+factor*bias_step;
        end
    else
        bContinue = 0;
    end

    Iterations = Iterations + 1;

    Ni = Ni + N;
    if Ni >= cfg_game.N
        % Ni = 1;
        break;
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
