function definput = arg_fastACI_experiment(definput)
% function definput = arg_fastACI_experiment(definput)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% definput.flags.plot = {'plot','no_plot'};

definput.keyvals.Ni  = []; % it will use i_current according to the stored experimental results
definput.keyvals.Nf  = []; % it will use the N from the experiment
definput.keyvals.file_model_decision_config  = [];