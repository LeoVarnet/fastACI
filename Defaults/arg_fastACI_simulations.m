function definput = arg_fastACI_simulations(definput)
% function definput = arg_fastACI_simulations(definput)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% definput.keyvals.file_model_decision_config  = [];
definput.keyvals.in_std = 0; % no standard deviation for the optimal detector
definput.keyvals.thres_for_bias = []; % no K value
definput.keyvals.fname_template_suffix = [];
definput.keyvals.model_cfg_script = []; % empty by default