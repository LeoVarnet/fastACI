function definput = arg_fastACI_experiment(definput)
% function definput = arg_fastACI_experiment(definput)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

definput.flags.idle = {'idle','no_idle'}; % No influence in anything at all
definput.keyvals.Ni  = []; % it will use i_current according to the stored experimental results
definput.keyvals.Nf  = []; % it will use the N from the experiment
% definput.keyvals.file_model_decision_config  = []; % commented on 13/04/2022
                                        % and moved to arg_fastACI_simulations.m
definput.keyvals.feedback = 1; % It can be overwritten in *_cfg
definput.keyvals.Language = 'EN'; % or 'EN' or 'FR'
definput.keyvals.randorder = 1; % DEAFAULT=1. if 0 the trials in data_passation will be ordinally sorted
definput.keyvals.sessionsN = 400; % default test block duration