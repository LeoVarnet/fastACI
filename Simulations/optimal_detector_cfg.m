function pars = optimal_detector_cfg(modelname,keyvals)

% modelname_short = strsplit(modelname,'_');
% modelname_short = modelname_short{1};

% dir_where = [fastACI_basepath 'Interim_results' filesep];
% if isempty(keyvals.file_model_decision_config)
%     outfile = [dir_where modelname_short '-optimal_detector.mat'];
% else
%     outfile = keyvals.file_model_decision_config;
% end

if ~isempty(keyvals.in_std)
    in_std = keyvals.in_std;
else
    in_std = input('Enter the value for in_std (MU) - press 0 if you don''t know:');
end
pars.in_std = in_std;

if ~isempty(keyvals.thres_for_bias)
    thres_for_bias = keyvals.thres_for_bias;
else
    thres_for_bias = input('Enter the value for thres_for_bias (MU) - press 0 if you don''t know:');
end
if ~isempty(keyvals.thres_for_bias_each_session)
    pars.thres_for_bias_each_session = keyvals.thres_for_bias_each_session;
end
pars.thres_for_bias = thres_for_bias;

in_std = pars.in_std;
in_var = in_std*in_std;
pars.in_var = in_var;
