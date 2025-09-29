function Save_model_calibration(fname,pars,dir_results)
% function Save_model_calibration(fname,pars,dir_results)
%
% See also:
%   g20211207_calibrating_the_model.m: first model+decision calibration 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fields_str = fieldnames(pars);

vars = [];
for i = 1:length(fields_str)
    vars = [vars '''' fields_str{i} ''','];
    
    exp2eval = sprintf('%s=pars.%s;',fields_str{i},fields_str{i});
    eval(exp2eval); % loading each variable into the local workspace
end
vars = vars(1:end-1); % deleting the last comma

if nargin < 3
    dir_results = [fastACI_basepath 'Interim_results' filesep];
end
if exist(dir_results,'dir')
    fname_out = [dir_results fname];
    if ~isoctave
        % If MATLAB is used:
        exp2eval = sprintf('save(''%s'',%s);',fname_out,vars);
    else
        % If GNU Octave is used:
        exp2eval = sprintf('save(''%s'',%s,''-mat7-binary'');',fname_out,vars);
    end
    eval(exp2eval);
else
    error('No %s directory was found, please create it and re-run this script',dir_results);
end

% Play_ready;

fprintf('%s: file successfully stored... \n',upper(mfilename),fname_out);
