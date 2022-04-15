function basepath = fastACI_file_template(experiment_full, modelname, type_decision, keyvals)
% function basepath = fastACI_file_template(experiment_full, modelname, type_decision, keyvals)
%
% See also fastACI_experiment.m (~L296), fastACI_experiment_constant.m, 
%     aci_detect.m (~L91)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

template_str = keyvals.fname_template_suffix;
if ~isempty(template_str)
    if ~strcmp(template_str(1),'-') % compares if the first string is '-' and adds it 
                              % if it is not. This improves the readability of the 
                              % template file name
        template_str = ['-' template_str];
    end
end

basepath = [fastACI_paths('dir_data') experiment_full filesep modelname filesep ...
    'template-' modelname '-' experiment_full '-trial-1' template_str '.mat'];
