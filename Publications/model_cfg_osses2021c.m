function p = model_cfg_osses2021c(run_str,modelname,templ_num,bStore_template)
% function p = model_cfg_osses2021c(run_str,modelname,templ_num,bStore_template)
%
% Former in-line name:
%   il_get_model_config_DAGA in publ_osses2021c_DAGA_1_sim.m
%
% Author: Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 4
    bStore_template = 0;
end
if nargin < 3
    templ_num = 10; % As used in the publication
end

p = [];
p.modelname = modelname;
if isempty(strfind('_',modelname)) % contains(modelname,'_') 
    modelname_script = modelname;
else
    % The case for relanoiborra2019:
    modelname_script = strsplit(modelname,'_');
    modelname_script = modelname_script{1};
end
p.modelname_script = modelname_script;
p.templ_num = templ_num;
p.bStore_template = bStore_template; % New option added on 17/09/2021, not 
                                     % relevant for the simulations here
p.in_std = 0; % No internal noise unless overwritten below
switch run_str
    case 'run-1'
        p.decision_script = 'aci_detect';
        p.template_script = 'model_template';
        p.template_every_trial = 0;
        p.det_lev = -6;
        p.type_decision = 'optimal_detector';
        p.thres_for_bias = 0;
        
    case 'run-3-m1p55'
        p.decision_script = 'aci_detect';
        p.template_script = 'model_template';
        p.template_every_trial = 0;
        p.det_lev = -6;
        p.type_decision = 'optimal_detector';
        p.thres_for_bias = -1.55;
        
    case 'run-3-p0p39'
        p.decision_script = 'aci_detect';
        p.template_script = 'model_template';
        p.template_every_trial = 0;
        p.det_lev = -6;
        p.type_decision = 'optimal_detector';
        p.thres_for_bias = 0.39;
        
    case 'run-3-p0p78'
        p.decision_script = 'aci_detect';
        p.template_script = 'model_template';
        p.template_every_trial = 0;
        p.det_lev = -6;
        p.type_decision = 'optimal_detector';
        p.thres_for_bias = 0.78;
        
    case 'run-4'
        p.decision_script = 'aci_detect';
        p.template_script = 'model_template';
        p.template_every_trial = 0;
        p.det_lev = -6;
        p.type_decision = 'relanoiborra2019_decision';
     
    case 'default'
        fprintf('\t%s: Default model parameters are being loaded\n',upper(mfilename));
        switch modelname
            case {'dau1997','relanoiborra2019'}
                p.decision_script = 'aci_detect';
                p.template_script = 'model_template';
                p.template_every_trial = 0;
                p.det_lev = -6;
                p.type_decision = 'optimal_detector';
                p.in_std = 0; % No calibration yet
                p.thres_for_bias = 0; 
                
            case 'osses2021'
                p.decision_script = 'aci_detect';
                p.template_script = 'model_template';
                p.template_every_trial = 0;
                p.det_lev = -6;
                p.type_decision = 'optimal_detector';
                p.in_std = 3.14; % MU, calibration on 7/12/2021
                p.thres_for_bias = 0; 
                
           case 'osses2022a'
                p.decision_script = 'aci_detect';
                p.template_script = 'model_template';
                p.template_every_trial = 0;
                p.det_lev = -6;
                p.type_decision = 'optimal_detector';
                p.in_std = 0; % No calibration yet
                p.thres_for_bias = 0; 
        end
        
    otherwise
        error('Run condition not recognised')
end

if ~isfield(p,'thres_for_bias')
    p.thres_for_bias = nan;
end

% p.bStore_template = 1;
