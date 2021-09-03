function def_sim = fastACI_set_simulation_config(modelname,def_sim)
% function def_sim = fastACI_set_simulation_config(modelname,def_sim)
%
% This script is related to model_cfg_replace.txt.
%
% def_sim.modelname = '$$modelname$$';
% def_sim.decision_script = '$$decision_script$$';
% def_sim.template_script = '$$template_script$$';
% def_sim.template_every_trial = $$template_every_trial$$;
% def_sim.templ_num = $$templ_num$$;
% def_sim.det_lev = $$det_lev$$;
% def_sim.type_decision = '$$type_decision$$';
% switch def_sim.type_decision
%   case 'optimal_detector'
%     def_sim.thres_for_bias = $$thres_for_bias$$;
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    def_sim = [];
end

def_sim.modelname = modelname;
bInput = input('Enter 0 to load defaults for the current auditory model, enter 1 to input parameter by parameter: ');

if bInput == 0
    switch modelname
        case 'king2019'
            def_sim.decision_script = 'king2019_detect'; 
            def_sim.template_script = 'king2019_template'; % this can be later automated
            def_sim.type_decision   = '';
            
        case 'osses2021'
            def_sim.decision_script = 'aci_detect';
            def_sim.template_script = 'model_template';
            def_sim.type_decision   = 'optimal_detector';
            
        otherwise
            error('Add model %s to the list',modelname);
    end
else
    % decision_script:
    script_options = {'aci_detect','king2019_detect'};
    Show_cell(script_options);
    idx = input('Choose the decision_script from the list above: ');
    def_sim.decision_script = script_options{idx};
    
    % template_script:
    script_options = {'model_template','king2019_template'};
    Show_cell(script_options);
    idx = input('Choose the template_script from the list above: ');
    def_sim.template_script = script_options{idx};
  
    % type_decision:
    switch def_sim.decision_script
        case 'aci_detect'
            script_options = {'optimal_detector','relanoiborra2019_decision'};
            Show_cell(script_options);
            idx = input('Choose the type_decision from the list above: ');
            def_sim.type_decision = script_options{idx};
            
        otherwise
            def_sim.type_decision = NaN;
    end
    
    switch type_decision
        case 'optimal_detector'
            def_sim.thres_for_bias = input('Enter the thres_for_bias (default=0): ');
        otherwise
            def_sim.thres_for_bias = 0;
    end
    
    def_sim.template_every_trial = input('template_every_trial (0=no - default, 1=yes): ');
    def_sim.templ_num = input('templ_num, number of averages for the template (1 or typically 10): ');
    def_sim.det_lev = input('det_lev, suprathreshold level to derive the template (task-dependenet): ');
end