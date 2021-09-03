function def_sim = osses2021_cfg

def_sim.modelname = 'osses2021';
def_sim.decision_script = 'aci_detect';
def_sim.template_script = 'model_template';
def_sim.template_every_trial = 0;
def_sim.templ_num = 10;
def_sim.det_lev = -6;
def_sim.type_decision = 'optimal_detector';
switch def_sim.type_decision
  case 'optimal_detector'
    def_sim.thres_for_bias = 0;
end
