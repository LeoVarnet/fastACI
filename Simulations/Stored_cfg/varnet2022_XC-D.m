function def_sim = king2019_cfg

def_sim.modelname = 'king2019';
def_sim.decision_script = 'king2019_detect'; % def_sim.decision_script = 'aci_detect';
def_sim.template_script = 'king2019_template'; % 'model_template';
def_sim.template_every_trial = 0;
def_sim.templ_num = 1;
def_sim.det_lev = NaN; % clean templates will be obtained
def_sim.type_decision = 'optimal_detector';
switch def_sim.type_decision
  case 'optimal_detector'
    optdet_params = optimal_detector_cfg(def_sim.modelname);
    def_sim.thres_for_bias = optdet_params.thres_for_bias;
    def_sim.in_var = optdet_params.in_var;
end
def_sim.bStore_template = 1;

dBFS = 100;
def_sim.subfs = 16000; % Hz
basef = 1000; % Hz
def_sim.basef = basef;
flow = basef;
fhigh = basef;
fmod = 4; 
modbank_Nmod = 1;
def_sim.modelpars = {'no_debug','no_phase_insens','basef',basef,'flow',flow, ...
    'fhigh',fhigh,'mflow',fmod,'mfhigh',fmod,'modbank_Nmod',modbank_Nmod,'dboffset',dBFS};

if exist('model_cfg.m','file')
    % Config: AFC toolbox
    global def
    
    model_cfg;
end
