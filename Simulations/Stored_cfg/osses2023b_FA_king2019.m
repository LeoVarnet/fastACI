function def_sim = king2019_cfg(keyvals)

% Parameters for the model (1 of 2)
def_sim.modelname = 'king2019';

% Parameters for the detector:
def_sim.decision_script = 'aci_detect';
def_sim.template_script = 'model_template';
def_sim.template_every_trial = 0;
def_sim.templ_num = 100;
def_sim.det_lev = 0; % 0 dB SNR for this tone-in-noise task
def_sim.type_decision = 'optimal_detector';
switch def_sim.type_decision
  case 'optimal_detector'
    optdet_params = optimal_detector_cfg(def_sim.modelname,keyvals);
    def_sim.thres_for_bias = optdet_params.thres_for_bias;
    def_sim.in_var = optdet_params.in_var;
end
def_sim.bStore_template = 1;

% Parameters for the model (2 of 2)
dBFS = 100; % full scale convention

modelpars = {};
basef = 500; % Hz, to match the basef from toneinnoise_ahumada1971 and toneinnoise_ahumada1975
flow = 80;
fhigh = 8000;
% modbank_Nmod = 5;
mflow  =   2; % Hz, modbank_fmin
mfhigh = 150;

modelpars{end+1} = 'no_debug'; % flag
modelpars{end+1} = 'phase_insens_hilbert'; % 'no_phase_insens'; % flag
% default values will be loaded for king2019:
modelpars(end+1:end+2) = {'compression_n',0.3};
modelpars(end+1:end+2) = {'basef', basef}; % keyval
modelpars(end+1:end+2) = {'flow' ,  flow}; % keyval
modelpars(end+1:end+2) = {'fhigh', fhigh}; % keyval
modelpars(end+1:end+2) = {'mflow' , mflow}; % keyval
modelpars(end+1:end+2) = {'mfhigh', mfhigh}; % keyval        
modelpars(end+1:end+2) = {'dboffset',100}; 
def_sim.modelpars = modelpars;

if exist('model_cfg.m','file')
    % Config: AFC toolbox
    global def
    
    model_cfg;
end
