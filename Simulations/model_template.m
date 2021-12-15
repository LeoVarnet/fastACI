function [templ_tar,templ_ref,cfg_sim] = model_template(cfg_game,data_passation,cfg_sim)
%MODEL_TEMPLATE
% 
% Based on casp_template, add-on of the AFC toolbox.
%
% 1. Description:
%       Generate a template for the optimal detector
%
%  This routine is called automatically by the AFC framework when any of
%  the CASP models is used as a test person. It generate the template
%  needed for the optimal detector.
%
% Author: Alejandro Osses, ENS, 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i_current    = data_passation.i_current;
store_expvar = data_passation.expvar(i_current);
n_stim       = data_passation.n_stim(i_current);

%%% Temporally setting the current variable to the suprathreshold level:
data_passation.expvar(i_current) = cfg_sim.det_lev; 

if isnan(cfg_sim.det_lev)
    data_passation.expvar(i_current) = 0; 
    disp('SNR of 0 dB will be used to assess the clean-speech template...')
    
    bClean_target = 1;
    bNoisy_target = ~bClean_target;
else
    bNoisy_target = 1;
    bClean_target = ~bNoisy_target; % i.e., noisy target, the default
    
end

fs = cfg_game.fs;
modelname = cfg_sim.modelname;
[modelpars,subfs] = model_params(modelname,fs);
cfg_sim.subfs = subfs;
%%%
str_stim = [];
cfg_local = cfg_game;

if isfield(cfg_local,'seeds_order')
    % During the template derivation different (random) seeds are used
    cfg_local = Remove_field(cfg_local,'seeds_order'); 
end
if isfield(cfg_local,'Rove_level')
    % No rove level allowed
    cfg_local = Remove_field(cfg_local,'Rove_level'); 
end

cfg_local.n_targets_sorted(n_stim) = 2; % target
eval(['str_stim = ' cfg_game.experiment '_user(cfg_local,data_passation);']); % calls user-function of the experiment
if bNoisy_target
    % default 
    signal1 = str_stim.tuser;
end
if bClean_target
    if isfield(str_stim,'stim_tone_alone')
        signal1 = str_stim.stim_tone_alone;
    else
        error('This experiment file does not seem to have the possibility to derive templates from clean signals...')
    end
end

cfg_local.n_targets_sorted(n_stim) = 1; % reference
eval(['str_stim = ' cfg_game.experiment '_user(cfg_local,data_passation);']); % calls user-function of the experiment
if bNoisy_target
    % default 
    signal2 = str_stim.tuser;
end
if bClean_target
    signal2 = str_stim.stim_tone_alone;
end
%%%

[ir_signal,params] = model_representation(signal1,modelname,modelpars);
ir_reference       = model_representation(signal2,modelname,modelpars);

[templ_ref,sizeIR] = Ensure_intrep_is_numeric(ir_reference);
templ_tar          = Ensure_intrep_is_numeric(ir_signal);

for i = 1:(cfg_sim.templ_num - 1)
  
    %%%
    str_stim = [];
    cfg_local = cfg_game;
    cfg_local.n_targets_sorted(n_stim) = 2; % target
    eval(['str_stim = ' cfg_game.experiment '_user(cfg_local,data_passation);']); % calls user-function of the experiment
    signal1 = str_stim.tuser;

    cfg_local.n_targets_sorted(n_stim) = 1; % reference
    eval(['str_stim = ' cfg_game.experiment '_user(cfg_local,data_passation);']); % calls user-function of the experiment
    signal2 = str_stim.tuser;
    %%%

    ir_signal    = model_representation(signal1,modelname,modelpars);
    ir_reference = model_representation(signal2,modelname,modelpars);
    
    ir_reference = Ensure_intrep_is_numeric(ir_reference);
    ir_signal    = Ensure_intrep_is_numeric(ir_signal);
    
    templ_tar  = templ_tar  + ir_signal;
    templ_ref  = templ_ref + ir_reference;
end

templ_tar = templ_tar / cfg_sim.templ_num;  % rescaling
templ_ref = templ_ref / cfg_sim.templ_num; % reference is built

% templ_target  = templ_target - reference; % template is built

% The following is a method to simulateously scale both templates to unit energy
[~,c] = Normalise_signal(1/sqrt(2)*[templ_tar; templ_ref],subfs);
templ_tar = c*templ_tar;
templ_ref = c*templ_ref;

disp('')

% restore work.expvaract to experiment_cfg startvar:
data_passation.expvar(i_current) = store_expvar;

disp('template calculation finished');

