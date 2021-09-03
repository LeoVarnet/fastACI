function [response,sim_work] = aci_detect(cfg_game,data_passation,cfg_sim,sim_work)
% ACI_DETECT
%
% Based on casp_detect.m, an add-on function from the AFC toolbox
%
%  This routine is called automatically by the AFC framework when any of
%  the CASP models are used as a test person.
%
%  The function returns the presentation interval selected by the model
%
% Options to be adjusted by the user:
%
%       def.templateeverytrial
%       def.template_script = 'casp_template', 'casp_template_piano'
%       type_processing:
%               % 'pianoinnoise' (if pianoinnoise experiment)
%               % 'pianoinnoise2' (if casp_template_piano)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
%%% New name            Old name
%   cfg_sim             simdef

%   def.experiment      def.expname
% warning('Modelling under construction: Not validated yet...')

% global def
% global work
def = [];
if nargin < 4
    sim_work = [];
end
sim_work = Ensure_field(sim_work,'templ_tar',[]);
sim_work = Ensure_field(sim_work,'templ_ref',[]);

% global set
 
def.experiment = cfg_game.experiment;

% if def.debug == 1
%   disp([work.vpname '_detect']);
% end
% 
% % This code only handles running references using the optimal detector
% 
% def = Ensure_field(def,'templateeverytrial',0);
% def = Ensure_field(def,'template_script','casp_template');
 
template_script = cfg_sim.template_script;

cfg_sim = Ensure_field(cfg_sim,'modelname','dau1997');
cfg_sim = Ensure_field(cfg_sim,'template_every_trial',0);
cfg_sim = Ensure_field(cfg_sim,'det_lev',-6); 
cfg_sim = Ensure_field(cfg_sim,'templ_num',10);

if ~isfield(cfg_sim,'type_decision')
    warning('Assinging the default decision device for model %s',cfg_sim.modelname);
    switch cfg_sim.modelname
        case 'dau1997'
            cfg_sim.type_decision = 'optimal_detector';
        case 'osses2021'
            cfg_sim.type_decision = 'optimal_detector';
        case 'relanoiborra2019'
            cfg_sim.type_decision = 'relanoiborra2019_decision';
    end
end
type_decision = cfg_sim.type_decision; % 'relanoiborra2019_decision'; 'optimal_detector'; % type_processing;
sim_work.type_decision = type_decision;

if isfield(cfg_sim,'thres_for_bias')
    thres_for_bias = cfg_sim.thres_for_bias;
else
    switch type_decision
        case 'optimal_detector'
            disp('Non calibrated thres_for_bias')
            thres_for_bias = 0;
    end
    % switch cfg_sim.modelname
    %     case 'osses2021'
    %         thres_for_bias = 0.37; 
    %     case 'king2019'
    %         thres_for_bias = 9.2685e-07; %1.2358e-06; 
    % end
end

% switch template_script
%     case 'casp_template_v1'
%         
%         template_idx = round(10*mod(def.version_nr,1)); 
%         version_nr = def.version_nr;
%         
%         % version_nr = 1; template_idx = 3; % 1
%         % version_nr = 2; 
%             % ver = 1   -> CCV P1-T1 (template_idx=1); CCV P1-T2 (template_idx=3)
%             % ver = 1.1 -> Same as 1 but subtracts the noise
%             % ver = 2   -> CCV P1-T1 and P1-T2 (ADD ALIGNMENT)
% end
                 
if (isempty(sim_work.templ_tar) == 1 || cfg_sim.template_every_trial == 1 )

    def.calculate_template = 1;
    % generate template if not existing
    switch template_script
        case 'model_template'
            [templ_tar,templ_ref,cfg_sim] = model_template(cfg_game,data_passation,cfg_sim); 
    end
    sim_work.templ_tar = templ_tar;
    sim_work.templ_ref = templ_ref;
    
    def.calculate_template = 0;
    
end

str_stim = [];
eval(['str_stim = ' cfg_game.experiment '_user(cfg_game,data_passation);']); % calls user-function of the experiment
tuser = str_stim.tuser;

type_processing = cfg_game.experiment;
 
fs = cfg_game.fs;
modelname = cfg_sim.modelname;
[modelpars,subfs,params_extra] = model_params(modelname,fs);
if ~isempty(params_extra)
    in_var = params_extra.in_var;
else
    in_var = 0;
    disp('No internal noise being used...')
end

switch type_processing
    case 'pianoinnoise'
        % see my codes from TU/e
        
    case 'pianoinnoise_v2'
        % see my codes from TU/e
        
    otherwise % normal calculation
        try
            [ir_signal,xx,mfc] = feval(modelname,tuser,modelpars{:}); % xx is unused
        catch
            ir_signal = feval(modelname,tuser,modelpars{:});
        end
        [sim_work.current_signal, xx, info] = Ensure_intrep_is_numeric(ir_signal); % xx is unused
end

switch type_decision
    case 'relanoiborra2019_decision'
        % current_signal = Ensure_intrep_is_numeric_set_back(sim_work.current_signal,info);
        templ_tar      = Ensure_intrep_is_numeric_set_back(sim_work.templ_tar,info);
        templ_ref      = Ensure_intrep_is_numeric_set_back(sim_work.templ_ref,info);
        
        decision_with_tar = relanoiborra2019_backend2Dst_ENS_version(ir_signal,templ_tar,subfs,mfc);
        decision_with_ref = relanoiborra2019_backend2Dst_ENS_version(ir_signal,templ_ref,subfs,mfc);
        
        if decision_with_tar.dfinal >= decision_with_ref.dfinal
            response = 2;
        else
            response = 1;
        end
        sim_work.thres_for_bias = 0; % thres_for_bias;
        sim_work.decision_var_mue2choose(data_passation.i_current,:) = [decision_with_tar.dfinal decision_with_ref.dfinal];
        
    case 'optimal_detector'
        mue_tar = optimal_detector(sim_work.current_signal,sim_work.templ_tar,subfs);
        mue_ref = optimal_detector(sim_work.current_signal,sim_work.templ_ref,subfs);
        
        %%% Cross-correlation with different time lags:
        % mu1    = il_get_optimal_detector(sim_work.current_signal,sim_work.templ_tar,subfs); 
        % mu1rev = il_get_optimal_detector(sim_work.templ_tar,sim_work.current_signal,subfs); 
        % mu1 = max([mu1 mu1rev]);
        % 
        % mu2    = il_get_optimal_detector(sim_work.current_signal,sim_work.templ_ref,subfs);
        % mu2rev = il_get_optimal_detector(sim_work.templ_ref,sim_work.current_signal,subfs);
        % mu2 = max([mu2 mu2rev]);
        
        mue2choose_nonoise = [mue_ref mue_tar];
        if in_var ~= 0
            sigma = sqrt(in_var);
            int_noise = normrnd(0,sigma,size(mue2choose_nonoise));
            mue2choose = mue2choose_nonoise+int_noise;
        else
            % if no internal noise:
            mue2choose = mue2choose_nonoise;
        end
                
        version_decision = 2; % as in previous research (osses2021)
        
        switch version_decision
            case 1
                [xx,response]= max(mue2choose);
            case 2
                diff_value = mue2choose(2)-mue2choose(1);
                
                if diff_value >= thres_for_bias
                    response = 2;
                else
                    response = 1;
                end
                sim_work.thres_for_bias = thres_for_bias;
        end
        sim_work.decision_var_mue2choose(data_passation.i_current,:) = mue2choose;
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function mue = il_get_optimal_detector(aktsig,template,subfs)
% 
% res = floor(1e-3*subfs);
% sample_max = 50*res; % from +/- 100 ms to +/- 50 ms, changed on 24/11/2017 at 15:52
% 
% N = size(template,1);
% 
% j = 1;
% for i = 1:res:sample_max
%     mue(j) = 1/subfs*sum(aktsig(1:N-i+1).*template(i:N));
%     j = j + 1;
% end
