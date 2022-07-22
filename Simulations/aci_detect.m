function [response,sim_work,cfg_sim] = aci_detect(cfg_game,data_passation,cfg_sim,sim_work,varargin) %keyvals)
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

definput.import={'fastACI_simulation_detect'}; % arg_fastACI_simulation_detect.m
[flags,keyvals]  = ltfatarghelper({},definput,varargin);

response = NaN; % initialisation

if nargin < 4
    sim_work = [];
end
sim_work = Ensure_field(sim_work,'templ_tar',[]);
sim_work = Ensure_field(sim_work,'templ_ref',[]);
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
        case {'osses2021','osses2022a'}
            cfg_sim.type_decision = 'optimal_detector';
        case 'relanoiborra2019'
            cfg_sim.type_decision = 'relanoiborra2019_decision';
        otherwise
            error('%s: Add the model and its default decision script to this list',upper(mfilename))
    end
end
type_decision = cfg_sim.type_decision; % 'relanoiborra2019_decision'; 'optimal_detector'; % type_processing;
sim_work.type_decision = type_decision;

if flags.do_bias_global || flags.do_no_bias_each_session
    if isfield(cfg_sim,'thres_for_bias')
        thres_for_bias = cfg_sim.thres_for_bias;
    else
        switch type_decision
            case 'optimal_detector'
                disp('Non calibrated thres_for_bias')
                thres_for_bias = 0;
        end
    end
end
if flags.do_no_bias_global || flags.do_bias_each_session
    idx_session = length(data_passation.resume_trial);
    if isfield(cfg_sim,'thres_for_bias_each_session')
        thres_for_bias = cfg_sim.thres_for_bias_each_session(idx_session);
        if (data_passation.i_current - data_passation.resume_trial(end) == 0) || data_passation.i_current == 1
            fprintf('%s: Using one bias value (thres_for_bias=%.4f) for each section of %.0f trials\n',upper(mfilename),thres_for_bias,cfg_game.sessionsN);
        end
    else
        error('%s: if do_bias_each_session, you need to specify the keyval ''thres_for_bias_each_session''.',upper(mfilename));
    end
end

if (isempty(sim_work.templ_tar) == 1 || cfg_sim.template_every_trial == 1 )

    % def.calculate_template = 1;
    % generate template if not existing
    if cfg_sim.template_every_trial == 1
        cfg_game_here = cfg_game;
        if exist([cfg_game.experiment '_online_user.m'],'file')
            cfg_game_here.experiment = [cfg_game.experiment '_online'];
        end
    else
        cfg_game_here = cfg_game;
    end
    switch template_script
        case 'model_template'
            [templ_tar,templ_ref,cfg_sim] = model_template(cfg_game_here,data_passation,cfg_sim,keyvals); 
            
        case 'model_template_update'
            error('%s: The template script ''model_template_update.m'' is only available in the private fastACI_sim repository',upper(mfilename));
    end
    sim_work.templ_tar = templ_tar;
    sim_work.templ_ref = templ_ref;
    
    % def.calculate_template = 0;
    if sim_work.bStore_template
        
        fname = fastACI_file_template(cfg_game.experiment_full,cfg_game.Subject_ID,cfg_sim.type_decision,keyvals);
        % fname = sprintf('%stemplate-%s-%s-trial-1',cfg_game.dir_results,cfg_game.Subject_ID,cfg_game.experiment_full);
        if exist(fname,'file')
            disp('A template was found on disk. Press ctrl+c to abort or press any button to continue (and overwrite) the previous template')
            pause
        end
        save(fname,'templ_tar','templ_ref','cfg_sim');
        
    end
    
end

str_stim = [];
eval(['str_stim = ' cfg_game.experiment '_user(cfg_game,data_passation);']); % calls user-function of the experiment
tuser = str_stim.tuser;

type_processing = cfg_game.experiment;
 
fs = cfg_game.fs;
modelname = cfg_sim.modelname;
[modelpars,subfs,params_extra] = model_params(modelname,fs,keyvals);
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
        N_Ch = size(tuser,2);
        if N_Ch == 1
            % Monaural input signal
            try
                [ir_signal,xx,mfc] = feval(modelname,tuser,modelpars{:}); % xx is unused
            catch
                ir_signal = feval(modelname,tuser,modelpars{:});
            end
            [sim_work.current_signal, xx, info] = Ensure_intrep_is_numeric(ir_signal); % xx is unused
        else
            % Monaural input signal
            sim_work.current_signal = [];
            for n = 1:N_Ch
                try
                    [ir_signal,xx,mfc] = feval(modelname,tuser(:,n),modelpars{:}); % xx is unused
                catch
                    ir_signal = feval(modelname,tuser(:,n),modelpars{:});
                end
                [ir_signal, xx, info] = Ensure_intrep_is_numeric(ir_signal); % xx is unused
                sim_work.current_signal = [sim_work.current_signal; ir_signal];
            end
            
        end
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
        
        if keyvals.maxtimelag_ms == 0
            % The value at lag=0
            mue_tar = optimal_detector(sim_work.current_signal,sim_work.templ_tar,subfs);
            mue_ref = optimal_detector(sim_work.current_signal,sim_work.templ_ref,subfs);

        else
            %%% Cross-correlation with different time lags:
            mue_forw  = il_get_optimal_detector(sim_work.current_signal,sim_work.templ_tar,subfs,keyvals.maxtimelag_ms); 
            mue_backw = il_get_optimal_detector(sim_work.templ_tar,sim_work.current_signal,subfs,keyvals.maxtimelag_ms); 
            ccf = [flip(mue_backw) mue_forw(2:end)]; % first element of mue_forw is at zero lag (removed, already in mue_backw)
            [mue_tar,idx_lag] = max(ccf);
            lag = -keyvals.maxtimelag_ms:1:keyvals.maxtimelag_ms; % so far only spaced at 1 ms;
            if lag(idx_lag) == 0
                % Nothing to do
            else
                fprintf('\t%s: A time lag different from 0 was chosen (lag=%.0f ms) for the target criteron\n',upper(mfilename),lag(idx_lag));
            end    
            
            mue_forw  = il_get_optimal_detector(sim_work.current_signal,sim_work.templ_ref,subfs,keyvals.maxtimelag_ms);
            mue_backw = il_get_optimal_detector(sim_work.templ_ref,sim_work.current_signal,subfs,keyvals.maxtimelag_ms);
            ccf = [flip(mue_backw) mue_forw(2:end)];
            [mue_ref,idx_lag] = max(ccf);
            if lag(idx_lag) == 0 % idx_lag == 1 || idx_lag == 1+length(mue_forw)
                % Nothing to do
            else
                fprintf('\t%s: A time lag different from 0 was chosen (lag=%.0f) for the reference criteron\n',upper(mfilename),lag(idx_lag));
            end
        end
        
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
function mue = il_get_optimal_detector(aktsig,template,subfs,maxtimelag_ms)
 
% Here, the second input signal ('template') is shifted backwards with respect
%   to the first input signal ('aktsig'), or, equivalently, that aktsig is 
%   shifted forward. In other words a positive time lag is applied to aktsig

res = floor(1e-3*subfs);
sample_max = maxtimelag_ms*res; % from +/- 100 ms to +/- 50 ms, changed on 24/11/2017 at 15:52

N = size(template,1);

j = 1;
for i = 1:res:sample_max+1
    % It starts always with zero lag:
    mue(j) = optimal_detector(aktsig(1:N-i+1),template(i:N),subfs);
    j = j + 1;
end
