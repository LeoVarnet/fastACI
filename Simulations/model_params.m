function [modelpars,subfs,params_extra] = model_params(modelname,fs,keyvals)
% function [modelpars,subfs,params_extra] = model_params(modelname,fs,keyvals)
%
% modelname = 'osses2021';
% fs = 44100;
% [modelpars,subfs] = model_params(modelname,fs);
%
% See also model_getparams.m (fastACI_sim repository):
% modelname = 'dau1997';
% fs = 44100;
% def = [];
% def.samplerate = fs;
% [modelpars,subfs] = model_getparams(modelname,def);
%
% Programmed by Alejandro Osses TU/e Eindhoven 2014-2018
% Extensions: 27/02/2023 - limiting 'fhigh' to fs/2 if fs is low enough
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    modelname = ''; % then standard parameters are loaded
end

params_extra = [];
 
% global simdef;
% fs    = def.samplerate;
modelpars = [];
modelpars{1} = fs;

subfs = [];

switch modelname
    case {'dau1997','dau1997_preproc'}
        
        subfs = 16000;
        % Custom parameters (comment to use the default parameters):
        % %%% Parameters up to AMT 0.11:
        % modelpars{end+1} = 'outerear';
        % modelpars{end+1} = 'middleear';
        % modelpars{end+1} = 'ihc_breebaart2001'; % 770-Hz LPF
        % modelpars{end+1} = 'adt_osses2020'; % overshoot limitation of 5 % OLD NAME: 'adt_dau_custom'
        % modelpars{end+1} = 'mfb_jepsen2008'; in_std = 0; in_var = in_std^2;
        % modelpars{end+1} = 'subfs'; modelpars{end+1} = subfs;
        
        %%% Parameters AMT 1.0:
        modelpars{end+1} = 'basef';
        modelpars{end+1} = [];
        
        modelpars{end+1} = 'no_outerear';
        modelpars{end+1} = 'no_middleear';
        modelpars{end+1} = 'ihc_dau1996'; % 1000-Hz LPF
        modelpars{end+1} = 'adt_dau1997'; % overshoot limitation of 10
        modelpars{end+1} = 'mfb_dau1997'; in_std = 0; in_var = in_std^2;
                
        params_extra.in_var = in_var;
        
    case 'king2019_old'
    % 
    %     pars = king2019_cfg(keyvals);
    %     for i = 1:length(pars.modelpars)
    %         modelpars{end+1} = pars.modelpars{i};
    %     end
    % 
    %     dBFS = 100;
    %     basef = 1000; % Hz
    %     flow = 80;
    %     fhigh = 8000;
    %     % modbank_Nmod = 5;
    %     mflow  =   2; % Hz, modbank_fmin
    %     mfhigh = 150;
    % 
    %     modelpars{end+1} = 'no_debug'; % flag
    %     modelpars{end+1} = 'phase_insens_hilbert'; % 'no_phase_insens'; % flag
    %     % default values will be loaded for king2019:
    %     modelpars(end+1:end+2) = {'compression_n',0.3};
    %     modelpars(end+1:end+2) = {'basef', basef}; % keyval
    %     modelpars(end+1:end+2) = {'flow' ,  flow}; % keyval
    %     modelpars(end+1:end+2) = {'fhigh', fhigh}; % keyval
    %     modelpars(end+1:end+2) = {'mflow' , mflow}; % keyval
    %     modelpars(end+1:end+2) = {'mfhigh', mfhigh}; % keyval
    %     % modelpars(end+1:end+2) = {'modbank_Nmod', modbank_Nmod}; % keyval
    %     modelpars(end+1:end+2) = {'dboffset', dBFS};
        
        disp('')
    case {'maxwell2020','maxwell2020_debug'}
        pars = maxwell2020_cfg;
        for i = 1:length(pars.modelpars)
            modelpars{end+1} = pars.modelpars{i};
        end
        subfs = fs;
        
    case {'relanoiborra2019','relanoiborra2019_preproc_debug'}
        pars = relanoiborra2019_cfg(keyvals);
        for i = 1:length(pars.modelpars)
            modelpars{end+1} = pars.modelpars{i};
        end
        subfs = pars.subfs;
        params_extra.in_var = pars.in_var;
        
    case 'osses2021'
        pars = osses2021_cfg(keyvals);
        for i = 1:length(pars.modelpars)
            modelpars{end+1} = pars.modelpars{i};
        end
        subfs = pars.subfs;
        if isfield(pars,'in_var')
            params_extra.in_var = pars.in_var;
        else
            params_extra.in_var = 0;
        end
        
    case 'osses2022a'
        pars = osses2022a_cfg(keyvals);
        for i = 1:length(pars.modelpars)
            modelpars{end+1} = pars.modelpars{i};
        end
        if isfield(pars,'subfs')
            subfs = pars.subfs;
        else
            subfs = [];
        end
        params_extra.in_var = pars.in_var;
        
    otherwise
        pars = []; % to be obtained in the next line...
        exp2eval = sprintf('pars = %s_cfg(keyvals);',modelname);
        eval(exp2eval); % evaluation of the *_cfg.m file
        
        if ~isfield(pars,'modelpars')
            fprintf('%s.m: The script %s_cfg.m does not provide the variable ''modelpars''. The default model configuration will be used\b',mfilename, modelname);
        else
            for i = 1:length(pars.modelpars)
                modelpars{end+1} = pars.modelpars{i};
            end
        end
end

if subfs > fs
    fprintf('\t%s.m: subfs is greater than fs, setting subfs to fs\n',mfilename);
    subfs = fs;
elseif isempty(subfs)
    subfs = fs; % assuming that there is no resampling in the auditory model
end
    

modelpars(end+1:end+2) = {'subfs',subfs};
if fs < 16000
    % Setting fhigh to fs/2, knowing that fs/2 should be by default greater than 8000 Hz
    modelpars(end+1:end+2) = {'fhigh',round(fs/2)}; 
end
