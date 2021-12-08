function [modelpars,subfs,params_extra] = model_params(modelname,fs)
% function [modelpars,subfs,params_extra] = model_params(modelname,fs)
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
        modelpars{end+1} = 'outerear';
        modelpars{end+1} = 'middleear';
        modelpars{end+1} = 'ihc_breebaart2001'; % 770-Hz LPF
        modelpars{end+1} = 'adt_osses2021'; % overshoot limitation of 5 % OLD NAME: 'adt_dau_custom'
        modelpars{end+1} = 'mfb_jepsen2008'; in_std = 0; in_var = in_std^2;
                
        params_extra.in_var = in_var;
        
    case 'king2019'
        
        if isunix
            dBFS = 100;
            basef = 1000; % Hz
            bIsLeo = 0;
            
        elseif iswindows
            dBFS = 100;
            basef = 1000; % Hz
            warning('Leo: update in this script the dBFS value...');
            
            bIsLeo = 1;
        end
        
        modelpars{end+1} = 'no_debug'; % flag
        if bIsLeo
            %%% Leo's parameters:
            modelpars{end+1} = 'no_phase_insens'; % flag
            modelpars(end+1:end+2) = {'basef',basef}; % keyval
            modelpars(end+1:end+2) = {'flow' ,basef};
            modelpars(end+1:end+2) = {'fhigh',basef};
            modelpars(end+1:end+2) = {'mflow' ,4};
            modelpars(end+1:end+2) = {'mfhigh',4};
        else
            % default values will be loaded for king2019:
            modelpars(end+1:end+2) = {'basef',basef}; % keyval
            modelpars(end+1:end+2) = {'flow' ,   80}; % keyval
            modelpars(end+1:end+2) = {'fhigh', 8000}; % keyval
            
            % modelpars(end+1:end+2) = {'modbank_Nmod', 5}; % keyval
        end
        modelpars(end+1:end+2) = {'dboffset',dBFS};
        
    case 'osses2021'
        subfs = 16000;
        modelpars(end+1) = {'LP_150_Hz_att'}; % as in Osses2021a, Fig. 14c
        pars = osses2021_cfg;
        params_extra.in_var = pars.in_var;
end

if isempty(subfs)
    subfs = Get_field_from_cell(modelpars,'subfs');
    
    eval(sprintf('var = arg_%s;',modelname)) % loading defaults from AMT

    if ~isempty(var.keyvals.subfs)
        subfs = var.keyvals.subfs;
        fprintf('%s: Using default ''subfs'' loaded from arg_%s\n',upper(mfilename),modelname);
    else
        subfs = fs;
    end

    modelpars(end+1:end+2) = {'subfs',subfs};
else
    modelpars(end+1:end+2) = {'subfs',subfs};
end
