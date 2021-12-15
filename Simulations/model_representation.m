function [ir_reference,params] = model_representation(referencestim,modelname,modelpars)
% function [ir_reference,params] = model_representation(referencestim,modelname,modelpars)
%
% 1. Description:
%       CASPREPRESENTATION  Generate an internal representation to be used
%       by the optimal detector
%
%  CASPTEMPLATE(target,reference,modelname,modelpars) generates the template
%  needed for the optimal detector. CASPTEMPLATE will run the model specified 
%  by modelname on the signals stored in target and reference and generate 
%  the template from this.
%
%  If target or reference is a matrix, each column will be considered a
%  signal, and averaging will be done. This is usefull for stochastic
%  signals.
%
%   Url: http://amtoolbox.sourceforge.net/doc/modelstages/casptemplate.php
%
% 2. Stand-alone example:
%   fs          = 44100; % sampling frequency of the waveforms insig1 and insig2supra
%   target      = insig1;
%   reference   = insig2supra;
%   [template,ir_reference] = casptemplate(target,reference,'dau1996preproc',{fs});
% 
% Based on casprepresentation.m
%
% Copyright (C) 2009-2014 Peter L. Soendergaard and Piotr Majdak.
% This file is part of AMToolbox version 0.9.5
%
% Last update by Alejandro Osses: 16/12/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<3
    modelpars={};
end;

nreferences = size(referencestim,2);

%% ----- Compute average internal representation of the references
switch modelname
    case {'dau1996_preproc','dau1996a_preproc'}
        [ir_reference,fc,subfs] = feval(modelname,referencestim(:,1),modelpars{:});
    case {'dau1997','king2019','osses2021','osses2022a','dau1997_preproc','jepsen2008_preproc'}
        % [ir_reference,fc,fcm,subfs] = feval(modelname,referencestim(:,1),modelpars{:});
        [ir_reference,fc,fcm] = feval(modelname,referencestim(:,1),modelpars{:});
        subfs = Get_field_from_cell(modelpars,'subfs');
        params.fcm = fcm;
    case {'verhulst2018_preproc','verhulst2018debug_preproc'}
        [ir_reference,fc,fcm,subfs] = feval(modelname,referencestim(:,1),modelpars{:});
        params.fcm = fcm;
    otherwise
        error('Add models to this list');
end

for ii=2:nreferences
  ir_reference = ir_reference + feval(modelname,referencestim(:,ii),modelpars{:});
end;

if ~iscell(ir_reference)
    ir_reference = ir_reference/nreferences;
else
    ir_reference = il_cell_divide(ir_reference,nreferences);
end

%OLDFORMAT

if nargout >= 2
    params.fc = fc;
    try
        params.subfs = subfs;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function inoutcell = il_cell_divide(inoutcell,num)

for i = 1:length(inoutcell)
    
    inoutcell{i} = inoutcell{i}/num;
    
end
