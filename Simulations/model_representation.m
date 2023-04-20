function [ir_reference,params] = model_representation(referencestim,modelname,modelpars)
% function [ir_reference,params] = model_representation(referencestim,modelname,modelpars)
%
% 1. Description:
% MODEL_REPRESENTATION Generates an internal representation to be used by
%       the optimal detector
%
%  If target or reference is a matrix, each column will be considered a
%  signal, and averaging will be done. This is usefull for stochastic
%  signals.
%
% Author: Alejandro Osses (16/12/2019), based on casprepresentation.m from AMT 0.9.5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<3
    modelpars={};
end

nreferences = size(referencestim,2);

%% ----- Compute average internal representation of the references
switch modelname
    case {'dau1996_preproc','dau1996a_preproc'}
        [ir_reference,fc,subfs] = feval(modelname,referencestim(:,1),modelpars{:});
    case {'dau1997','king2019','osses2021','osses2022a','relanoiborra2019_preproc_debug', ...
          'maxwell2020','maxwell2020_debug'} % 'dau1997_preproc','jepsen2008_preproc'}
        % [ir_reference,fc,fcm,subfs] = feval(modelname,referencestim(:,1),modelpars{:});
        [ir_reference,fc,fcm] = feval(modelname,referencestim(:,1),modelpars{:});
        subfs = Get_field_from_cell(modelpars,'subfs');
        params.fcm = fcm;
    otherwise
        ir_reference = feval(modelname,referencestim(:,1),modelpars{:});
        subfs = Get_field_from_cell(modelpars,'subfs');
        if isempty(subfs)
            if max( size(ir_reference) ) >= referencestim(:,1)
                subfs = modelpars{1}; % subfs equal to sampling frequency
            else
                error('Validate here the automatic assessment of subfs');
            end
        end
        fc = [];
        
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
