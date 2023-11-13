function cfg_inout = speechACI_Audika_legacy(cfg_inout,bShow)
% function cfg_crea = speechACI_Audika_legacy(cfg_crea,bShow)
%
% Function comparable to *_cfg.m functions from AFC toolbox
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    bShow = 1;
end

fields_updated = [];

old_name = 'dir_main';
new_name = 'dir_data_experiment';
date_changed = '24-06-2021';
if isfield(cfg_inout,old_name)
    exp2eval = sprintf('cfg_inout.%s=cfg_inout.%s;',new_name,old_name);
    eval(exp2eval);
    exp2eval = sprintf('cfg_inout = Remove_field(cfg_inout,''%s'');',old_name);
    eval(exp2eval);
    fields_updated{end+1} = sprintf('\t%s renamed to %s (old name not valid after %s)',old_name,new_name,date_changed);
end

old_name = 'dir_speech';
new_name = 'dir_target';
date_changed = '24-06-2021';
if isfield(cfg_inout,old_name)
    exp2eval = sprintf('cfg_inout.%s=cfg_inout.%s;',new_name,old_name);
    eval(exp2eval);
    exp2eval = sprintf('cfg_inout = Remove_field(cfg_inout,''%s'');',old_name);
    eval(exp2eval);
    fields_updated{end+1} = sprintf('\t%s renamed to %s (old name not valid after %s)',old_name,new_name,date_changed);
end

if bShow
    if ~isempty(fields_updated)
        Show_cell(fields_updated);
    end
end