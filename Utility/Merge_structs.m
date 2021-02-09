function cfg_out = Merge_structs(cfg_in,cfg_out)
% function cfg_out = Merge_structs(cfg_in,cfg_out)
%
% Author: Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fn = fieldnames(cfg_in);
for i = 1:length(fn)
   cfg_out.(fn{i}) = cfg_in.(fn{i}); % adding recently created fields to cfg_out
end