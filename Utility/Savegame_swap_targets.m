function [savegame_file_out] = Savegame_swap_targets(savegame_file,varargin)
% function [savegame_file_out] = Savegame_swap_targets(savegame_file,varargin)
%
% Stand-alone example:
%    dir_where = '/home/alejandro/Documents/Databases/data/fastACI_data/modulationACI/S1/';
%    savegame_file = [dir_where 'savegame_2021_2_24_14_39.mat'];
%    savegame_file_out = Savegame_swap_targets(savegame_file);
%
% Programmed by Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% From argument function:
definput.import={'fastACI_getACI'}; % arg_fastACI_getACI.m
[flags,keyvals]  = ltfatarghelper({},definput,varargin);

var = load(savegame_file);
%% 1. Reading the experimental data (*.mat file):
[cfg_game, data_passation, ListStim] = Convert_ACI_data_type(savegame_file,keyvals);

tars_old = [1 2];
tars_new = [2 1];
switch cfg_game.experiment
    case 'AM_revcorr'
        if cfg_game.N_target ~= 2
            error('Only validated for experiments with two targets')
        end
        if length(cfg_game.response_correct_target) ~= 2
            error('Only validated for experiments with two targets')
        end
        
        cfg_game.target_names     = cfg_game.target_names(tars_new);
        cfg_game.response_names   = cfg_game.response_names(tars_new);
        cfg_game.n_targets_sorted = il_change_idxs(cfg_game.n_targets_sorted,tars_old,tars_new);
        
        % data_passation
        data_passation.n_targets  = il_change_idxs(data_passation.n_targets,tars_old,tars_new);
        data_passation.n_responses= il_change_idxs(data_passation.n_responses,tars_old,tars_new);
end
%%%
[path,name,ext]=fileparts(which(savegame_file));
if isempty(path)
    [path,name,ext]=fileparts( savegame_file );
end

if ~isempty(keyvals.dir_out)
    dir_out = keyvals.dir_out;
else
    dir_out = [path filesep]; 
end
fname_out = [name '-swap-tar'];

switch cfg_game.experiment
    case 'AM_revcorr'
        save([dir_out fname_out],'cfg_game','data_passation','ListStim');
end

disp('')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varinout = il_change_idxs(varinout,tars_old,tars_new)

idx1_old = find(varinout == tars_old(1));
idx2_old = find(varinout == tars_old(2));

varinout(idx1_old) = tars_new(1);
varinout(idx2_old) = tars_new(2);

disp('')