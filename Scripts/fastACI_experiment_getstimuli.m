function cfg_game = fastACI_experiment_getstimuli(cfg_game)
% function cfg_inout = speechACI_Logatome_init(cfg_game)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
experiment = cfg_game.experiment;

init_file = [experiment '_init'];
exp2eval = sprintf('%s(cfg_game);',init_file);
eval(exp2eval);

disp('')
%%% 