function [list_signals, list_target_signals] = Get_n_signals(cfg_game)
% function [list_signals, list_target_signals] = Get_n_signals(cfg_game)

list_signals = [];
list_target_signals = [];
for i=1:cfg_game.N_signal
    % for N_signal == 2 ==> first half equal to one, second half equal to two
    N_conditions = cfg_game.N_noise; % TODO: this is still not the most suitable naming...
    list_signals        = [list_signals i*ones(1,ceil(N_conditions))];%
    list_target_signals = [list_target_signals cfg_game.response_correct_target(i)*ones(1,ceil(N_conditions))];
end