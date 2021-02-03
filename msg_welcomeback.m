%% Welcome back message

clc
fprintf('\n   *** WELCOME BACK! *** \n\n');
fprintf(['   Loading ' cfg_game.load_name '...\n']);
fprintf(['   Progress: stimulus # ' num2str(i-1) ' of ' num2str(cfg_game.N) '.\n\n']);
fprintf('   Press any key to start.\n')
pause
clc