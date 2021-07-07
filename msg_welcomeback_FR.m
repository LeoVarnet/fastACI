%% Welcome back message

clc
fprintf('\n   *** REPRISE DE L''EXPERIENCE ! *** \n\n');
fprintf(['   Chargement ' cfg_game.load_name '...\n']);
fprintf(['   Progression : ecoute numero ' num2str(i_current-1) ' sur un total de ' num2str(cfg_game.N) '.\n\n']);
fprintf('   Appuyez sur n''importe quelle touche pour commencer.\n')

% fprintf('\n   *** WELCOME BACK! *** \n\n');
% fprintf(['   Loading ' cfg_game.load_name '...\n']);
% fprintf(['   Progress: stimulus # ' num2str(i_current-1) ' of ' num2str(cfg_game.N) '.\n\n']);
% fprintf('   Press any key to start.\n')
pause
clc