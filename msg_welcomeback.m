%% Welcome back message
clc
switch cfg_game.Language
    case 'EN'
        fprintf('\n   *** WELCOME BACK! *** \n\n');
        fprintf('   Loading %s...\n',cfg_game.load_name);
        fprintf(['   Progress: stimulus # %.0f of %.0f.\n\n'],i_current-1,cfg_game.N);
        fprintf('   Press any key to start.\n')
    case 'FR'
        fprintf('\n   *** REPRISE DE L''EXPERIENCE ! *** \n\n');
        fprintf('   Chargement %s...\n',cfg_game.load_name);
        fprintf('   Progression : ecoute numero %.0f sur un total de %.0f.\n\n',i_current-1,cfg_game.N);
        fprintf('   Appuyez sur n''importe quelle touche pour commencer.\n')
end
pause
clc