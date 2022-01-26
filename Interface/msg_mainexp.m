%%  Main experiment + instructions message

clc
% display_instructions
msg_instructions(cfg_game);
switch cfg_game.Language
    case 'EN'
        fprintf('   Press any key to start.\n')
    case 'FR'
        fprintf('   Appuyez sur n''importe quelle touche pour commencer.\n')
end
pause
clc