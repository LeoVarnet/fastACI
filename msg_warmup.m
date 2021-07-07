%% Warmup + instructions message
clc
msg_instructions;
fprintf('\n'); 

switch cfg_game.Language
    case 'EN'
        fprintf('\t The experiment starts with a warm-up phase. Contrary to the test phase, here you will be\n');
        fprintf('\t allowed to replay the stimulus and the target tones, and the correct answer will be\n');
        fprintf('\t indicated after each trial.\n\n');
        fprintf('\t Press any key to start.\n');
    
    case 'FR'
        fprintf('\t L''experience debute par une phase d''echauffement. Contrairement a la vraie phase de \n');
        fprintf('\t test, vous aurez ici la possibilite de rejouer le son precedent ou les %s et %s sans\n',cfg_game.response_names{1},cfg_game.response_names{2});
        fprintf('\t bruit de fond.\n\n');
        fprintf('\t Appuyez sur n''importe quelle touche pour commencer.\n');
end        
pause
clc