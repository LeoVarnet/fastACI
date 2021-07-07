%% Warmup + instructions message

clc
msg_instructions;
fprintf('\n'); 
fprintf('   L''experience debute par une phase d''echauffement. Contrairement a la vraie phase de \n')
fprintf(['   test, vous aurez ici la possibilite de rejouer le son precedent ou les ' cfg_game.response_names{1} ' et ' cfg_game.response_names{2} ' sans\n'])
fprintf('   bruit de fond.\n\n')
fprintf('   Appuyez sur n''importe quelle touche pour commencer.\n')

% fprintf('   The experiment starts with a warm-up phase. Contrary to the test phase, here you will be\n')
% fprintf('   allowed to replay the stimulus and the target tones, and the correct answer will be\n')
% fprintf('   indicated after each trial.\n\n')
% fprintf('   Press any key to start.\n')
pause
clc