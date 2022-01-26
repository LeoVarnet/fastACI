%% Warmup + instructions message
clc
msg_instructions(cfg_game);
fprintf('\n'); 

switch cfg_game.Language
    case 'EN'
        fprintf('\t The experiment starts with a warm-up phase. Contrary to the test phase, here you will be\n');
        fprintf('\t allowed to replay the stimulus and the target tones, and the correct answer will be\n');
        fprintf('\t indicated after each trial.\n\n');
        
        fprintf('\t Press any key to start.\n');
    
    case 'FR'
        fprintf('   L''exp\351rience d\351bute par une phase d''\351chauffement pour vous familiariser avec les sons\n')
        fprintf('   et la t\342che. Lorsque vous vous sentirez pr\352t, appuyez sur 6 pour passer \340 l''exp\351rience.\n\n');
        
        fprintf('   Appuyez sur n''importe quelle touche pour commencer.\n');
end        
pause
clc