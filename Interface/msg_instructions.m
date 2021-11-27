%% Display instructions
fprintf('\n   *** INSTRUCTIONS *** \n\n');
switch cfg_game.Language
    case 'EN'
        fprintf('   You will be presented with tones in a background noise. Your task is to indicate on\n');
        fprintf('   each trial whether the tone was fluctuating ("modulated tone") or stable ("pure tone").\n');
        fprintf('   The whole experiment (%.0f trials) is divided into %.0f sessions of %.0f trials each, but\n',cfg_game.N,cfg_game.N/cfg_game.sessionsN,cfg_game.sessionsN);
        fprintf('   you are allowed to take a break at any moment.\n');
        fprintf('   The difficulty of the task will automatically adjust to your performance level throughout\n');
        fprintf('   the experiment. Therefore, the modulation to be detected should rapidly become thiner.\n');
        fprintf('   The experiment includes the same number of pure and modulated tones. If you cannot hear\n');
        fprintf('   one of the two categories anymore, maybe it is time for a break!\n');
        
    case 'FR'
        fprintf('   Vous allez entendre des sons de parole dans un bruit de fond. Votre t\342che consiste \340\n');
        fprintf('   indiquer apr\350s chaque \351coute si le son se terminait par un "%s" (appuyez sur 1) ou \n',cfg_game.response_names{1});
        fprintf('   un "%s" (appuyez sur 2).\n',cfg_game.response_names{2});
        fprintf('   L''exp\351rience comporte %.0f \351coutes au total, et elle est divis\351e en %.0f sessions de \n',cfg_game.N,cfg_game.N/cfg_game.sessionsN); 
        fprintf('   %.0f minutes environ. Vous pouvez prendre une pause en cours de session si n\351cessaire.\n',cfg_game.sessionsN*15/400);
        fprintf('   Le volume de la voix s''ajustera automatiquement en fonction de vos performances tout\n');
        fprintf('   au long de l''exp\351rience. Il est donc normal que la t\342che vous semble rapidement plus\n');
        fprintf('   difficile. Si vous ne parvenez plus \340 distinguer les sons de parole derri\350re le bruit, \n');
        fprintf('   vous pouvez r\351pondre au hasard.\n')
        fprintf('   L''exp\351rience comporte le m\352me nombre de "%s" et de "%s" au total. Dans l''id\351al, vous \n',cfg_game.response_names{1},cfg_game.response_names{2});
        fprintf('   devriez donc r\351pondre "%s" 50 %% du temps. Pour vous aider, votre pourcentage de r\351ponses\n',cfg_game.response_names{1});
        fprintf('   "%s" sera affich\351 tout au long de l''exp\351rience, ainsi que les r\351ponses correctes.  Si\n',cfg_game.response_names{1});
        fprintf('   vous n''entendez plus l''un ou l''autre des deux sons pendant une longue p\351riode, il est\n');
        fprintf('   peut-\352tre temps de faire une pause...\n')
end 