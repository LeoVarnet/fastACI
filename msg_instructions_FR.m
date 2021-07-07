%% Display instructions
fprintf('\n   *** INSTRUCTIONS *** \n\n');

fprintf('   Vous allez entendre des sons de parole dans un bruit de fond. Votre tache consiste a\n')
fprintf(['   indiquer apres chaque ecoute si le son se terminait par un ' cfg_game.response_names{1} ' ou un ' cfg_game.response_names{2} '.\n'])
fprintf('   L''experience (%.0f ecoutes) est divisee en %.0f sessions de %.0f ecoutes chacune, mais\n',cfg_game.N,cfg_game.N/cfg_game.sessionsN,cfg_game.sessionsN)
fprintf('   vous pouvez prendre une pause en cours de session si necessaire.\n')
fprintf('   Le niveau de bruit s''ajustera automatiquement en fonction de vos performances tout\n')
fprintf('   au long de l''experience. Il est donc normal que la tache vous semble rapidement plus\n')
fprintf('   difficile. Si vous ne parvenez plus a distinguer le son de parole derriere le bruit, \n')
fprintf('   vous pouvez repondre au hasard.\n')
fprintf(['   L''experience comporte le meme nombre de ' cfg_game.response_names{1} ' et de ' cfg_game.response_names{2} '. Si vous n''entendez plus l''une\n'])
fprintf('   ou l''autre des deux sons pendant une longue periode, il est peut-etre temps de faire\n')
fprintf('   une pause...\n')
    
% fprintf('   You will be presented with tones in a background noise. Your task is to indicate on\n')
% fprintf('   each trial whether the tone was fluctuating ("modulated tone") or stable ("pure tone").\n')
% fprintf('   The whole experiment (%.0f trials) is divided into %.0f sessions of %.0f trials each, but\n',cfg_game.N,cfg_game.N/cfg_game.sessionsN,cfg_game.sessionsN)
% fprintf('   you are allowed to take a break at any moment.\n')
% fprintf('   The difficulty of the task will automatically adjust to your performance level throughout\n')
% fprintf('   the experiment. Therefore, the modulation to be detected should rapidly become thiner.\n')
% fprintf('   The experiment includes the same number of pure and modulated tones. If you cannot hear\n')
% fprintf('   one of the two categories anymore, maybe it is time for a break!\n')
    