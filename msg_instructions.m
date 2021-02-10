%% Display instructions
fprintf('\n   *** INSTRUCTIONS *** \n\n');
fprintf('   You will be presented with tones in a background noise. Your task is to indicate on\n')
fprintf('   each trial whether the tone was fluctuating ("modulated tone") or stable ("pure tone").\n')
fprintf('   The whole experiment (%.0f trials) is divided into %.0f sessions of %.0f trials each, but\n',cfg_game.N,cfg_game.N/cfg_game.sessionsN,cfg_game.sessionsN)
fprintf('   you are allowed to take a break at any moment.\n')
fprintf('   The difficulty of the task will automatically adjust to your performance level throughout\n')
fprintf('   the experiment. Therefore, the modulation to be detected should rapidly become thiner.\n')
fprintf('   The experiment includes the same number of pure and modulated tones. If you cannot hear\n')
fprintf('   one of the two categories anymore, maybe it is time for a break!\n')
    