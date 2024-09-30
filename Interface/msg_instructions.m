function msg_instructions(cfg_game)
% function msg_instructions(cfg_game)
% 
% Description: Shows the instructions of the corresponding experiment. If the
%    experiment contains a *_instructions.m script, then the text in that
%    script will be shown, otherwise, generic text will be displayed.
%
% This script is called from msg_mainexp.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Display instructions
fprintf('\n   *** INSTRUCTIONS *** \n\n');

script2run = [cfg_game.experiment '_instructions'];
if exist(script2run,'file')
    eval(script2run)
else
    %%% Generic message:
    switch cfg_game.Language
        case 'EN'
            fprintf('   You will be presented with tones in a background noise. Your task is to indicate on\n');
            fprintf('   each trial whether the tone was fluctuating ("modulated tone") or stable ("pure tone").\n');
            fprintf('   The whole experiment (%.0f trials) is divided into %.0f sessions of %.0f trials each, but\n',cfg_game.N_trials,cfg_game.N_trials/cfg_game.sessionsN,cfg_game.sessionsN);
            fprintf('   you are allowed to take a break at any moment.\n');
            fprintf('   The difficulty of the task will automatically adjust to your performance level throughout\n');
            fprintf('   the experiment. Therefore, the modulation to be detected should rapidly become thiner.\n');
            fprintf('   The experiment includes the same number of pure and modulated tones. If you cannot hear\n');
            fprintf('   one of the two categories anymore, maybe it is time for a break!\n');

        case 'FR'
            fprintf('   Vous allez entendre des sons de parole dans un bruit de fond. Votre t\342che consiste \340\n');
            fprintf('   indiquer apr\350s chaque \351coute si le son se terminait par un "%s" (appuyez sur 1) ou \n',cfg_game.response_names{1});
            fprintf('   un "%s" (appuyez sur 2).\n',cfg_game.response_names{2});
            fprintf('   L''exp\351rience comporte %.0f \351coutes au total, et elle est divis\351e en %.0f sessions de \n',cfg_game.N_trials,cfg_game.N_trials/cfg_game.sessionsN); 
            fprintf('   %.0f minutes environ. Vous pouvez prendre une pause en cours de session si n\351cessaire.\n',cfg_game.sessionsN*15/400);
            fprintf('   Le volume de la voix s''ajustera automatiquement en fonction de vos performances tout\n');
            fprintf('   au long de l''exp\351rience. Il est donc normal que la t\342che vous semble rapidement plus\n');
            fprintf('   difficile. Si vous ne parvenez plus \340 distinguer la voix derri\350re le bruit, vous\n');
            fprintf('   pouvez r\351pondre au hasard.\n')
            fprintf('   L''exp\351rience comporte le m\352me nombre de "%s" et de "%s" au total. Dans l''id\351al, vous \n',cfg_game.response_names{1},cfg_game.response_names{2});
            fprintf('   devriez donc r\351pondre "%s" 50 %% du temps. Pour vous aider, votre pourcentage de r\351-\n',cfg_game.response_names{1});
            fprintf('   ponses "%s" et "%s" sera affich\351 si votre biais devient trop important.\n', cfg_game.response_names{1},cfg_game.response_names{2});
            fprintf('   Apr\350s chaque tentative, la r\351ponse correcte vous sera indiqu\351e. Si vous n''entendez \n');
            fprintf('   plus l''un ou l''autre des deux sons pendant une longue p\351riode, il est peut-\352tre \n');
            fprintf('   temps de faire une pause...\n')
    end 
end
