% Script: segmentation_instructions.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

time_per_trial = 2.2; % seconds per trial
dur_estimated = cfg_game.sessionsN*time_per_trial/60; 
switch cfg_game.Language
    case 'EN'
        % Text by Alejandro, based on the French text by Leo (below):
        fprintf('\t You will be presented with sentences that can be either "c''est l''%s" or "c''est la %s".\n',cfg_game.response_names{1},cfg_game.response_names{2});
        fprintf('\t Your task is to indicate on each trial whether if the last word is "%s" (press 1) or \n',cfg_game.response_names{1});
        fprintf('\t "%s" (press 2). \n',cfg_game.response_names{2});
        fprintf('\t The whole experiment has a total of %.0f sounds, divided into %.0f sessions of approximately\n',cfg_game.N,cfg_game.N/cfg_game.sessionsN);
        fprintf('\t %.0f minutes. Please note that you can take a break at any moment.\n',dur_estimated);
        fprintf('\t It is normal that the voice sometimes sounds unnatural. If the answer is not obvious,\n');
        fprintf('\t please make a guess or answer at random\n');
        fprintf('\t The experiment includes the same number of %s and %s words. This means that ideally \n',cfg_game.response_names{1},cfg_game.response_names{2});
        fprintf('\t you should choose %s approximately 50%% of the times. In order to help you, your percen-\n',cfg_game.response_names{1});
        fprintf('\t tage of resonses "%s" will be shown if your bias becomes too important.\n',cfg_game.response_names{1});
        fprintf('\t After every response, the correct response  will be shown. If you can''t hear anymore which\n');
        fprintf('\t of the two words was played for a long period, maybe it is time to have a break!\n');

    case 'FR'
        % Text by Leo:
        fprintf('\t Vous allez entendre une voix prononcer la phrase "c''est l''%s" ou "c''est la %s".\n',cfg_game.response_names{1},cfg_game.response_names{2});
        fprintf('\t Votre t\342che consiste \340 indiquer apr\350s chaque \351coute si la phrase se terminait par\n');
        fprintf('\t "l''%s" (appuyez sur 1) ou "la %s" (appuyez sur 2).\n',cfg_game.response_names{1},cfg_game.response_names{2});
        fprintf('\t L''exp\351rience comporte %.0f \351coutes au total, et elle est divis\351e en %.0f sessions de \n',cfg_game.N,cfg_game.N/cfg_game.sessionsN); 
        fprintf('\t %.0f minutes environ. Vous pouvez prendre une pause en cours de session si n\351cessaire.\n',dur_estimated);
        fprintf('\t Il est normal que la t\342che vous semble parfois difficile. Si vous ne parvenez pas \340\n');
        fprintf('\t identifier le dernier mot, vous pouvez r\351pondre au hasard\n');
        fprintf('\t L''exp\351rience comporte le m\352me nombre de "%s" et de "%s" au total. Dans l''id\351al, \n',cfg_game.response_names{1},cfg_game.response_names{2});
        fprintf('\t vous devriez donc r\351pondre "%s" 50 %% du temps. Pour vous aider, votre pourcentage\n',cfg_game.response_names{1});
        fprintf('\t de r\351ponses "%s" et "%s" sera affich\351 si votre biais devient trop important.\n', cfg_game.response_names{1},cfg_game.response_names{2});
        fprintf('\t Apr\350s chaque tentative, la r\351ponse correcte vous sera indiqu\351e. Si vous n''entendez \n');
        fprintf('\t plus l''un ou l''autre des deux mots pendant une longue p\351riode, il est peut-\352tre \n');
        fprintf('\t temps de faire une pause...\n')
end 
