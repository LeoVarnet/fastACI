% Script: speechACI_Logatome_instructions.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

time_per_trial = 2; % seconds per trial
dur_estimated = cfg_game.sessionsN*time_per_trial/60; 
switch cfg_game.Language
    case 'EN'
        % Text by Alejandro, based on the French text by Leo (below):
        fprintf('\t You will be presented with words in a background noise. Your task is to indicate on\n');
        fprintf('\t each trial whether if the word contains the syllable %s (enter 1) or \n',cfg_game.response_names{1});
        fprintf('\t contains the syllable %s (enter 2). \n',cfg_game.response_names{2});
        fprintf('\t The whole experiment has a total of %.0f sounds, divided into %.0f sessions of\n',cfg_game.N,cfg_game.N/cfg_game.sessionsN);
        fprintf('\t approximately %.0f minutes. Please note that you can take a break at any moment.\n',dur_estimated);
        fprintf('\t The level of the voice (the ''volume'') will be automatically adjusted according \n');
        fprintf('\t to your performance throughout the experiment. Therefore, it is normal that \n');
        fprintf('\t the words will rapidly become difficult to understand. If you can''t understand the \n');
        fprintf('\t words you can answer at random.\n');
        fprintf('\t The experiment includes the same number of %s and %s words. This means that ideally\n',cfg_game.response_names{1},cfg_game.response_names{2});
        fprintf('\t you should choose %s approximately 50%% of the times. In order to help you, \n',cfg_game.response_names{1});
        fprintf('\t your percentage of resonses "%s" will be shown if your bias becomes too important.\n',cfg_game.response_names{1});
        fprintf('\t After every response, the correct response  will be shown. If you can''t hear\n');
        fprintf('\t anymore which of the two words was played for a long period, \n');
        fprintf('\t maybe it is time to have a break!\n');

    case 'FR'
        % Text by Leo:
        fprintf('\t Vous allez entendre des sons dans un bruit de fond. Votre t\342che consiste \340\n');
        fprintf('\t indiquer apr\350s chaque \351coute si le son contenait un bip ("%s", appuyez sur 1) ou \n',cfg_game.response_names{1});
        fprintf('\t non ("%s", appuyez sur 2).\n',cfg_game.response_names{2});
        fprintf('\t L''exp\351rience comporte %.0f \351coutes au total, et elle est divis\351e en %.0f sessions de \n',cfg_game.N,cfg_game.N/cfg_game.sessionsN); 
        fprintf('\t %.0f minutes environ. Vous pouvez prendre une pause en cours de session si n\351cessaire.\n',dur_estimated);
        fprintf('\t L''exp\351rience comporte le m\352me nombre de "%s" et de "%s" au total. Dans l''id\351al, vous \n',cfg_game.response_names{1},cfg_game.response_names{2});
        fprintf('\t devriez donc r\351pondre "%s" 50 %% du temps.',cfg_game.response_names{1});
        fprintf('\t Pour vous aider, un son de r\351f\351rence contenant un ton clairement audible sera jou\351 toutes les 10 \351coutes')
       % fprintf('\t ponses "%s" et "%s" sera affich\351 si votre biais devient trop important.\n', cfg_game.response_names{1},cfg_game.response_names{2});
       % fprintf('\t Apr\350s chaque tentative, la r\351ponse correcte vous sera indiqu\351e. Si vous n''entendez \n');
      %  fprintf('\t plus l''un ou l''autre des deux sons pendant une longue p\351riode, il est peut-\352tre \n');
      %  fprintf('\t temps de faire une pause...\n')
end 