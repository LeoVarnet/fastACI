% Script: segmentation_instructions.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

time_per_trial = 2.2; % seconds per trial
dur_estimated = cfg_game.sessionsN*time_per_trial/60; 
switch cfg_game.Language
    case 'EN'
        % Text by Alejandro, based on the French text by Leo (below):
        fprintf('\t You will be presented with sentences that can be either "c''est l''%s" or "c''est la %s".\n',cfg_game.response_names{1},cfg_game.response_names{2});
        fprintf('\t Your task is to indicate on each trial whether the order is "c''est l''%s" the "c''est la %s" (press 1) or \n',cfg_game.response_names{1},cfg_game.response_names{2});
        fprintf('\t "c''est la %s" the "c''est l''%s"  (press 2). \n',cfg_game.response_names{2},cfg_game.response_names{1});
        fprintf('\t The whole experiment has a total of %.0f sounds, divided into %.0f sessions of approximately\n',cfg_game.N,(cfg_game.N/cfg_game.intervalnum)/cfg_game.sessionsN);
        fprintf('\t %.0f minutes. Please note that you can take a break at any moment.\n',dur_estimated);
        fprintf('\t It is normal that the voice sometimes sounds unnatural. If the answer is not obvious,\n');
        fprintf('\t please make a guess or answer at random\n');

    case 'FR'
        % Text by Leo:
        fprintf('\t Vous allez entendre une voix prononcer la phrase "c''est l''%s" ou "c''est la %s".\n',cfg_game.response_names{1},cfg_game.response_names{2});
        fprintf('\t Votre t\342che consiste \340 indiquer apr\350s chaque \351coute si l''odre des sons était \n');
        fprintf('\t "c''est l''%s" puis "c''est la %s" (appuyez sur 1) ou "c''est la %s" puis "c''est l''%s" (appuyez sur 2).\n',cfg_game.response_names{1},cfg_game.response_names{2},cfg_game.response_names{2},cfg_game.response_names{1});
        fprintf('\t L''exp\351rience comporte %.0f \351coutes au total, et elle est divis\351e en %.0f sessions de \n',cfg_game.N,(cfg_game.N/cfg_game.intervalnum)/cfg_game.sessionsN); 
        fprintf('\t %.0f minutes environ. Vous pouvez prendre une pause en cours de session si n\351cessaire.\n',dur_estimated);
        fprintf('\t Il est normal que la t\342che vous semble parfois difficile. Si vous ne parvenez pas \340\n');
        fprintf('\t identifier le dernier mot, vous pouvez r\351pondre au hasard\n');
end 
