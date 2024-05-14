function [] = results_ahumada1975(participantname)

expename = 'replication_ahumada1975';
dir_data = fastACI_dir_data();
dir_savegame = [dir_data expename filesep participantname filesep 'Results' filesep];
D = dir([dir_savegame 'savegame*.mat']);
load([dir_savegame D(1).name])

percentcorrect = mean(data_passation.is_correct);
bias = mean(data_passation.n_responses==2);
nbr_tones = sum(data_passation.n_targets==1);
nbr_noises = sum(data_passation.n_targets==2);
nbr_hits = sum(data_passation.n_targets==1 & data_passation.is_correct==1);
nbr_miss = sum(data_passation.n_targets==1 & data_passation.is_correct==0);
nbr_falsealarms = sum(data_passation.n_targets==2 & data_passation.is_correct==0);
nbr_correctrejections = sum(data_passation.n_targets==2 & data_passation.is_correct==1);
nbr_trials = length(data_passation.is_correct);

fprintf(['\n*** R\351sultats participant ' participantname ' ***\n'])
fprintf(['\nPourcentage de r\351ponses correctes : ' num2str(100*percentcorrect) ' %%.\n'])
fprintf(['Biais en faveur de bruit seul : ' num2str(100*bias) ' %%.\n'])
fprintf(['\nNombre total d''\351coutes : ' num2str(nbr_trials) '.\n'])
fprintf(['- Nombre de stimuli avec un ton : ' num2str(nbr_tones) '.\n'])
fprintf(['\tNombre de tons correctement d\351tect\351s : ' num2str(nbr_hits) '.\n'])
fprintf(['\tNombre de tons manqu\351s : ' num2str(nbr_miss) '.\n'])
fprintf(['- Nombre de stimuli bruit seul : ' num2str(nbr_noises) '.\n'])
fprintf(['\tNombre de bruits seuls correctement d\351tect\351s : ' num2str(nbr_correctrejections) '.\n'])
fprintf(['\tNombre de bruits seuls manqu\351s : ' num2str(nbr_falsealarms) '.\n'])
end