function publ_osses2022b_quick(file_savegame)
% function publ_osses2022b_quick(file_savegame)
%
% Author: Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[cfg_game, data_passation, ListStim] = Convert_ACI_data_type(file_savegame);

expvar = data_passation.expvar;
expvar_description = cfg_game.expvar_description;
Me = [median(expvar),prctile(expvar,25),min(expvar),prctile(expvar,75),max(expvar)];
fprintf('\texpvar: Median (min, perc25, perc75, max)=%.1f (%.1f, %.1f, %.1f, %.1f) %s\n',Me,expvar_description);

bias_to_1 = length(find(data_passation.n_responses==1));
bias_to_2 = length(find(data_passation.n_responses==2));
N = length(data_passation.n_responses);
fprintf('\ttotal trials so far=%.0f.\n\t\tTimes 1 was chosen=%.0f (%.1f per cent)\n\t\tTimes 2 was chosen=%.0f (%.1f per cent)\n', ...
    N, ...
    bias_to_1,100*bias_to_1/N, ...
    bias_to_2,100*bias_to_2/N);

idxi = max(data_passation.resume_trial(end),1);
idxf = length(expvar);
idxs = idxi:idxf;

figure; 
plot(expvar(idxs),'bo-'); grid on
xlabel('Trial number');
ylabel(sprintf('expvar %s',expvar_description));
title(sprintf('Staircase of the last session (i=%.0f to i=%.f)',idxi,idxf));

disp('')
