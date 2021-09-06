function load_name = Get_savenames(dir_results, exp2filter, dir_results_completed)
% function load_name = Get_savenames(dir_results, exp2filter, dir_results_completed)
%
% Programmed by Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ListSavegame = dir([dir_results exp2filter]);
index_savegame = 0;
bytes = 0;

index_all = 1:length(ListSavegame); % index of all the files that were found
for j=index_all
    if ListSavegame(j).bytes > bytes
        index_savegame=j;
        bytes = ListSavegame(j).bytes; % looks for the largest MAT file
    end
end
load_name = [dir_results ListSavegame(index_savegame).name];

% --- Now we will move the old (completed sessions) to the 
%     subject's folder.
index_all(index_savegame) = [];
if ~isempty(index_all)
    for j=1:length(index_all)
        movefile([dir_results           ListSavegame(index_all(j)).name], ... % src
                 [dir_results_completed ListSavegame(index_all(j)).name]);
    end
end