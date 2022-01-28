function [bAre_stored_locally, dir_subject] = publ_varnet2013_FrontHumNeurosci_0_checkdata
% function [bAre_stored_locally, dir_subject] = publ_varnet2013_FrontHumNeurosci_0_checkdata
%
% This script checks whether the data related to varnet2013 (Varnet, Knoblauch, 
%     Meunier, and Hoen, 2013, Front. Hum. Neurosci.) are located locally and 
%     otherwise, it requests the user to download the data.
%
% Author: Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dir_subject = [];

dir_where = fastACI_dir_data;
experiment = 'speechACI_varnet2013';

dir_where = [dir_where experiment filesep];
dirs2check = {'varnet2013_S01','Sujet_Leo_S1','Sujet_Léo_S1'};

bAre_stored_locally = ones(size(dirs2check));

for i = 1:length(dirs2check)
    
    dir2check = [dir_where dirs2check{i} filesep];
    if ~exist(dir2check,'dir')
        bAre_stored_locally(i) = 0;
    else
        bAre_stored_locally(i) = 1;
    end
    
end

if sum(bAre_stored_locally) == 0
    fprintf('%s: directory %s not found on disk\n',upper(mfilename),dir2check);
    fprintf('\t Next steps to perform: \n');
    % fprintf('\t 1. Go to https://doi.org/10.5281/zenodo.5483835 (reference osses2021c_data) and download the two zip files\n');
    fprintf('\t 1. Contact the fastACI team to request the folder ''Sujet_Leo_S1'' and copy the zip file to your disk.')
    fprintf('\t 2. Uncompress the zip file and locate it into %s\n',dir_where);
    fprintf('\t 3. Re-run the current script until you don''t get errors \n');

    % if ~exist(dir_where,'dir')
    %     mkdir(dir_where)
    % end
else
    idx = find(bAre_stored_locally==1,1,'first');
    dir_subject = [dir_where dirs2check{idx} filesep];

    fprintf('%s: directory %s successfully found on disk\n',upper(mfilename),dir2check);
    fprintf('\t If you wish, you can now run either publ_osses2021c_DAGA_1_sim or publ_osses2021c_DAGA_2_figs... \n');
end