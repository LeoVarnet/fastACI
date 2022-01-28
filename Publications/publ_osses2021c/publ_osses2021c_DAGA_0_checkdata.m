function bAre_stored_locally = publ_osses2021c_DAGA_0_checkdata
% function bAre_stored_locally = publ_osses2021c_DAGA_0_checkdata
%
% This script checks whether the data related to osses2021c (Osses and Varnet,
%     2022, DAGA) is located locally and otherwise, it requests the user to
%     download the data.
%
% Author: Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dir_where = fastACI_dir_data;
experiment = 'speechACI_varnet2013';

dir_where = [dir_where experiment filesep];
dirs2check = {'osses2021c_S01', ...
              'osses2021c_S02'};

bAre_stored_locally = ones(size(dirs2check));

for i = 1:length(dirs2check)
    
    dir2check = [dir_where dirs2check{i} filesep];
    if ~exist(dir2check,'dir')
        fprintf('%s: directory %s not found on disk\n',upper(mfilename),dir2check);
        fprintf('\t Next steps to perform: \n');
        fprintf('\t 1. Go to https://doi.org/10.5281/zenodo.5483835 (reference osses2021c_data) and download the two zip files\n');
        fprintf('\t 2. Uncompress the zip files and locate them into %s\n',dir_where);
        fprintf('\t 3. Re-run the current script until you don''t get errors \n');
        
        if ~exist(dir_where,'dir')
            mkdir(dir_where)
        end
        
        bAre_stored_locally = 0;
        return;
    else
        bAre_stored_locally = 1;
    end
    
end

if bAre_stored_locally == 1
    fprintf('%s: directory %s successfully found on disk\n',upper(mfilename),dir2check);
    fprintf('\t If you wish, you can now run either publ_osses2021c_DAGA_1_sim or publ_osses2021c_DAGA_2_figs... \n');
end
        