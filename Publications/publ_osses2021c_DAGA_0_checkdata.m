function bAre_stored_locally = publ_osses2021c_DAGA_0_checkdata
% function bAre_stored_locally = publ_osses2021c_DAGA_0_checkdata
%
% Author: Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dir_where = fastACI_dir_data;
experiment = 'speechACI_varnet2013';
dirs2check = {[experiment filesep 'osses2021c_S01' filesep], ...
              [experiment filesep 'osses2021c_S02' filesep]};

bAre_stored_locally = ones(size(dirs2check));

for i = 1:length(dirs2check)
    
    dir2check = [dir_where dirs2check{i}];
    if ~exist(dir2check,'dir')
        fprintf('%s: directory %s successfully not found on disk\n',upper(mfilename),dir2check);
        fprintf('\t Next steps to perform: \n');
        fprintf('\t 1. Go to https://doi.org/10.5281/zenodo.5483835 (reference osses2021c_data) and download the two zip files\n');
        fprintf('\t 2. Uncompress the zip files and locate them into %s\n',dir_where);
        fprintf('\t 3. Re-run the current script until you don''t get errors \n');
        
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
        