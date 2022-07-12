function Save_crosspred_reduced(file_crosspred)
% function Save_crosspred_reduced(file_crosspred)
%
% 1. Description:
%       Read and reduces a cross-prediction file (Crosspred.mat).
%
% 2. Stand-along example (from an absolute path):
%       % Generating the reduced file:
%       file_crosspred = '/home/alejandro/Documents/Databases/data/fastACI_data_z_tmp/20220705_sim_Q1_osses2022a/Results-run-1-S01-v5-N-0010/ACI-osses2022a-speechACI_Logatome-white-nob-gt-l1glm/Crosspred.mat';
%       Save_crosspred_reduced(file_crosspred);
%
%       % Checking that the reduction was successful:
%       [crosspred,outs_performance] = Read_crosspred(file_crosspred);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bUpdated = 0; 
crosspred = [];
info_toolbox = [];
load(file_crosspred);

if isfield(crosspred,'yhat_train')
    crosspred = Remove_field(crosspred,'yhat_train');
    bUpdated = 1;
end
if isfield(crosspred,'yhat_test')
    crosspred = Remove_field(crosspred,'yhat_test');
    bUpdated = 1;
end

if bUpdated == 1
    fname = [file_crosspred(1:end-4) '_orig.mat'];
    movefile(file_crosspred,fname);
    
    fname = file_crosspred;
    save(fname,'crosspred','info_toolbox');
    fprintf('%s: Reduced Crosspred.mat file was created...\n',upper(mfilename));
else
    fprintf('%s: The Crosspred file is already in its reduced form, nothing was done...\n',upper(mfilename));
end


    
disp('')