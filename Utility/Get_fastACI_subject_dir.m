function [dir_res, bExist] = Get_fastACI_subject_dir(experiment,Subject_ID)
% function [subjectname,dirname_new] = Get_fastACI_subject_dir(experiment,Subject_ID)
%
% Programmed by Alejandro Osses, ENS, France 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dir_res = [fastACI_dir_data experiment filesep Subject_ID filesep];

bExist = exist(dir_res,'dir');
if ~bExist
    mkdir(dir_res);
    bExist = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end