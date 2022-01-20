function publ_osses2022b_preregistration_0_init_participants
% function publ_osses2022b_preregistration_0_init_participants
%
% Author: Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dir_here = fileparts(mfilename('fullpath')); dir_here = [dir_here filesep];
dir2add = [dir_here 'publ_osses2022b' filesep];
addpath(dir2add)

bInit_only = 1;
Subjects = {'S01','S02','S03','S04','S05','S06','S07','S08','S09','S10','S11','S12','S13','S14'};
N_subjects = length(Subjects);

for i = 1:N_subjects
    Subj = Subjects{i};
    f20220119_all_sessions_latin_square(Subj,bInit_only);
end
        
rmpath(dir2add);