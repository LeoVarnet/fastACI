function publ_osses2022b_preregistration_0_init_participants
% function publ_osses2022b_preregistration_0_init_participants
%
% Author: Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bCopy_from_crea = 1;
bCreate_from_scratch = ~bCopy_from_crea;

Subjects = {'S01','S02','S03','S04','S05','S06','S07','S08','S09','S10'}; % ,'S11','S12','S13','S14'};

if bCreate_from_scratch
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
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if bCopy_from_crea

    noise_types = {'white','sMPSv1p3','bumpv1p2_10dB'};
    noise_types_label = {'WN','MPS','BP'};

    %%% First, we check whether the participants' data are stored on disk:
    for i = 1:length(Subjects)
        for j = length(noise_types):-1:1 % reversed order (in case we need to remove one noise condition)
            outs = publ_osses2022b_utils(Subjects{i},noise_types{j},'Get_filenames');
            if outs.bGenerate_stimuli == 1
                warning('Re-generate waveforms first...');

                dir_crea = [fastACI_basepath 'Publications' filesep 'publ_osses2022b' filesep 'data_' Subjects{i} filesep '0-init' filesep];
                file_crea = Get_filenames(dir_crea,['cfgcrea*' noise_types{j} '.mat']);
                if length(file_crea) == 1
                    [~,bReproducible] = fastACI_experiment_init_from_cfg_crea([dir_crea file_crea{1}]);

                    if bReproducible == 0
                        warning('The sounds for Subject %s, cond=%s, don''t seem reproducible. Skippping this condition.',Subjects{i},noise_types{j});
                        noise_types(j) = [];
                        noise_types_label(j) = [];
                    end
                end

            end
        end
    end
    
end