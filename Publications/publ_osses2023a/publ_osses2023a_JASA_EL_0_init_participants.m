function publ_osses2023a_JASA_EL_0_init_participants
% function publ_osses2023a_JASA_EL_0_init_participants
%
% Author: Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc

bCopy_from_crea = 1;

Subjects = {'S001','S002','S003','S004','S005','S006','S007','S008','S009','S010','S011','S012','S013','S014','S015','S016', ... % LAMI
            'S021','S022','S023','S024','S025','S026','S027','S028','S029','S030','S031','S032','S033','S034','S035','S080','S081','S082', ...  % LAPEL
            'S040','S041','S042','S060','S061', ... % LACROCH
            'S050','S051','S052','S053','S054', ... % LALARM
            'S070','S071','S072'}; % LAMI_SHIFTED 

if bCopy_from_crea

    for i = 1:length(Subjects)
        dir_crea = [fastACI_basepath 'Publications' filesep 'publ_osses2023a' filesep 'data_' Subjects{i} filesep '0-init' filesep];
        file_crea = Get_filenames(dir_crea,'cfgcrea*.mat');
        
        if ~isempty(file_crea)
            fprintf('%s.m: Generating the waveforms for participant %s\n',mfilename,Subjects{i});
            fprintf('\tCreation file: %s\n',file_crea{1});
            Generate_stimuli_from_cfg_crea([dir_crea file_crea{1}]);
        end
        
    end
    
end