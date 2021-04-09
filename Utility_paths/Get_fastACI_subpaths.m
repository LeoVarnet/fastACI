function misc = Get_fastACI_sim_subpaths(type) 
% function misc = Get_fastACI_sim_subpaths(type)
%   
% 1. Description:
% 
%   Get local subpaths.
% 
%   type:
%       'tb_Loudness_v12'
%       'tb_NMT' (Nucleus MATLAB Toolbox)
%   
%   Where:
%       db = database
%       ex = experiment
%       lx = LaTeX
%       tb = toolbox
%
% % Example 1: Subpaths under NMT toolbox:
%       misc_sub = Get_UGent_subpaths('tb_NMT');
%
% Programmed by Alejandro Osses, TUe-UGent-ENS, 2014-2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% path = Get_UGent_paths(type);
% 
% if strcmp(type,'Experiments_AFC')
%     misc.afcfrontends = [path 'afcfrontends' filesep];
%     misc.experiments  = [path 'experiments'  filesep];
%     misc.incrementdetection = [misc.experiments 'incrementdetection' filesep];  
%     misc.modulationincrement_ewert2004    = [misc.experiments 'modulationincrement_ewert2004'    filesep];
%     misc.modulationdetection_pieper2016   = [misc.experiments 'modulationdetection_pieper2016'   filesep];
%     misc.modulationdetection_verhulst2018 = [misc.experiments 'modulationdetection_verhulst2018' filesep];
%     misc.pianoinnoise = [misc.experiments 'pianoinnoise' filesep];
%     misc.tfs1         = [misc.experiments 'tfs1'         filesep];
%     misc.toneinnoise_FE = [misc.experiments 'toneinnoise_FE' filesep];
%      
% % elseif strcmp(type,'Psychoacoustics')
% %     
% %     % misc.CASP4AMT = [path 'CASP_Jepsen4AMT' filesep];
% %     misc.Loudness_CF  = [path 'Loudness_CF'  filesep];
% %     misc.Loudness_MGB = [path 'Loudness_MGB' filesep];
%      
% elseif strcmp(type,'tb_AM_AddOns')
%      
%     misc.defaults    = [path 'defaults'    filesep];
%     misc.demos       = [path 'demos'       filesep];
%     misc.experiments = [path 'experiments' filesep];
%     misc.general     = [path 'general'     filesep];
%     misc.mex         = [path 'mex'         filesep];
%     misc.modelstages = [path 'modelstages' filesep];
%     misc.models      = [path 'models'      filesep];
%      
% elseif strcmp(type,'tb_NMT')
%      
% %     misc.Filterbank        = [path 'Matlab' filesep 'Filterbank' filesep];
% %     misc.FrontEnd          = [path 'Matlab' filesep 'FrontEnd'   filesep];
% %     misc.FTM               = [path 'Matlab' filesep 'FTM'        filesep];
% %     misc.Implant           = [path 'Matlab' filesep 'Implant'    filesep];
% %     misc.LoudnessGrowth    = [path 'Matlab' filesep 'LoudnessGrowth' filesep];
%     misc.Processing        = [path 'Matlab' filesep 'Processing' filesep]; % Ensure_field is here
% %     misc.Sequence          = [path 'Matlab' filesep 'Sequence'   filesep]; % Collate_into_sequence.m
% %     misc.Strategy          = [path 'Matlab' filesep 'Strategy'   filesep]; 
%     misc.Utility           = [path 'Matlab' filesep 'Utility'    filesep]; % From_dB; To_dB; Matlab_version.m
%      
% elseif strcmp(type,'tb_Plot4papers_JU')
%     misc.Utility           = [path 'Utility'  filesep];
%      
% elseif strcmp(type,'Reports')
%     misc.Reports2021 = [path '2021' filesep];    
%     misc.Reports2020 = [path '2020' filesep];
%     misc.Reports2019 = [path '2019' filesep];
%     misc.Reports2018 = [path '2018' filesep];
%     
% elseif strcmp(type,'Verhulstetal_current')    
%     misc.AddOns = [path 'AddOns' filesep];
%     
% else
    misc = [];
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
