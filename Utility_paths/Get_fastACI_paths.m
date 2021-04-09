function misc = Get_UGent_paths(type) 
% function misc = Get_UGent_paths(type)
%   
%   Get local paths as structured at the TU/e Dell laptop (since 01/05/2014).
%   If 'type' is specified, then 'misc' will be a character containing the
%   that specific path. See example 2 (below).
% 
%   types available:
%       'db_fastl2007'
%       'db_voice_of_dragon'
%       'db_speechmaterials' - location X-drive
%       'db_speechmaterials_local' - local back-up of speech materials
%       'db_voice_of_dragon'
%       'Music_acoustics'
%       'praat'     - location of Praat (terminal: whereis praat)
%       'tb_AM'     - Auditory Modelling Toolbox
%       'tb_NMT'    - Nucleus MATLAB Toolbox    
%   
%   Where:
%       db = database
%       ex = experiment
%       lx = LaTeX
%       tb = toolbox
% 
% % Example 1: to get all important directories
%       misc = Get_UGent_paths;
%
% % Example 2: to get directory where database of the voice of the dragon is:
%       misc = Get_UGent_paths('db_voice_of_dragon');
%
% Programmed by Alejandro Osses, Hearing Technology @ WAVES, UGent, Belgium, 2018
% Created on    : 05/10/2018
% Last update on: 05/10/2018
% Last use on   : 30/09/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Alejandro's paths

if ~isunix
    which_pc = 1; % Win7
else
    which_pc = 2; % Unix-based
end

switch which_pc
    case 1
        
%         % % Win7, TU/e:
%         misc.username           =  'aosses';
%         if info.bOnline
%             misc.MATLAB         = 'G:\MATLAB\';
%         else
%             misc.MATLAB         = 'D:\MATLAB_git\';
%         end
%         misc.Databases          = 'D:\Databases\';
%         misc.outputs            = 'D:\Output\';
%         
%         misc.SVN_KUL            = ['D:\SVN-KU-Leuven\alejandro' filesep];
%         misc.CIs                =  [misc.MATLAB 'CIs' filesep];
% 
%         misc.ex_APEX_results    = [misc.SVN_KUL 'Meas' filesep 'Experiments' filesep 'Results_XML' filesep]; % work at KUL
%         misc.praat              = 'C:\praat5376_win32\praatcon.exe';
        
    case 2 
            
        misc.username           = 'alejandro';
        home_path               = [filesep 'home' filesep misc.username filesep];
        
        misc.MATLAB             = [home_path 'Documents/MATLAB/MATLAB_UGent/'];
        % misc.Databases          = '~/Documents/Databases/';

        misc.outputs            = ['~/Documents/MATLAB/outputs' filesep];
        misc.outputs_full       = [home_path 'Documents/MATLAB/outputs' filesep];
        % misc.praat              =  '/usr/bin/praat'; % whereis praat
        % misc.praat_scripts      = ['~/Documents/Praat_svn/'];
    
end

misc.DSP                = [misc.MATLAB   'DSP'             filesep]; 
% misc.Piano              = [misc.DSP      'Piano_analysis'  filesep];
% misc.Experiments        = [misc.MATLAB   'Experiments'     filesep]; 
misc.Experiments_AFC    = [misc.MATLAB   'Experiments_AFC' filesep];
misc.Experiments_WAE    = [misc.MATLAB   'Experiments_WAE' filesep]; % added on 15/04/2020
% misc.F0_extraction      = [misc.MATLAB   'F0_extraction'   filesep];
% misc.Filterbank         = [misc.MATLAB   'Filterbank'      filesep];
% misc.language           = 'NL'; % other possibilities: 'EN'
% misc.Localisation       = [misc.MATLAB   'Localisation'    filesep];
% misc.praat_scripts      = [misc.MATLAB   'Praat'           filesep];
misc.Stim_generation    = [misc.MATLAB   'Stim_generation' filesep];
misc.Plotting           = [misc.MATLAB   'Plotting'        filesep];
% misc.Psychoacoustics    = [misc.MATLAB   'Psychoacoustics' filesep];
misc.Reports            = [misc.MATLAB   'Reports'         filesep];
% misc.Reports_KUL        = [misc.MATLAB   'Reports_KUL'     filesep];
% misc.Room_acoustics     = [misc.MATLAB   'Room_acoustics'  filesep];
% misc.Loudness_CF        = [misc.Psychoacoustics 'Loudness_CF'  filesep];
% misc.Loudness_MGB       = [misc.Psychoacoustics 'Loudness_MGB' filesep];
% misc.Roughness_Daniel   = [misc.Psychoacoustics 'Roughness_Daniel'   filesep];
% misc.Speech             = [misc.MATLAB   'Speech'          filesep];
% misc.ltass              = [misc.Speech   'LTASS'           filesep];
% misc.ICRA               = [misc.Speech   'ICRA'            filesep];
% misc.ICRA_Tobias        = [misc.Speech   'ICRA_Tobias'     filesep];
% misc.ICRA_Tobias_Tools  = [misc.ICRA_Tobias 'Tools'       filesep];
% misc.Misc               = [misc.MATLAB   'tb_Misc'         filesep];
misc.Music_acoustics    = [misc.MATLAB   'Music_acoustics' filesep];
misc.Stats              = [misc.MATLAB   'Stats'           filesep];
misc.Text               = [misc.MATLAB  'Text'             filesep];
misc.tb_AFC             = [misc.MATLAB   'tb_AFC'          filesep];
misc.tb_AM              = [misc.MATLAB   'tb_AM'           filesep];
misc.tb_AM_git          = '/home/alejandro/Documents/MATLAB/amtoolbox-code/';
misc.tb_AM_AddOns       = [misc.MATLAB   'tb_AM_AddOns'    filesep];
misc.tb_APEX            = [misc.MATLAB   'tb_APEX'         filesep];
misc.tb_APEX_tools      = [misc.tb_APEX  'tools'           filesep];
misc.tb_NMT             = [misc.MATLAB   'tb_NMT_4.31'     filesep];
misc.tb_Plot4papers_JU  = [misc.MATLAB   'tb_Plot4papers_JU' filesep]; % Developed by Jaime Undurraga
misc.tb_UR              = [misc.MATLAB   'tb_UR_EAR_v2_1'    filesep];

misc.Verhulstetal_debug        = [misc.MATLAB 'Verhulstetal_debug'         filesep];

if nargin == 1
    % These folders will not be added if Start_UGent is run
    misc.Verhulstetal2017Model     = [misc.Verhulstetal_debug 'v110_2017'      filesep]; %'Verhulstetal2017Model' filesep];
    misc.Verhulstetal2018Model     = [misc.Verhulstetal_debug 'v110_2018'      filesep]; % It is theoretically the same as 'v110_2017'
    misc.Verhulstetal2018Model_corr= [misc.Verhulstetal_debug 'v115_2018_corr' filesep];
    % misc.Verhulstetal2018Model_ic_Carney = [misc.MATLAB 'Verhulstetal2018Model_ic_Carney' filesep]; % commented on 1/07/2019

    misc.Verhulstetal_v110         = [misc.Verhulstetal_debug 'v110_2018'      filesep];
    misc.Verhulstetal_v120         = [misc.Verhulstetal_debug 'v120_2018'      filesep];
    misc.Verhulstetal_v120_corr    = [misc.Verhulstetal_debug 'v120_2018_corr' filesep]; % created on 22/10/2019
    misc.Verhulstetal_v120_corr_py = [misc.Verhulstetal_debug 'v120_2018_corr_py' filesep]; % created on 14/04/2020
    misc.Verhulstetal2019Model     = [misc.Verhulstetal_debug 'v125_2019'      filesep];
    misc.Verhulstetal_v125         = misc.Verhulstetal2019Model;
end

misc.Verhulstetal_current      = [misc.MATLAB 'Verhulstetal_current' filesep];

if (nargin==1)
    if (~isfield(misc,type))
        error('Invalid type');
    end
    misc=misc.(type);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
