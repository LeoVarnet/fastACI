function startup_fastACI
% function startup_fastACI
%
% Programmed by Alejandro Osses, ENS 2021-
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dir_fastACI = [fileparts(which(mfilename)) filesep];

bAMT = 1;
bFastACI_exp = 1;

paths = il_get_subpaths(bFastACI_exp,dir_fastACI);

addpath(paths.Utility);

Add_paths(paths);

cur_dir = [pwd filesep];
if bAMT
    dir_AMT = [userpath filesep 'amtoolbox-code' filesep]; % dir_AMT = '/home/alejandro/Desktop/AMT_now/amtoolbox-full-1.0/';
    
    if ~exist(dir_AMT,'dir')
        dir_AMT = uigetdir(fileparts(which(mfilename)),'Please indicate the folder where AMT 1.0 is located');
        dir_AMT = [dir_AMT filesep];
    end
    
    if exist(dir_AMT,'dir')
        addpath([dir_AMT filesep]);
        amt_start; % amt_start_ENS ;

        addpath([ltfatbasepath 'signals' filesep]); % noise.m (pink noise)

        rmpath([dir_AMT 'legacy' filesep]);
    else
        if ~exist('amt_start','file')
            fprintf('%s: AMT Toolbox 1.0 was not added to the path, please modify this script (%s.m)\n',upper(mfilename),mfilename);
            fprintf('    to include the correct path of the AMT toolbox in your local computer\n');
        end
    end
end
cd(cur_dir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp([mfilename '.m: startup file for the fastACI toolbox'])

function paths = il_get_subpaths(bFastACI_exp,dir_fastACI)

if bFastACI_exp  == 0
    error('%s: bFastACI_exp must be 1)',upper(mfilename));
end

if bFastACI_exp
    % addpath([dir_fastACI     'Utility_paths' filesep]);
    
    paths.fast_ACI        = dir_fastACI;
    paths.Base            = [dir_fastACI 'Base'            filesep];
    paths.Defaults        = [dir_fastACI 'Defaults'        filesep];
    paths.Experiments     = [dir_fastACI 'Experiments'     filesep];
    paths.Publications    = [dir_fastACI 'Publications'    filesep];
    paths.Interface       = [dir_fastACI 'Interface'       filesep];
    paths.legacy          = [dir_fastACI 'legacy'          filesep];
    paths.modulationACI        = [paths.Experiments 'modulationACI'        filesep];
    paths.modulationACI_seeds  = [paths.Experiments 'modulationACI_seeds'  filesep];
    paths.Praat           = [dir_fastACI 'Praat'           filesep];
    paths.speechACI_Logatome   = [paths.Experiments 'speechACI_Logatome'   filesep];
    paths.speechACI_varnet2013 = [paths.Experiments 'speechACI_varnet2013' filesep];
    paths.speechACI_varnet2015 = [paths.Experiments 'speechACI_varnet2015' filesep];
    paths.Scripts         = [dir_fastACI 'Scripts'         filesep];
    paths.Simulations     = [dir_fastACI 'Simulations'     filesep];
    paths.Stats           = [dir_fastACI 'Stats'           filesep];
    paths.Stim_generation = [dir_fastACI 'Stim_generation' filesep];
    paths.tb_ACI          = [dir_fastACI 'tb_ACI'          filesep];
end
paths.Plotting = [dir_fastACI 'Plotting'        filesep];
paths.Utility  = [dir_fastACI 'Utility' filesep]; % This folder is always added