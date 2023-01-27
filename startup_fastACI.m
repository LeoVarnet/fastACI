function startup_fastACI
% function startup_fastACI
%
% Programmed by Alejandro Osses, ENS 2021-
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bAMT         = 1;
bFastACI_exp = 1;
bFastACI_sim = 1; % only if found... (Leo's or Alejandro's)

dir_fastACI  = [fileparts(which(mfilename)) filesep];

paths = il_get_subpaths(bFastACI_exp,dir_fastACI);

addpath(paths.Utility);
Add_paths(paths);

cur_dir = [pwd filesep];
if bAMT
    try
        dir_AMT = fastACI_dir_amtoolbox; % [userpath filesep 'amtoolbox-code' filesep];
    catch me
        fastACI_set_amtoolbox;
        dir_AMT = fastACI_dir_amtoolbox;
    end
    
    if exist(dir_AMT,'dir')
        addpath([dir_AMT filesep]);
        amt_start;

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

try
    % Checks if dir_data exists...
    fastACI_dir_data;
catch me
    % ... and creates dir_data if it does not exist
    fastACI_set_dirdata;
end

% This is an optional toolbox, to be able to run the Intellitest:
bAFC = exist('fastACI_dir_AFCtoolbox','file');
if bAFC
    bAFC = input('The optional AFC toolbox can be added, do you want that (1=yes, 0=no): ');
    
    if bAFC
        paths = il_get_subpaths_AFC;
        Add_paths(paths);
        
        dir_dst = [paths.afc 'base' filesep];
        file_src = [fastACI_basepath 'tb_AFC_AddOns' filesep 'procedures' filesep 'afc_proto_to_copy.m'];
        info_afc_proto_custom = dir(file_src);
        
        file_dst = [dir_dst 'afc_proto.m'];
        info_afc_proto        = dir(file_dst);
        
        if info_afc_proto.bytes == info_afc_proto_custom.bytes
            % Nothing to do
        else
            movefile(file_dst,[file_dst(1:end-2) '_orig.m'])
            copyfile(file_src,file_dst);
        end
        disp('')
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        file_src = [paths.afc 'scripts' filesep 'rms.m'];
        if exist(file_src,'file')
            % This will avoid the shadowing of the MATLAB rms function
            file_dst = [paths.afc 'scripts' filesep 'rms_orig.m'];
            movefile(file_src,file_dst);
        end
    end
end

fastACI_path = [fileparts(fileparts(fastACI_basepath)) filesep 'fastACI_sim' filesep];
if exist(fastACI_path,'dir')
    fprintf('fastACI_sim toolbox found: only the tb_fastACI_AddOns folder will be added...\n');
    addpath([fastACI_path 'MATLAB' filesep 'tb_fastACI_AddOns' filesep]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp([mfilename '.m: startup file for the fastACI toolbox'])

function paths = il_get_subpaths(bFastACI_exp,dir_fastACI)

if bFastACI_exp == 1
    % addpath([dir_fastACI     'Utility_paths' filesep]);
    
    paths.fast_ACI             = dir_fastACI;
    paths.Base                 = [dir_fastACI       'Base'                 filesep];
    paths.Calibration          = [dir_fastACI       'Calibration'          filesep];
    paths.Defaults             = [dir_fastACI       'Defaults'             filesep];
    paths.Experiments          = [dir_fastACI       'Experiments'          filesep];
    paths.Publications         = [dir_fastACI       'Publications'         filesep];
    paths.publ_osses2021c      = [paths.Publications 'publ_osses2021c'     filesep];
    paths.publ_osses2022b      = [paths.Publications 'publ_osses2022b'     filesep];
    paths.publ_osses2022c      = [paths.Publications 'publ_osses2022c'     filesep];
    paths.publ_varnet2013      = [paths.Publications 'publ_varnet2013'     filesep];
    paths.publ_varnet2022a     = [paths.Publications 'publ_varnet2022a'    filesep];
    paths.pres_osses2022_02    = [paths.Publications 'pres_osses2022_02'   filesep];
    paths.Interface            = [dir_fastACI       'Interface'            filesep];
    paths.legacy               = [dir_fastACI       'legacy'               filesep];
    paths.Local                = [dir_fastACI       'Local'                filesep];
    paths.modulationACI        = [paths.Experiments 'modulationACI'        filesep];
    paths.modulationFM         = [paths.Experiments 'modulationFM'         filesep];
  % paths.modulationACI_seeds  = [paths.Experiments 'modulationACI_seeds'  filesep]; % removed on 25/03/2022
    paths.Praat                = [dir_fastACI       'Praat'                filesep];
    paths.speechACI_Logatome   = [paths.Experiments 'speechACI_Logatome'   filesep];
    paths.speechACI_varnet2013 = [paths.Experiments 'speechACI_varnet2013' filesep];
    paths.speechACI_varnet2015 = [paths.Experiments 'speechACI_varnet2015' filesep];
    paths.segmentation         = [paths.Experiments 'segmentation'         filesep];
    paths.Scripts              = [dir_fastACI       'Scripts'              filesep];
    paths.Simulations          = [dir_fastACI       'Simulations'          filesep];
    paths.Model                = [paths.Simulations 'Model'                filesep];
    paths.Model_extra          = [paths.Model       'relanoiborra2019_preproc_debug' filesep];
    paths.Model_arg            = [paths.Simulations 'Model_arg'            filesep];
    paths.Model_stages         = [paths.Simulations 'Model_stages'         filesep];
    paths.Model_decisions      = [paths.Simulations 'Model_decisions'      filesep];
    paths.Stats                = [dir_fastACI       'Stats'                filesep];
    paths.Stats_tb_GLM_penalty = [paths.Stats       'tb_GLM_penalty'       filesep];
    paths.Stats_tb_optim_legacy= [paths.Stats       'tb_optim_legacy'      filesep];
    paths.Stim_generation      = [dir_fastACI       'Stim_generation'      filesep];
    paths.tb_ACI               = [dir_fastACI       'tb_ACI'               filesep];
    paths.tb_AFC_AddOns        = [dir_fastACI       'tb_AFC_AddOns'        filesep];
    paths.tb_AFC_AddOns_exp    = [paths.tb_AFC_AddOns 'experiments'        filesep];
    paths.tb_AFC_AddOns_proc   = [paths.tb_AFC_AddOns 'procedures'         filesep];
else
    error('%s: bFastACI_exp must be 1)',upper(mfilename));
end
paths.Plotting = [dir_fastACI 'Plotting'        filesep];
paths.Utility  = [dir_fastACI 'Utility' filesep]; % This folder is always added

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function paths = il_get_subpaths_AFC

path = fastACI_dir_AFCtoolbox;

paths.afc    = path;
paths.addons = [path 'addons' filesep];
paths.base   = [path 'base'   filesep];
paths.calibration = [path 'calibration' filesep];
paths.experiments = [path 'experiments' filesep];
paths.gui    = [path 'gui' filesep];
paths.models = [path 'models' filesep];
paths.scripts = [path 'scripts' filesep];
% addpath([fileparts(which(mfilename)), filesep, 'soundmexpro', filesep, 'bin'])
% addpath([fileparts(which(mfilename)), filesep, 'procedures'])
% addpath([fileparts(which(mfilename)), filesep, 'sessions'])  
