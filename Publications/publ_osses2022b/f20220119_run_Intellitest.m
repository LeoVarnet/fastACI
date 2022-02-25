function f20220119_run_Intellitest(subjectname, Condition)
% function f20220119_run_Intellitest(subjectname, Condition)
%
% Author: Alejandro Osses
%
% Original name: f20220119_run_Intellitest
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc

global global_vars

if nargin == 0
    subjectname = input('Enter the ID of the subject to be tested (e.g., ''S09''): '); % 'SAO'; %'SAO'; % modelname = 'king2019';
end
if nargin < 2
    Conditions = {'varspeech16','varspeech16_mod8'};
    Show_cell(Conditions);
    bInput = input('Choose the test condition from above: ');
    Condition = Conditions{bInput};
end

available_hardware = {'Petite-Cabine', 'Grande-Cabine'};
if ~ismac
    available_hardware{end+1} = 'Sony-MDR-Alejandro';
    available_hardware{end+1} = 'Sony-WH-Alejandro';
end
Show_cell(available_hardware);                    
bInput = input('Choose your hw config: ');
hardware_cfg = available_hardware{bInput};

lvl_target   = 79; % -18 dBFS peak level => -21 dBFS rms 
switch hardware_cfg
    case 'Sony-MDR-Alejandro'
        % My own old headphones
        lvl_from_SLM = 90.6; % for the left headphone
        
    case 'Sony-WH-Alejandro' 
        % NR headphones
        lvl_from_SLM = 91.7; % for the left headphone
        
    case 'Petite-Cabine'
        lvl_from_SLM = 82.3; % for the left headphone
        
    case 'Grande-Cabine'
        lvl_from_SLM = 82.1; % for  the left headphone, preamp zc 0032, Id No 23156 
end
dBFS = 100+(lvl_target-lvl_from_SLM);

global_vars.dBFS = dBFS;

global_vars.result_path = [fastACI_basepath 'Publications' filesep 'publ_osses2022b' filesep ...
    'data_' subjectname filesep '1-experimental_results' filesep 'intellitest' filesep];

afc_main('Intellitest',subjectname,Condition)

