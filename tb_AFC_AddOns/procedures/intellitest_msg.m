% intellitest_msg - based on olsa_msg.m, which is the english message 
%   definition file Version 1.30.0, last modified 16.04.2013 09:36 from AFC toolbox
%------------------------------------------------------------------------------

msg=struct(...
    'measure_msg','Beginning measurement',	...
    'correct_msg','--- CORRECT ---',			...
    'false_msg','--- WRONG ---',				...
    'maxvar_msg','Maximum level reached',	...
    'minvar_msg','Minimum level reached' ...
    );

msg.ready_msg    = {'Select the words from', ...
    'the matrix and press OK to continue'};
msg.start_msg    = {'You have started a new measurement.', ...
    'Press any key to continue.'};
msg.next_msg     = {'End of Run.', ...
    'Press "s" for a new run or "e" to end.'};
msg.finished_msg = {'Experiment Done.', ...
    'Press "e" to end.'};

msg.experiment_windetail = 'Experiment: %s';
msg.measurement_windetail = 'Measurement %d of %d';
msg.measurementsleft_windetail = '%d of %d measurements left';

msg.okButtonString = 'OK';
msg.allCorrectButtonString = 'Alles verstanden';

msg.startButtonString = 'Start';
msg.endButtonString = 'End';

% List copied from ../Speech/dutch/Matrix/woordenmatrix.xls
msg.buttonString = {...
'ababa'   ,	'agaga', 'amama', 'assassa'; ...
'achacha' , 'ajaja', 'anana', 'atata'; ...
'adada'   , 'akaka', 'apapa', 'avava'; ...
'afafa'   , 'alala', 'arara', 'azaza'};

% file:///home/alejandro/Documents/Databases/Toolbox/20211011-Intellitest/VCVCVs/ababa1a.wav, ababa2a.wav, ababa3a.wav
% file:///home/alejandro/Documents/Databases/Toolbox/20211011-Intellitest/VCVCVs/achacha1a.wav, achacha2a.wav, achacha3a.wav
% file:///home/alejandro/Documents/Databases/Toolbox/20211011-Intellitest/VCVCVs/adada1a.wav, adada2a.wav, adada3a.wav
% file:///home/alejandro/Documents/Databases/Toolbox/20211011-Intellitest/VCVCVs/afafa1a.wav, afafa2a.wav, afafa3a.wav
% file:///home/alejandro/Documents/Databases/Toolbox/20211011-Intellitest/VCVCVs/agaga1a.wav, agaga2a.wav, agaga3a.wav
% file:///home/alejandro/Documents/Databases/Toolbox/20211011-Intellitest/VCVCVs/ajaja1a.wav, ajaja2a.wav, ajaja3a.wav
% file:///home/alejandro/Documents/Databases/Toolbox/20211011-Intellitest/VCVCVs/akaka1a.wav, akaka2a.wav, akaka3a.wav
% file:///home/alejandro/Documents/Databases/Toolbox/20211011-Intellitest/VCVCVs/alala1a.wav, alala2a.wav, alala3a.wav
% file:///home/alejandro/Documents/Databases/Toolbox/20211011-Intellitest/VCVCVs/amama1a.wav, amama2a.wav, amama3a.wav
% file:///home/alejandro/Documents/Databases/Toolbox/20211011-Intellitest/VCVCVs/anana1a.wav, anana2a.wav, anana3a.wav
% file:///home/alejandro/Documents/Databases/Toolbox/20211011-Intellitest/VCVCVs/apapa1a.wav, apapa2a.wav, apapa3a.wav
% file:///home/alejandro/Documents/Databases/Toolbox/20211011-Intellitest/VCVCVs/arara1a.wav, arara2a.wav, arara3a.wav
% file:///home/alejandro/Documents/Databases/Toolbox/20211011-Intellitest/VCVCVs/assassa1a.wav, assassa2a.wav, assassa3a.wav
% file:///home/alejandro/Documents/Databases/Toolbox/20211011-Intellitest/VCVCVs/atata1a.wav, atata2a.wav, atata3a.wav
% file:///home/alejandro/Documents/Databases/Toolbox/20211011-Intellitest/VCVCVs/avava1a.wav, avava2a.wav, avava3a.wav
% file:///home/alejandro/Documents/Databases/Toolbox/20211011-Intellitest/VCVCVs/azaza1a.wav, azaza2a.wav, azaza3a.wav

% eof
