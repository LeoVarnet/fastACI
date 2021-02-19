function cfg_inout = modulationACI_init(cfg_inout)
% function cfg_crea = modulationACI_init(cfg_crea)
%
% Function comparable to *_cfg.m functions from AFC toolbox
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dir_where  = [cfg_inout.path cfg_inout.folder_name filesep]; 
        
fs         = cfg_inout.fs;
dur        = cfg_inout.stim_dur;
noise_type = cfg_inout.noise_type;
N_noise    = cfg_inout.N_noise;
N_signal   = cfg_inout.N_signal;
N_total    = N_noise*N_signal; 

disp('The new stimuli will be generated using the following parameters:')
fprintf('\tNoise type=%s\n',noise_type);
fprintf('\tDuration=%.2f s\n',dur);
fprintf('\tSampling frequency=%.2f Hz\n',fs);
fprintf('\tA total of %.0f stimuli will be generated: %.0f ''noises'', %.0f ''conditions''\n',N_total,N_noise,N_signal);
fprintf('\tTarget folder: %s\n',dir_where);
disp('Press any button to continue, press ctrl+C to abort')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
pause()

if ~isdir(dir_where)
    mkdir(dir_where);
    bGenerate_stimuli = 1;
elseif ~isempty(dir(dir_where))
    warning('%s: Folder %s is not empty.',upper(mfilename),dir_where);
    bInit_participant_only = input('Do you to initialise a participant without generating new sounds? (1=yes; 0=no): ');
    bGenerate_stimuli = ~bInit_participant_only;

    if bGenerate_stimuli
        error('%s: Folder %s is not empty. Choose a new ''folder_name'' or remove manually the existing content in the target folder',upper(mfilename),dir_where);
    end
end

%% Noise samples generation
for i=1:N_total
    clc
    if bGenerate_stimuli
        fprintf('Creating noise stimulus # %.0f of %.0f\n',i,N_total);

        str_stim = [];
        str_inout = []; 
        str_inout.istarget = 0;
        str_inout.expvar   = 0; % idle
        str2eval = sprintf('str_stim=%s_user(str_inout,cfg_crea);',cfg_inout.experiment);
        eval(str2eval);
        noise_cal = str_stim.stim_noise_alone;
    end
    % Ensuring that the noise iteration number has four characters:
    stimnumber=num2str(i);
    if i<1000; stimnumber=['0' stimnumber]; end
    if i<100;  stimnumber=['0' stimnumber]; end
    if i<10;   stimnumber=['0' stimnumber]; end
    % Ensuring 'End'

    fname = ['Noise_' stimnumber '.wav'];
    ListStim(i).name = fname;
    if bGenerate_stimuli
        audiowrite([dir_where fname], noise_cal, fs);
    end
end
cfg_inout.ListStim = ListStim;