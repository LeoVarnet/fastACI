function cfg_inout = replication_ahumada1975_init(cfg_inout)
% function cfg_crea = replication_ahumada1975_init(cfg_crea)
%
% Function comparable to *_init.m functions from AFC toolbox
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 0
    help toneinnoise_ahumada1975_init;
    return;
end

% cfg_inout = Check_cfg_crea_dirs(cfg_inout); % updates the folders to the local computer if
%                                             % crea comes from another computer
dir_subject = Get_fastACI_subject_dir(cfg_inout.experiment_full,cfg_inout.Subject_ID);
if ~exist(dir_subject,'dir')
    mkdir(dir_subject);
end

dBFS       = cfg_inout.dBFS;
fs         = cfg_inout.fs; 
%%% 1. Speech sounds:
% dir_target = [dir_subject 'tone-alone' filesep];
% lvl_target = cfg_inout.lvl_target;
ramp_dur   = cfg_inout.ramp_dur_noise;
N_ramp     = round(ramp_dur*fs);
 
% if isfield(cfg_inout,'dir_target')
%     if ~exist(cfg_inout.dir_target,'dir')
%         cfg_inout.dir_target = dir_target;
%     end
% else
%     cfg_inout.dir_target = dir_target;
% end
 
% if exist(dir_target,'dir')
%     bGenerate_stimuli = 0;
% else
%     % The original speech sounds are zero padded:
%     bGenerate_stimuli = 1;
% 
%     % If you are in this part of the code, 'dir_target' does not exist
%     mkdir(dir_target);
% end

% %%% Target does not need to be created
% if bGenerate_stimuli
%        
% else
%     
% end

% dur_target = size(insig,1)/fs;
% cfg_inout.dur_target = dur_target;
% 
% dur_target = cfg_inout.dur_target;
noise_type = cfg_inout.Condition;
N_presentation = cfg_inout.N_presentation;
N_target   = cfg_inout.N_target;
N          = cfg_inout.N; % N_presentation*N_target; 
%%% 1. Background noises:
dir_noise  = [dir_subject 'NoiseStims-' noise_type filesep];
cfg_inout.dir_noise = dir_noise;
lvl_noise  = cfg_inout.lvl_noise;
 
disp('The new stimuli will be generated using the following parameters:')
fprintf('\tNoise type=%s\n',noise_type);
% fprintf('\tDuration=%.2f s\n',dur_target);
fprintf('\tSampling frequency=%.2f Hz\n',fs);
fprintf('\tA total of %.0f stimuli will be generated: %.0f ''noises'', %.0f ''targets''\n',N,N_presentation,N_target);
fprintf('\tTarget folder: %s\n',dir_noise);
disp('Press any button to continue, press ctrl+C to abort')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
pause()
 
bGenerate_stimuli = 0;
if ~exist(dir_noise,'dir')
    bGenerate_stimuli = 1;
    mkdir(dir_noise);
end
 
cfg_inout.dir_data_experiment = [fileparts(dir_subject(1:end-1)) filesep];

if bGenerate_stimuli
    [cfg_inout,s_current, bGenerate_stimuli] = Check_seeds_and_initialise(cfg_inout);
end

% Initialises the seed numbers on the fly:
seed_number_max = 4*cfg_inout.N; % arbitrary number, seeds will range from 0 to 4*N
method = 'perm'; % method = 'unif';
  
s_current = rng; % gets current seed
  
type_seed = 'shuffle';
try
    rng(type_seed); % seed based on the computer's clock
catch me
    rng('default');
    rng(type_seed); % seed based on the computer's clock
end

switch method
    case 'unif'
        error('Not tested recently')
    case 'perm'
        numbers = randperm(seed_number_max);
        seed_numbers = numbers(1:cfg_inout.N); % takes only the first 'N' numbers
end

cfg_inout.seeds_order = seed_numbers; % to be used sequentially
cfg_inout.seeds_order_method = method;
 
if cfg_inout.randorder == 1
    cfg_inout.stim_order = randperm(cfg_inout.N); 
else
    error('Not tested recently')
end

%%% Seed set back
rng(s_current);
%%% 

%% Noise samples generation
ListStim = [];
for i=1:N
    if bGenerate_stimuli
        if i == 1
            N_samples = round(cfg_inout.dur_noise*fs);
             
            % Generates the cosine ramp to be applied to the sounds:
            rp = ones([N_samples 1]); 
            rp(1:N_ramp)         = rampup(N_ramp);
            rp(end-N_ramp+1:end) = rampdown(N_ramp);
        end
          
        seed_number = cfg_inout.seeds_order(i);
        rng(seed_number); % insig = randn(N_samples,1); rng(seed_number)
          
        insig = Generate_noise(N_samples,noise_type,fs);
          
        insig = scaletodbspl(insig,lvl_noise,dBFS);
        insig = rp.*insig;
   
        fname_part1 = 'Noise'; 
         
        % Ensuring that the noise iteration number has four characters:
        stimnumber=num2str(i);
        if i<10000
            stimnumber=['0' stimnumber];
        end
        if i<1000
            stimnumber=['0' stimnumber];
        end
        if i<100
            stimnumber=['0' stimnumber];
        end
        if i<10
            stimnumber=['0' stimnumber];
        end
        % Ensuring 'End'
        fprintf(['Generating noise #' stimnumber '\n']);

        fname = [fname_part1 '_' stimnumber '.wav'];
    else
        if nargout >= 1
            fname = files{i};
        end
    end
    
    if nargout >= 1
        ListStim(i).name = fname;
    end
     
    if bGenerate_stimuli
        audiowrite([dir_noise fname],insig,fs);
    end
end
cfg_inout.ListStim = ListStim;
