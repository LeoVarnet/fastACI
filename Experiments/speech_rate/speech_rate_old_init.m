function cfg_inout = segmentation_init(cfg_inout)
% function cfg_inout = segmentation_init(cfg_inout)
%
% 1. Description:
%      This should only be run once, to create the wave files.
%      The experiment segmentationDebug offers a more extended choice of
%        conditions to run. The experiment segmentationDebug is hosted in
%        the fastACI_sim (private) repository.
% Author: Leo Varnet
% Author: Alejandro Osses (documentation and code stylisation)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
dir_speech_orig = [fastACI_basepath 'Stimuli' filesep 'speechrate' filesep];

dBFS       = cfg_inout.dBFS;
dur_ramp   = cfg_inout.dur_ramp; 
 
cfg_inout = Ensure_field(cfg_inout,'Condition','');
cfg_inout = Check_cfg_crea_dirs(cfg_inout); % updates the folders to the local computer if
                                            % crea comes from another computer
if ~exist(cfg_inout.dir_data_experiment,'dir')
    error('Please define an existing cfg_inout.dir_data_experiment (Not found: %s)',cfg_inout.dir_data_experiment);
end
dir_subject = [cfg_inout.dir_data_experiment cfg_inout.Subject_ID filesep];
if ~exist(dir_subject,'dir')
    mkdir(dir_subject);
end

%% 1. Speech sounds:
dir_target = cfg_inout.dir_target;
 
if exist(dir_target,'dir')
    bGenerate_stimuli = Check_if_dir_is_empty(dir_target,'*.wav'); 
else
    % The original speech sounds are zero padded:
    bGenerate_stimuli = 1;

    % If you are in this part of the code, 'dir_speech' does not exist
    mkdir(dir_target);
end
  

files = {'6BPA_tiger_rages_80_100.wav'};

if bGenerate_stimuli
    
    for i = 1:length(files)
        [insig,fs]  = audioread([dir_speech_orig files{i}]);
                 
        if i == 1
            if fs ~= cfg_inout.fs
                error('Sounds do not have the same sampling frequency as specified in the *._set.m file');
            end
        end
        lvls_before(i) = rmsdb(insig)+dBFS;
        % insig = scaletodbspl(insig,lvl_target,dBFS); % level adjustment before the silence is added
        lvls_after(i) = rmsdb(insig)+dBFS;
         
        sil = zeros(round(dur_ramp*fs),1);
        insig  = [sil; insig; sil];
        
        audiowrite([dir_target files{i}],insig,fs);
    end
else
    files = Get_filenames(dir_target,'*.wav');
    [insig,fs] = audioread([dir_target files{1}]); % reading one speech file
    if fs ~= cfg_inout.fs
        error('Sounds do not have the same sampling frequency as specified in the %s_set.m file');
    end
end
N_samples = length(insig);
  
%%% 2. Noise sounds:
dir_noise = cfg_inout.dir_noise;
if exist(dir_noise,'dir')
    % if dir is empty then it will generate the noises:
    bGenerate_stimuli = Check_if_dir_is_empty(dir_noise,'*.wav'); 
else
    bGenerate_stimuli = 1;
end
 
if bGenerate_stimuli == 1
    if ~strcmp(dir_noise(end),filesep) % if last character is \ or / (it should be the case always)
        dir_noise = [fileparts(dir_noise(1:end-1)) filesep];
    end
end
  
if bGenerate_stimuli
    [cfg_inout,s_current, bGenerate_stimuli] = Check_seeds_and_initialise(cfg_inout);
end

if bGenerate_stimuli == 0
    dir_noise = cfg_inout.dir_noise;
    
    files = Get_filenames(dir_noise,'*.wav');
else
    % If you are in this part of the code, 'dir_noise' does not exist:
    mkdir(dir_noise);
     
    s_current = rng; % gets current seed
end
 
ListStim = [];
if bGenerate_stimuli
    Nsamp_in_seg = 0.1*fs;
    Nseg = floor(N_samples/Nsamp_in_seg);

    cfg_inout.scalevec = nan(Nseg,cfg_inout.N);      % initialisation
end
for i = 1:cfg_inout.N
    if bGenerate_stimuli         
        seed_number = cfg_inout.seeds_order(i);
        rng(seed_number); % insig = randn(N_samples,1); rng(seed_number)
         
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Stimulus generation:
        if i == 1
            files = Get_filenames(dir_target,'*.wav');
        end
   %     if i < cfg_inout.N/2
            % loading L'amie
            idx_target = 1;
    %    else
            % loading La mie
    %        idx_target = 2;
    %    end
        [insig,fs] = audioread([dir_target files{idx_target}]);
        
        if fs ~= cfg_inout.fs
            error('Sounds do not have the same sampling frequency as specified in the %s_set.m file');
        end

        %%% Random modification, parameters by Leo

        switch lower(cfg_inout.Condition)
            case 'lami_shifted'
                do_shift = 1;
            otherwise
                do_shift = 0;
        end

        % % f0vec is the random vector of f0 shifts (in cents). Burr et al. used a
        % % S.D. of 100 cents, clipped at +/-2.2 S.D.
        % f0vec = il_randn_clipped(Nseg+1,1,2.2); % Nseg+1 because we want one f0 shift at each segment edge
        % f0vec = 100*f0vec;
        % % timevec is the random vector of time shifts (in sec). I chose a S.D. of
        % % 15 ms (remember that each segment is ~800 ms long and that very compressed
        % % segments sound unnatural), clipped at +/-2.2 S.D.
        % timevec = [0; il_randn_clipped(Nseg-1,1,2.2); 0]; % Nseg-1 because we want one time shift at each segment edge except the onset and offset
        % timevec = 0.015*timevec;

        scalevec = [il_randn_clipped(Nseg,1,1)]; % Nseg-1 because we want one time shift at each segment edge except the onset and offset
        scalevec = 1.5+scalevec;

        % Calling WORLD. 
        switch i
            case {1, cfg_inout.N/2}
                % The following lines are to store the 'original signals' but
                %   still subjected to the il_WorldSynthesiser. This is useful
                %   for simulations, because then these signals will have 
                %   the same length as the processed sounds. In this sense
                %   you need to set factor to 0.
                factor = 1; % 0.2
                if factor == 0 || factor == 1
                    suff_here = '';
                elseif factor < 1
                    suff_here = ['-0p' num2str(100*factor)];
                end
                %scalevec = 0.8*ones([1, length(timevec)-1]);
                outsig = il_WorldSynthesiser(insig, fs, scalevec);
                dir_target_new = [dir_target(1:end-1) '-idle' filesep]; mkdir(dir_target_new);
                fname_here = [dir_target_new files{idx_target}(1:end-4) '-idle' suff_here '.wav'];
                audiowrite(fname_here,outsig,fs);
        end
        outsig = il_WorldSynthesiser(insig, fs, scalevec);
        cfg_inout.scalevec(:,i) = scalevec;
        %cfg_inout.f0vec(:,i)   = f0vec;
        
        fname_part1 = 'Speech-proc'; 
        number = i;
 
        % Ensuring that the noise iteration number has four characters:
        stimnumber=num2str(number);
        if number<10000
            stimnumber=['0' stimnumber];
        end
        if number<1000
            stimnumber=['0' stimnumber];
        end
        if number<100
            stimnumber=['0' stimnumber];
        end
        if number<10
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
        audiowrite([dir_noise fname],outsig,fs);
    end
end
if bGenerate_stimuli
    try
        fsound = [fastACI_basepath 'Stimuli' filesep 'ready.wav'];
        [fsound,fs_sound] = audioread(fsound);
        sound(fsound,fs_sound);
    end
end
if isfield(cfg_inout,'bRove_level')
    
    if ~isfield(cfg_inout,'Rove_level')
        max_rove = cfg_inout.Rove_range; % dB
        cfg_inout.Rove_level = (rand(1,cfg_inout.N)-0.5)*2*max_rove;
    end
end
 
if nargout >= 1
    cfg_inout.ListStim = ListStim;
end

%%% Seed set back
if bGenerate_stimuli
    rng(s_current);
end
%%% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = il_WorldSynthesiser(insig, fs, scale_param, do_shift)
% function y = il_WorldSynthesiser(insig, fs, scale_param)
%
% WorldSynthesizer_Leo(x, fs, 1, 1, 1);
%
% scale_param is the scale factor
%
% Adapted from WorldSynthesizer_Leo.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<6
    do_shift = 0; % if do_shift = 1, then the first segment edge is not zero
end

if ~exist('Harvest.m','file')
    script2look = 'fastACI_dir_world';
    fprintf('%s.m: WORLD toolbox not found, trying to find the script %s.m...\n',mfilename,script2look);
    
    try
        WORLD_path = []; % trick to debug
        eval(sprintf('WORLD_path=%s;',script2look));
        addpath(WORLD_path);
    catch
        error('Please install the WORLD toolbox')
    end
end
f0_parameter = Harvest(insig, fs);
spect_params = CheapTrick(insig, fs, f0_parameter);
src_params   = D4C(insig, fs, f0_parameter);

% Generate the vector of modifications
fs_world = 1/(f0_parameter.temporal_positions(2)-f0_parameter.temporal_positions(1));
Nsamples = length(f0_parameter.temporal_positions);
Ninput = length(scale_param);
Nseg = Ninput;
Nsamp_in_seg = 0.1*fs_world;%floor(Nsamples/Nseg);

time_vect = src_params.temporal_positions;

% segment edge version
    % for i_segment=1:Nseg
    %     idx_i = (i_segment-1)*Nsamp_in_seg+1;
    %     idx_f = i_segment*Nsamp_in_seg+1; 
    %     num_steps_in = Nsamp_in_seg+1;
    %     % 
    %     % in1 = f0_param(i_segment);
    %     % in2 = f0_param(i_segment+1);
    %     % f0_vect_segment = linspace(in1,in2,num_steps_in);
    %     % f0_vect(idx_i:idx_f) = f0_vect_segment; % the last element (at idx_f)
    %     %                 % will be overwritten by the same value of the next 
    %     %                 % segment (idx_i) except for the last segment.
    %     in1 = src_params.temporal_positions((i_segment-1)*Nsamp_in_seg+1)+time_param(i_segment);
    %     in2 = src_params.temporal_positions(i_segment*Nsamp_in_seg)+time_param(i_segment+1);
    %     time_vect(idx_i:idx_f) = linspace(in1,in2,num_steps_in);
    % end


    % scale factor version
    in2=0;
    for i_segment=1:Nseg
        idx_i = (i_segment-1)*Nsamp_in_seg+1;
        idx_f = i_segment*Nsamp_in_seg+1;
        num_steps_in = Nsamp_in_seg+1;

        % original duration of the segment
        segment_dur = src_params.temporal_positions(i_segment*Nsamp_in_seg) - src_params.temporal_positions((i_segment-1)*Nsamp_in_seg+1);
        in1 = in2; % end of previous segment %src_params.temporal_positions((i_segment-1)*Nsamp_in_seg+1)+time_param(i_segment);
        in2 = in1 + segment_dur*scale_param(i_segment);%src_params.temporal_positions(i_segment*Nsamp_in_seg)+time_param(i_segment+1);
        time_vect(idx_i:idx_f) = linspace(in1,in2,num_steps_in);
    end

% f0 kept
%source_parameter.f0 = source_parameter.f0 .* 2.^(f0_vect/1200);
% f0 trajectory kept but mean f0 neutralised
% origf0 = source_parameter.f0; origf0(origf0==0) = nan; origf0 = mean(origf0,'omitnan');
% source_parameter.f0 = source_parameter.f0 - origf0 + meanf0; % neutralise mean
% source_parameter.f0 = source_parameter.f0 .* 2.^(f0_vect/1200);
% f0 flattened and mean f0 neutralised
% src_params.f0 = meanf0 .* 2.^(f0_vect/1200);

% fft_size = (size(spect_params.spectrogram, 1) - 1) * 2;
% w = (0 : fft_size - 1) * fs / fft_size;
% w2 = (0 : fft_size / 2) * fs / fft_size / spec_param;
% for i = 1 : size(spect_params.spectrogram, 2)
%   tmp = [spect_params.spectrogram(:, i); spect_params.spectrogram(end - 1 : -1 : 2, i)];
%   spect_params.spectrogram(:, i) = interp1(w, tmp, w2, 'linear', 'extrap');
% end

src_params.temporal_positions = time_vect;

y = Synthesis(src_params, spect_params);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = il_randn_clipped(N,M,threshold)

y = randn(N,M);
while any(any(abs(y)>threshold))
    y(abs(y)>threshold) = randn;
end
