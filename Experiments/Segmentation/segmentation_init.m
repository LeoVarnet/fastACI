function cfg_inout = segmentation_init(cfg_inout)
% function cfg_inout = speechLAMIv2_init(cfg_inout)
%
% This should only be run once, to create the wave files.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
dir_speech_orig = [fastACI_basepath 'Stimuli' filesep 'LAMI' filesep];
% dir_speech_orig = '/home/alejandro/Documents/Documenten/04-Ontvangen/20220512-Leo-LAMIE/'; warning('Temporal solution...')
% if ~exist(dir_speech_orig,'dir')
%     dir_speech_orig = [pwd filesep];
% end
dBFS       = cfg_inout.dBFS;
% lvl_target = cfg_inout.SPL;
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
  
if bGenerate_stimuli
   % Speaker_ID = cfg_inout.Cond_extra_2;
    %files = Get_filenames(dir_speech_orig,['*.wav']);%files = Get_filenames(dir_speech_orig,[Speaker_ID '*.wav']);
    switch lower(cfg_inout.Condition)
        case {'lami', 'lami_shifted'}
            files = {'l_amie.wav', 'la_mie.wav'};
            meanf0 = mean([216.1, 205.6]);
        case {'lapel', 'lapel_shifted'}
            files = {'l_appel.wav', 'la_pelle.wav'};
            meanf0 = mean([215.2, 217.7]);
        case 'lapesanteur'
            files = {'l_apesanteur.wav', 'la_pesanteur.wav'};
            meanf0 = mean([212.4, 208]);
        case 'latension'
            files = {'l_attension.wav', 'la_tension.wav'};
            meanf0 = mean([220, 206.5]);
        case 'lacroch'
            files = {'l_accroche.wav', 'la_croche.wav'};
            meanf0 = mean([210.8, 204.3]);
        case 'lalarm'
            files = {'l_alarme.wav', 'la_larme.wav'};
            meanf0 = mean([215.6, 200.1]);
        case 'alapel1'
            files = {'v1_l_appel.wav', 'v1_la_pelle.wav'};
            meanf0 = mean([194.9, 183.3]);
        case 'alapel2'
            files = {'v2_l_appel.wav', 'v2_la_pelle.wav'};
            meanf0 = mean([191, 180]);
        case 'alapel3'
            files = {'v3_l_appel.wav', 'v3_la_pelle.wav'};
            meanf0 = mean([192.1, 177.5]);
        otherwise
            error('Stimuli condition not recognized. The third argument of fastACI_experiment ''segmentation'' should be one of the following: ''LAMI'', ''LAMI_shifted'', ''LAPEL'',  ''LAPEL_shifted'', ''LAPESANTEUR'', ''LATENSION'', ''LACROCH'', ''LALARM'', ''ALAPEL1'', ''ALAPEL2'', ''ALAPEL3''\n')
    end
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
        
        % idx = strfind(files{i},'_eq');
        % if ~isempty(idx)
        %     fprintf('\t%s: suffix ''eq'' found and being removed from the speech sample...\n',mfilename)
        %     files{i} = [files{i}(1:idx-1) '.wav'];
        % end
         
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

    cfg_inout.timevec = nan(Nseg+1,cfg_inout.N);      % initialisation
    cfg_inout.f0vec   = nan(size(cfg_inout.timevec)); % initialisation
end
for i = 1:cfg_inout.N
    if bGenerate_stimuli         
        seed_number = cfg_inout.seeds_order(i);
        rng(seed_number); % insig = randn(N_samples,1); rng(seed_number)
         
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Stimulus generation:
        files = Get_filenames(dir_target,'*.wav');
        if i < cfg_inout.N/2
            % loading L'amie
            [insig,fs] = audioread([dir_target files{1}]);
        else
            % loading La mie
            [insig,fs] = audioread([dir_target files{2}]); 
        end
        
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

        % f0vec is the random vector of f0 shifts (in cents). Burr et al. used a
        % S.D. of 100 cents, clipped at +/-2.2 S.D.
        f0vec = il_randn_clipped(Nseg+1,1,2.2); % Nseg+1 because we want one f0 shift at each segment edge
        f0vec = 100*f0vec;
        % timevec is the random vector of time shifts (in sec). I chose a S.D. of
        % 15 ms (remember that each segment is ~800 ms long and that very compressed
        % segments sound unnatural), clipped at +/-2.2 S.D.
        timevec = [0; il_randn_clipped(Nseg-1,1,2.2); 0]; % Nseg-1 because we want one time shift at each segment edge except the onset and offset
        timevec = 0.015*timevec;

        % Calling WORLD. The fourth parameter is the spectrum modification. 
        % We (at ENS) will not play with this dimension
        outsig = il_WorldSynthesiser(insig, fs, f0vec, 1, timevec, do_shift, meanf0);
        cfg_inout.timevec(:,i) = timevec;
        cfg_inout.f0vec(:,i)   = f0vec;
        
        %%% plot and play the result
        % figure;
        % sgram(outsig, fs)
        % sound(outsig, fs);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % lvls_before_noise(i) = rmsdb(insig)+dBFS;
        % % insig = scaletodbspl(insig,lvl_target,dBFS);
        % lvls_after_noise(i) = rmsdb(insig)+dBFS;

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

function y = il_WorldSynthesiser(x, fs, f0_param, spec_param, time_param, do_shift, meanf0)
% function y = il_WorldSynthesiser(x, fs, f0_param, spec_param, time_param)
%
% WorldSynthesizer_Leo(x, fs, 1, 1, 1);
%
% f0_param is the pitch shift value in cents
% time_param is the time shift value in sec
%
% Adapted from WorldSynthesizer_Leo.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<6
    do_shift = 0; % if do_shift = 1, then the first segment edge is not zero
end
if nargin<7
    meanf0 = 210; % Desired mean f0
end

f0_parameter = Harvest(x, fs);
spectrum_parameter = CheapTrick(x, fs, f0_parameter);
source_parameter = D4C(x, fs, f0_parameter);

% Generate the vector of modifications
fs_word = 1/(f0_parameter.temporal_positions(2)-f0_parameter.temporal_positions(1));
Nsamples = length(f0_parameter.temporal_positions);
Ninput = length(f0_param);
Nseg = Ninput-1;
Nsamp_in_seg = 0.1*fs_word;%floor(Nsamples/Nseg);

f0_vect = ones(1,Nsamples);
time_vect = source_parameter.temporal_positions;
if ~do_shift
    % normal
    for i_segment=1:Nseg
        f0_vect((i_segment-1)*Nsamp_in_seg+1:i_segment*Nsamp_in_seg+1) = linspace(f0_param(i_segment),f0_param(i_segment+1),Nsamp_in_seg+1);
        time_vect((i_segment-1)*Nsamp_in_seg+1:i_segment*Nsamp_in_seg+1) = linspace(source_parameter.temporal_positions((i_segment-1)*Nsamp_in_seg+1)+time_param(i_segment),source_parameter.temporal_positions(i_segment*Nsamp_in_seg)+time_param(i_segment+1),Nsamp_in_seg+1);
    end
else
    % shifted
    shift_in_samp = Nsamples-Nseg*Nsamp_in_seg;
    for i_segment=1:Nseg
        f0_vect(shift_in_samp+(i_segment-1)*Nsamp_in_seg:shift_in_samp+i_segment*Nsamp_in_seg) = linspace(f0_param(i_segment),f0_param(i_segment+1),Nsamp_in_seg+1);
        time_vect(shift_in_samp+(i_segment-1)*Nsamp_in_seg:shift_in_samp+i_segment*Nsamp_in_seg) = linspace(source_parameter.temporal_positions(shift_in_samp+(i_segment-1)*Nsamp_in_seg)+time_param(i_segment),source_parameter.temporal_positions(shift_in_samp+i_segment*Nsamp_in_seg)+time_param(i_segment+1),Nsamp_in_seg+1);
    end
end

% f0 kept
%source_parameter.f0 = source_parameter.f0 .* 2.^(f0_vect/1200);
% f0 trajectory kept but mean f0 neutralized
origf0 = source_parameter.f0; origf0(origf0==0) = nan; origf0 = mean(origf0,'omitnan');
source_parameter.f0 = source_parameter.f0 - origf0 + meanf0; % neutralize mean
source_parameter.f0 = source_parameter.f0 .* 2.^(f0_vect/1200);
% f0 flattened and mean f0 neutralized
%source_parameter.f0(source_parameter.f0==0) = nan;
%source_parameter.f0 = meanf0 .* 2.^(f0_vect/1200);

fft_size = (size(spectrum_parameter.spectrogram, 1) - 1) * 2;
w = (0 : fft_size - 1) * fs / fft_size;
w2 = (0 : fft_size / 2) * fs / fft_size / spec_param;
for i = 1 : size(spectrum_parameter.spectrogram, 2)
  tmp = [spectrum_parameter.spectrogram(:, i); spectrum_parameter.spectrogram(end - 1 : -1 : 2, i)];
  spectrum_parameter.spectrogram(:, i) = interp1(w, tmp, w2, 'linear', 'extrap');
end

source_parameter.temporal_positions = time_vect;%source_parameter.temporal_positions .* time_vect

y = Synthesis(source_parameter, spectrum_parameter);

%audiowrite(output_filename, y, fs);
function y = il_randn_clipped(N,M,threshold)
y = randn(N,M);
while any(any(abs(y)>threshold))
    y(abs(y)>threshold) = randn;
end