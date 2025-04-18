function [Data_matrix, cfg_inout] = fastACI_getACI_dataload(cfg_inout, ListStim, cfg_game, data_passation, varargin)
% [Data_matrix, tf] = fastACI_getACI_dataload(cfg, ListStim) 
%
% Load stimuli data and generates appropriate matrix Data_matrix for the
% calculation of classification image
% Data_matrix is a N_TrialsLoad X DimCI matrix (Data_matrix(i,f,t) or
% Data_matrix(i,t) or Data_matrix(i,f), where i is the stimulus
% presentation number as specified in ordre_aleatoire ; tf is a structure
% array containing the time and/or frequency index
%
% varargin : 
% matrix_crea(cfg, ListStim, 'dimonly', 'yes') does not calculate
% Data_matrix
% matrix_crea(cfg, ListStim, 'zscore', 'no') does not perform a z-scoring
% on the resulting Data_matrix
%
% - zscore option removed from this script on 8/04/2021, as that would be
%       a preprocessing of the data and not a 'data load'
%
% - cfg_game is only used if the option WithSNR is requested...
% - data_passation is only used if the option WithSNR is requested...
%
% Authors: Leo Varnet and Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% TODO:
%     - Less or equal in assignment of frequencies
%     - There are frequencies were no bin is counted
%          step_df = fs/NFFT = 5.7422 Hz
%          min(diff(f_speclogedges)) = 1.29 Hz BUT step_df should be less...

if ~isfield(cfg_inout,'keyvals') || ~isfield(cfg_inout,'flags')
    definput.import={'fastACI_getACI'};
    [cfg_inout.flags,cfg_inout.keyvals]  = ltfatarghelper([],definput,varargin);
    % [cfg_inout.flags,cfg_inout.keyvals]  = ltfatarghelper({},definput,varargin);

    fprintf('%s: Default values are being loaded',upper(mfilename));
end

dimonly = 0;
numerase = 0;

path = [];
subpath = [];

TF_type = cfg_inout.flags.TF_type;

if isfield(cfg_inout,'FolderInBruit')
    cfg_inout.dir_noise = cfg_inout.FolderInBruit;
    warning('Old variable naming for ''dir_noise'' (old name: %s)','FolderInBruit');
    cfg_inout = rmfield(cfg_inout,'FolderInBruit');
end

while ~isempty(varargin)
    switch varargin{1}
        case 'dimonly'
            dimonly = isyes(varargin{2});
        otherwise
            error(['error : unrecognised argument : ' varargin{1}]);
    end
    varargin=varargin(3:end); % removes the already-read varargins
end

if isfield(cfg_inout,'N_trialselect')
    N_trialselect = cfg_inout.N_trialselect;
else
    N_trialselect = length(ListStim);
end

if ~isfield(cfg_inout,'idx_trialselect')
    cfg_inout.idx_trialselect = 1:N_trialselect;
end

switch TF_type
    case {'lyon','noise_logspect'}
        % This parameter is only needed for lyon and noise_logspect
        cfg_inout = Ensure_field(cfg_inout,'decimation',1);
end

WithSNR    = cfg_inout.keyvals.apply_SNR;
WithSignal = cfg_inout.keyvals.add_signal;

if ~strcmp(cfg_inout.dir_noise(end),filesep)
    cfg_inout.dir_noise = [cfg_inout.dir_noise filesep];
end

if WithSignal
    % Only checked if the targets need to be loaded
    if ~strcmp(cfg_inout.dir_target(end),filesep)
        cfg_inout.dir_target = [cfg_inout.dir_target filesep];
    end 
end

% -------------------------------------------------------------------------
% Extract time and/or frequency index
WavFile = [cfg_inout.dir_noise ListStim(1).name];
if exist(WavFile,'file')
    [bruit,fs] = audioread(WavFile);
else
    error('Sounds not found on disk, redefine dir_noise and/or dir_target')
    % speechACI_Logatome_init
end
bruit=mean(bruit,2);

%%%
if WithSignal
    
    fname_wav = Get_filenames(cfg_inout.dir_target,'*.wav');
    for i = 1:cfg_inout.N_target
        % Checking if the target names coincide
        if ~strfind(fname_wav{i},cfg_inout.target_names{i})
            error('Target sound that is being loaded (%s) does not match the corresponding target name (%s)',fname_wav{i},cfg_inout.target_names{i});
        end
        
        fname_full = [cfg_inout.dir_target fname_wav{i}];
        if exist(fname_full,'file')
            [insig_target(:,i),fs] = audioread(fname_full);
        else
            error('Sounds not found on disk, redefine dir_target')
        end
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gets 't' and 'f' for each T-F representation:
switch TF_type
    case {'spect','tf'}
        if strcmp(cfg_inout.flags.TF_type,'tf')
            TF_type = 'spect';
            warning('%s: flag ''tf'' will be soon deprecated, use ''spect'' instead...',upper(mfilename));
        end
        cfg_inout = Ensure_field(cfg_inout,'spect_overlap',0); 
        cfg_inout = Ensure_field(cfg_inout,'spect_Nwindow',512); 
        cfg_inout = Ensure_field(cfg_inout,'spect_NFFT',512);
        cfg_inout = Ensure_field(cfg_inout,'spect_unit','dB');
        
        Nwindow = cfg_inout.spect_Nwindow;
        overlap = cfg_inout.spect_overlap;
        NFFT    = cfg_inout.spect_NFFT;
        spect_unit = cfg_inout.spect_unit;
        [~,f,t,~] = spectrogram(bruit, Nwindow, Nwindow*overlap,NFFT,fs);
        
        t_correction = ((Nwindow/fs)/2)/2; % first time sample will be at 1/fs
        t = t-t_correction;
        
    case 'lyon'
        path = fileparts(which('data_load'));
        path = fileparts(path); % one level up
        path = [path filesep 'tb_AuditoryToolbox' filesep];
        addpath(path);

        subpath = [path 'src' filesep];
        addpath(subpath);

        Tcoch = 1/(fs/cfg_inout.decimation);
        t=(Tcoch:Tcoch:(length(bruit)/fs))-Tcoch/2;
        [~, f] = DesignLyonFilters(fs,cfg_inout.earQ,cfg_inout.stepfactor);
        f=f(end:-1:1);
        
    case 'noise_logspect'
        
        %%% Define as defaults:
        cfg_inout = Ensure_field(cfg_inout,'Nwindow',512); 
        Nwin = cfg_inout.Nwindow;
        Noverlap = (7/8)*Nwin;%0;%
        
        f_step = fs/512;
        f_spec = 0:f_step:fs-f_step; % Frequencies if we do a 512-point FFT
        N_log_bins = sum(f_spec>=cfg_inout.f_limits(1) & f_spec<=cfg_inout.f_limits(2));
        
        NFFT = 512*15; % Nfft parameter must be chosen to ensure that there are enough frequency bins for logarithmic resampling
        logspect_unit = cfg_inout.logspect_unit;
        
        % prepare resampling

        fcut = cfg_inout.f_limits;
        
        if fcut(1) == 0
            warning('Lower fcut frequency must me non-zero: a logarithmic spacing needs to be used');
            fcut(1) = 80;
        end
        
        logf_bin = linspace(log10(fcut(1)),log10(fcut(2)),N_log_bins); % 256
        logf_edges = mean([logf_bin(1:end-1);logf_bin(2:end)]);
        logf_edges = [2*logf_edges(1)-logf_edges(2), logf_edges, 2*logf_edges(end)-logf_edges(end-1)];

        f_speclogbin   = 10.^(logf_bin);
        f_speclogedges = 10.^(logf_edges);
        f = f_speclogbin; 
        f(1) = round(f(1)); f(end) = round(f(end)); % trick to always include first and last bin
        
        [~,f_spec,t] = spectrogram(bruit,Nwin,Noverlap,NFFT,fs); 
        for i_freq=1:length(f_speclogedges)-1
            freqidx = f_spec>=f_speclogedges(i_freq) & f_spec< f_speclogedges(i_freq+1);
            T(:,i_freq) = freqidx;freqidx/sum(freqidx);
        end

    case 'gammatone'
        basef = 8000;
        [fhigh, idx_min] = min([8000 fs/2]);
        if idx_min == 2
            warning('fhigh is limited to fs/2=%.0f Hz',fs/2);
        end
        flags_gamma = {'basef',basef,'flow',40,'fhigh',fhigh,'bwmul',cfg_inout.bwmul, ...
            'dboffset',100,'no_adt','binwidth',cfg_inout.binwidth, ...
            'no_outerear','no_middleear'};
        switch cfg_inout.glmfct
            case 'l1glm'
                if fhigh < 8000
                    flags_gamma{end+1} = 'fc_nextpow2';
                else
                    flags_gamma{end+1} = 'no_fc_nextpow2'; % this is the default
                end 
        end
        [outsig,f,t,extras] = Gammatone_proc(bruit,fs,flags_gamma{:});
        f = f(:); % column array
        
        cfg_inout.gamma_flags       = flags_gamma;
        cfg_inout.gamma_main_script = 'Gammatone_proc';
        
    case 'modulationACI_proc'
        basef = 1000;
        flags_gamma = {'basef',basef,'flow',950,'fhigh',1050,'bwmul',1, ...
            'dboffset',100,'no_adt','binwidth',cfg_inout.binwidth, ...
            'no_outerear','no_middleear', 'hilbert'};
        [outsig,f,t,extras] = Gammatone_proc(bruit,fs,flags_gamma{:});
        f = f(:); % column array
        
        cfg_inout.gamma_flags       = flags_gamma;
        cfg_inout.gamma_main_script = 'Gammatone_proc';
        
    case 'adapt'
        basef = 8000;
        flags_gamma = {'basef',basef,'flow',40,'fhigh',8000,'bwmul',cfg_inout.bwmul, ...
            'dboffset',100,'adt','binwidth',cfg_inout.binwidth, ...
            'no_outerear','no_middleear'};
        [outsig,f,t,extras] = Gammatone_proc(bruit,fs,flags_gamma{:});
        f = f(:); % column array
        
        cfg_inout.adapt_flags       = flags_gamma;
        cfg_inout.adapt_main_script = 'Gammatone_proc';
        
    case 'king2019'
        error('Under construction')
        
    case 'dau1997'
        error('Under construction')
        
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg_inout.f_limits_idx = find(f>=cfg_inout.f_limits(1) & f<=cfg_inout.f_limits(2));
f = f(cfg_inout.f_limits_idx);

cfg_inout.t_limits_idx = find(t>=cfg_inout.t_limits(1) & t<=cfg_inout.t_limits(2));
t = t(cfg_inout.t_limits_idx);
switch cfg_game.experiment
    case 'localisationILD'
        % Adding an exception:
        cfg_inout.t_limits_idx = [cfg_inout.t_limits_idx cfg_inout.t_limits_idx+length(t)];
        t = [t t+max(t)];
end

N_t = length(t);
N_f = length(f); 

if isfield(cfg_inout,'N')
    % Makes sure that N_TrialsLoad is less than N (is not the case in 
    %     case the experiment was not completed)
    N_trialselect = min(N_trialselect,cfg_inout.N);
end
    
Data_matrix = zeros(N_trialselect, N_f, N_t);

cfg_inout.f=f;
cfg_inout.N_f = N_f;
cfg_inout.t=t;
cfg_inout.N_t = N_t;
% End: Extract time and/or frequency index
% -------------------------------------------------------------------------

% Creating data matrix
if dimonly == 0
    fprintf('\n%% Creating data matrix %%\n');
    
    %%%
    if exist(cfg_inout.dir_noise,'dir')
        dir_noise = cfg_inout.dir_noise;
    else
        dir_noise = [cd filesep cfg_inout.dir_noise filesep];
    end
    if ~strcmp(dir_noise(end),filesep)
        dir_noise = [dir_noise filesep];
    end
    %%%

    for i=1:N_trialselect
        n_stim = data_passation.n_stim(cfg_inout.idx_trialselect(i));%cfg_inout.stim_order(cfg_inout.idx_trialselect(i));
        fprintf(repmat('\b',1,numerase));
        msg=sprintf('Loading stim no %.0f of %.0f\n', i,N_trialselect);
        fprintf(msg);
        numerase=numel(msg);
        
        file2load = [dir_noise ListStim(n_stim).name];
        WavFile = strcat(file2load);
        [bruit,fs] = audioread(WavFile);
        bruit=mean(bruit,2);
        
        if WithSignal
            switch cfg_inout.experiment
                case 'speechACI_Logatome'
                    if ~strcmp( cfg_inout.Condition,'bump' )
                        warning('Only the speechACI_Logatome with bump noises has been so far validated for this option')
                    end
                otherwise
                    warning('Only one experiment validated so far for this option (speechACI_Logatome)')
            end
            idx_here = cfg_game.n_targets_sorted(n_stim);
            trial = bruit + insig_target(:,idx_here);
        else
            % This is the default
            switch cfg_game.experiment
                case 'localisationILD'
                    data_passation_here = data_passation;
                    data_passation_here.i_current = i;
                    str_stim = localisationILD_user(cfg_game,data_passation_here);
                    
                    trial = str_stim.tuser(:); % collated in time
                    
                otherwise
                    trial = bruit;
            end
        end
        
        if WithSNR
            % error('WithSNR: Not validated yet')
            % % Pbruit = mean(bruit.^2);
            % % A=sqrt((Pbruit)*10^(ListStim(n_stim).RSB/10));
            % % bruit = bruit/A;
            if i == 1
                warning('Only validated for AM revcorr')
                data_pass_here = data_passation;
                cfg_game_here = modulationACI_set(cfg_game);
            end
            data_pass_here.i_current = i;
            str_stim = modulationACI_user(cfg_game_here,data_pass_here);
            
            trial = str_stim.tuser; % bruit + target_here;
        end
        
        switch TF_type
            % -------------------------------------------------------------
            case 'spect'
                [~,f_here,t_here,p_bruit] = spectrogram(trial,Nwindow,Nwindow*overlap,NFFT,fs);
                % cfg_inout.freq_analysis_index = find(f>=cfg_inout.freq_analysis(1) & f<=cfg_inout.freq_analysis(2));
                % cfg_inout.time_analysis_index = find(t>=cfg_inout.time_analysis(1) & t<=cfg_inout.time_analysis(2));
                p_bruit=p_bruit(cfg_inout.f_limits_idx,cfg_inout.t_limits_idx);

                switch spect_unit
                    case {'dB','db'}
                        p_bruit=10*log10(abs(p_bruit));
                    case 'linear'
                        % Nothing to do
                end
                
                Data_matrix(i, :, : ) = p_bruit;
            % -------------------------------------------------------------
            case 'lyon'
                if WithSignal
                    error('WithSignal: Not validated yet')
                    % WavFile = strcat(cd, '/', cfg_inout.FolderInSignal, '/', ListStim(n_stim).signal, '.wav');
                    % signal = audioread(WavFile);
                    % signal=mean(signal,2);
                    % [ stimulus,~] = Addition_RSB(signal, bruit, ListStim(n_stim).RSB);
                    % bruit=stimulus;
                end
                
                [bruit_coch,f_here] = LyonPassiveEar(trial,fs,cfg_inout.decimation,cfg_inout.earQ,cfg_inout.stepfactor,1,1);
                bruit_coch=bruit_coch(end:-1:1,:);
                f_here = f_here(end:-1:1);
                
                bruit_coch=bruit_coch(cfg_inout.f_limits_idx,cfg_inout.t_limits_idx);
                f_here = f_here(cfg_inout.f_limits_idx);
                
                Data_matrix(i, :, : ) = bruit_coch;
                
            % -------------------------------------------------------------    
            case 'noise_logspect'

                [~,f_here,t_here,p] = spectrogram(trial,Nwin,Noverlap,NFFT,fs); 
                %p=p(find(f_spec>=fcut(1) & f_spec<=fcut(end)),find(t_spec>=tcut(1) & t_spec<=tcut(end)),:);
                if i == 1
                    warning('Converting to dB')
                end
                
                switch logspect_unit
                    case {'dB','db'}
                        p=10*log10(abs(p));
                    case 'linear'
                        % Nothing to do
                end
     
                pdB_logf2 = T'*p(:,cfg_inout.t_limits_idx);
                
                %resampling
                Data_matrix(i,:,:)=pdB_logf2(:,1:cfg_inout.decimation:end);
 
                if i == 1
                    % Makes sure that the time variable is decimated
                    cfg_inout.t=cfg_inout.t(1:cfg_inout.decimation:end);
                end
                
            % -------------------------------------------------------------    
            case {'gammatone','adapt','modulationACI_proc'}
                switch cfg_game.experiment
                    case 'localisationILD'
                        outsig = Gammatone_proc(trial,fs,flags_gamma{:});
                        outsig = transpose(outsig); % permute(outsig,[2 1]); % put time in the second dimension and frequency in the first one
                        
                    otherwise
                        outsig = Gammatone_proc(trial,fs,flags_gamma{:});
                        outsig = transpose(outsig); % permute(outsig,[2 1]); % put time in the second dimension and frequency in the first one
                end
                
                % outsig=20*log10(abs(outsig));
                % idx = find(isnan(outsig));
                
                Data_matrix(i, :, : ) = outsig(cfg_inout.f_limits_idx, cfg_inout.t_limits_idx);
                        
        end % end switch
    end % end for
    
    fprintf('\n');
    
else
    error('dimonly: Not validated yet')
    Data_matrix =[];
end
   
if ~isempty(path)
    rmpath(path);
end
if ~isempty(subpath)
    rmpath(subpath);
end

end
