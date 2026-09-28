function [Data_matrix,cfg_inout] = toneinnoise_ahumada1975_dataload(cfg_inout, ListStim, cfg_game, data_passation)

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

TF_type = 'spect';

if isfield(cfg_inout,'FolderInBruit')
    cfg_inout.dir_noise = cfg_inout.FolderInBruit;
    warning('Old variable naming for ''dir_noise'' (old name: %s)','FolderInBruit');
    cfg_inout = rmfield(cfg_inout,'FolderInBruit');
end

if isfield(cfg_inout,'N_trialselect')
    N_trialselect = cfg_inout.N_trialselect;
else
    N_trialselect = length(ListStim);
end

if ~isfield(cfg_inout,'idx_trialselect')
    cfg_inout.idx_trialselect = 1:N_trialselect;
end

% WithSNR    = cfg_inout.keyvals.apply_SNR;
% WithSignal = cfg_inout.keyvals.add_signal;

if ~strcmp(cfg_inout.dir_noise(end),filesep)
    cfg_inout.dir_noise = [cfg_inout.dir_noise filesep];
end

% if WithSignal
%     % Only checked if the targets need to be loaded
%     if ~strcmp(cfg_inout.dir_target(end),filesep)
%         cfg_inout.dir_target = [cfg_inout.dir_target filesep];
%     end 
% end

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
% if WithSignal
%     
%     fname_wav = Get_filenames(cfg_inout.dir_target,'*.wav');
%     for i = 1:cfg_inout.N_target
%         % Checking if the target names coincide
%         if ~strfind(fname_wav{i},cfg_inout.target_names{i})
%             error('Target sound that is being loaded (%s) does not match the corresponding target name (%s)',fname_wav{i},cfg_inout.target_names{i});
%         end
%         
%         fname_full = [cfg_inout.dir_target fname_wav{i}];
%         if exist(fname_full,'file')
%             [insig_target(:,i),fs] = audioread(fname_full);
%         else
%             error('Sounds not found on disk, redefine dir_target')
%         end
%     end
%     
% end

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
        cfg_inout = Ensure_field(cfg_inout,'spect_unit','energy');
        
        Nwindow = cfg_inout.spect_Nwindow;
        overlap = cfg_inout.spect_overlap;
        NFFT    = cfg_inout.spect_NFFT;
        spect_unit = cfg_inout.spect_unit;
        [~,f,t,~] = spectrogram(bruit, Nwindow, Nwindow*overlap, NFFT, fs);
        
        t_correction = ((Nwindow/fs)/2)/2; % first time sample will be at 1/fs
        t = t-t_correction;
        
        T = [0:0.1:0.4];
        for idx_T = 1:length(T)
            t_idx(idx_T,:) = t>T(idx_T) & t<T(idx_T)+0.1;
        end
        f_idx = [ 14:20 ];
        F = f(f_idx);
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
N_T = length(T);
N_F = length(F); 

if isfield(cfg_inout,'N')
    % Makes sure that N_TrialsLoad is less than N (is not the case in 
    %     case the experiment was not completed)
    N_trialselect = min(N_trialselect,cfg_inout.N);
end
    
Data_matrix = zeros(N_trialselect, N_F, N_T);

cfg_inout.f=F;
cfg_inout.N_f = N_F;
cfg_inout.t=T;
cfg_inout.N_t = N_T;
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
        n_stim = cfg_inout.stim_order(cfg_inout.idx_trialselect(i));
        fprintf(repmat('\b',1,numerase));
        msg=sprintf('Loading stim no %.0f of %.0f\n', i,N_trialselect);
        fprintf(msg);
        numerase=numel(msg);
        
        file2load = [dir_noise ListStim(n_stim).name];
        WavFile = strcat(file2load);
        [bruit,fs] = audioread(WavFile);
        bruit=mean(bruit,2);
        
        trial = bruit;
        
        
                [~,f_here,t_here,p_bruit] = spectrogram(trial,Nwindow,Nwindow*overlap,NFFT,fs);
                % cfg_inout.freq_analysis_index = find(f>=cfg_inout.freq_analysis(1) & f<=cfg_inout.freq_analysis(2));
                % cfg_inout.time_analysis_index = find(t>=cfg_inout.time_analysis(1) & t<=cfg_inout.time_analysis(2));
                
                p_bruit=p_bruit(cfg_inout.f_limits_idx,cfg_inout.t_limits_idx);
                
                p_bruit = p_bruit(f_idx,:);
                p_bruit = p_bruit*t_idx';
                
                
                switch spect_unit
                    case {'dB','db'}
                        p_bruit=10*log10(abs(p_bruit));
                    case 'linear'
                        % Nothing to do
                    case 'energy'
                        p_bruit = p_bruit.^2;
                end
                
                Data_matrix(i, :, : ) = p_bruit;
            % -------------------------------------------------------------
       
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

