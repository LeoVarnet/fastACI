function publ_varnet2022a_utils(Subject_ID,type_action)
% function publ_varnet2022a_utils(Subject_ID,type_action)
%
% - I don't understand why '[tone, fs] = audioread([dir_target 'nontarget.wav']);' in il_noisetone_
% - New: level adjustment for targets before loading noises (for i = 1:length(fname_targets))
% - I don't understand why, when deriving the templates, normalised signamls (by RMS)
%   are used + 30-Hz LPF is not applied.
% - Why 0.5 in computeCI?: CIrand(:,:,i_rand) = 0.5 * (CIrand_H + CIrand_FA - CIrand_M - CIrand_CR); %%% From Leo's
% - Why 0.5 in computeCI?: CIboot
%
% Author: Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

experiment = 'modulationACI';
dir_subj = [fastACI_dir_data experiment filesep Subject_ID filesep];

file = Get_filenames(dir_subj,'save*.mat');
if length(file) == 1
    savefilename = [dir_subj file{1}];
else
    error('Problem finding savegame file...')
end

ListStim = []; % : [3000×1 struct]
cfg_game = []; % : [1×1 struct]
data_passation = []; % : [1×1 struct]
load(savefilename);
n_response = data_passation.n_response;
n_signal   = data_passation.N_signal(1:length(n_response));
m          = data_passation.m(1:length(n_response));
                
switch type_action
    case 'Get_Behavior'
        
        M_bins = -20.5:1:-0.5;
        n_window = 100;
        [N_m,m_edge] = histcounts(m, M_bins);

        for i_m = 1:length(m_edge)-1
            m_bin(i_m)=0.5*(m_edge(i_m)+m_edge(i_m+1));
            idx_m = find(m>m_edge(i_m) & m<m_edge(i_m+1));
            H(i_m) = sum(n_signal(idx_m)==2 & n_response(idx_m)==2);
            M(i_m) = sum(n_signal(idx_m)==2 & n_response(idx_m)==1);
            CR(i_m) = sum(n_signal(idx_m)==1 & n_response(idx_m)==1);
            FA(i_m) = sum(n_signal(idx_m)==1 & n_response(idx_m)==2);
        end

        for i = 0:(length(n_response)/n_window)-1
            m_windowed(i+1) = mean(m(i*n_window+1:(i+1)*n_window));
            response_windowed = n_response(i*n_window+1:(i+1)*n_window);
            signal_windowed = n_signal(i*n_window+1:(i+1)*n_window);
            bias_windowed(i+1) = mean(response_windowed);
            PC_targetpresent(i+1) = mean(response_windowed(signal_windowed==2))-1;
            PC_targetabsent(i+1) = 2-mean(response_windowed(signal_windowed==1));
            % RT_windowed(i+1) = mean(RT(i*n_window+1:(i+1)*n_window));
        end

        trialnum = 1:n_window:length(m);
        save([dir_subj 'Behavior'],'m_windowed','PC_targetpresent','PC_targetabsent', ...
            'bias_windowed','trialnum','N_m','m_edge', 'H', 'M', 'CR', 'FA');
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    case 'Get_CIt'
        bLeo = 1;
        bAlejandro = ~bLeo;
        
        % Analyse noise in bands
        basef = 1000;
        BW_ERB = 1;
        if bLeo
            fcut = ERB2f(f2ERB(basef)+BW_ERB/2*[-1 1]);
        end
        if bAlejandro
            fcut = audtofreq(freqtoaud(basef)+BW_ERB/2*[-1 1]); % here there is a slight difference (ERB2f)
        end

        Nchannel = length(fcut)-1; 
        undersampling = 100;
        fcut_noiseE = 30; % Hz
        dBFS =  93.6139; % based on cal signal which is 72.6 dB using Sennheiser HD650
        SNR  = -10; % dB
        lvl  =  65; % dB
        
        dir_target = [dir_subj 'TargetStims' filesep];
        dir_noise  = [dir_subj 'NoiseStims'  filesep];

        % ListSignals(2).name = 'target.wav';ListSignals(1).name = 'nontarget.wav';
        fname_targets = {'target.wav','nontarget.wav'};
        
        ListSignals = [];
        for i = 1:length(fname_targets)
            ListSignals(i).name = fname_targets{i};
            
            %signal = create_AM(cfg_game.fc, cfg_game.fm, 10^(mean(m)/10)*(i-1), cfg_game.stim_dur, cfg_game.fs)';
            signal = audioread([dir_target fname_targets{i}]);
            lvl_signal_before = rmsdb(signal)+dBFS;
            stim_normal = scaletodbspl(signal,lvl,dBFS); % stim_normal = dBlvl(signal,cfg_game.SPL);
            lvl_signal_after = rmsdb(stim_normal)+dBFS;
            
            %%% New
            lvl_diffe = lvl_signal_after - lvl_signal_before;
            if lvl_diffe > 1
                % Stores the normalised signal:
                audiowrite([dir_target fname_targets{i}], stim_normal, fs);
            else
                fprintf('%s: the level difference before/after calibration is too small (delta=%.1f dB), no new files will be stored\n', ...
                    upper(mfilename),lvl_diffe);
            end
            %%%
        end

        if bLeo
            method = 'butter';
        end
        if bAlejandro
            method = 'Gammatone_proc';
        end

        opts = [];
        opts.method = method;
        opts.dBFS = dBFS;
        opts.SNR  = SNR;
        opts.lvl  = lvl;
        opts.basef= basef;
        idx_order = data_passation.n_stim;
        noise_E = il_noisetone_converter(dir_target,dir_noise, ListStim, idx_order, ...
                fcut, undersampling, fcut_noiseE, opts);

        % noise_E = 1x360x3000

        disp('')
 
        % Analyse targets and compute ideal template
        idx_order = [1 2];
        [signal_E] = il_noisetone_converter([],dir_target,ListSignals,idx_order,fcut,undersampling,fcut_noiseE,opts);
        ideal_template = signal_E(:,:,2)-signal_E(:,:,1);

        % Ideal template in noise: Deleted but you can find it back in the original Script3_AnalysisComplex.m
        %% Compute CI and target-present and target-absent sub-CI
 
        n_rand = 200;
        n_boot = 200;
        n_trials = length(n_response);

        %subselection of trials
        trials2analyze = 1:n_trials;
        is_norm = 0;
        [CI, CIrand, CIboot, ResponseMatrix, CI2, CI2rand, CI2boot, CI1, CI1rand, CI1boot] = computeCI(n_signal(trials2analyze),n_response(trials2analyze), noise_E(:,:,trials2analyze), n_rand, n_boot, is_norm);

        tE=(1:size(CI,2))/(cfg.fs/cfg.save_undersmpl);

        CIrand_ci = [prctile(CIrand,CIlim(2),3); prctile(CIrand,CIlim(1),3)];
        CI1rand_ci = [prctile(CI1rand,CIlim(2),3); prctile(CI1rand,CIlim(1),3)];
        CI2rand_ci = [prctile(CI2rand,CIlim(2),3); prctile(CI2rand,CIlim(1),3)];
        save([dir_subj 'CIt'],'CI','CI1','CI2','CIrand_ci','CI1rand_ci','CI2rand_ci','tE','ideal_template')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function undersmplE = il_noisetone_converter(dir_target,dir_noise,ListStim,n_stim,fcut,undersampling,fEcut,opts)
% New extra input parameter: method

method = opts.method;
dBFS   = opts.dBFS; 
SNR    = opts.SNR; % dB
lvl    = opts.lvl; % dB
Nchannel = length(fcut)-1;

basef = opts.basef; % Hz

% ATTENTION BRICOLAGE %
if ~isempty(dir_target)
    [tone, fs] = audioread([dir_target 'nontarget.wav']);
else
    tone = [];
end

for i_trial=1:length(n_stim)
    fprintf(['stim #' num2str(i_trial) '\n'])
    [noise, fs] = audioread([dir_noise ListStim(n_stim(i_trial)).name ]);
    
    % ATTENTION BRICOLAGE %
    if ~isempty(tone)
        %   Creating the trial (but always using the non-modulated sound)
        S = il_Addition_RSB(tone,noise,SNR);
        S = scaletodbspl(S,lvl,dBFS); % same as dBlvl(S,lvl) with 72.6 as lvl_ref
    else
        if i_trial == 1
            warning('Alejandro: I do not agree')
        end
        S=noise/rms(noise);
    end
    
    for i_channel = 1:Nchannel
        switch method
            case 'butter'
                % Filter design:
                [B, A]   = butter(2, fcut(i_channel:i_channel+1)/(fs/2));
                % Applying the filter:
                Sout = abs(filterfilter(B, A, S)); % absolute value is full wave rectification
                % Sout = max(filterfilter(B, A, S),0); % This would be half-wave rectification
                if ~isempty(tone)
                     % Regular processing:
                     [B2, A2] = butter(1, fEcut/(fs/2));
                     E = filterfilter(B2, A2, Sout);
                else
                    if i_trial == 1 && i_channel == 1
                        warning('Alejandro: I do not agree')
                        E = Sout;
                    end
                end
                E = E(1:undersampling:end);
                 
             case 'Gammatone_proc'
                 binwidth = undersampling/fs; % 'undersampling' samples
                 flags_gamma = {'basef',basef,'flow',fcut(i_channel),'fhigh',fcut(i_channel+1), ...
                     'bwmul',1,'dboffset',dBFS,'no_adt','binwidth',binwidth, ...
                    'no_outerear','no_middleear', 'hilbert'};
                 [E, fc, t, outs] = Gammatone_proc(S,fs,flags_gamma{:});
         end
         undersmplE(i_channel, :, i_trial) = E;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ somme, A ] = il_Addition_RSB(S, B, SNR)
% [ somme, A ] = Addition_RSB(S, B, RSB) 
% additionne un son S à un bruit B avec un RSB donné (en dB), et renvoie 
% le facteur A tel que somme = A*S + B;

Ps = mean(S.^2);
Pb = mean(B.^2);

A=sqrt((Pb/Ps)*10^(SNR/10));
somme = A*S + B;