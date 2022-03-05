function publ_varnet2022a_utils(Subject_ID,type_action,flags,keyvals)
% function publ_varnet2022a_utils(Subject_ID,type_action,flags,keyvals)
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

% cfg_game = []; % : [1×1 struct]
% data_passation = []; % : [1×1 struct]
% load(savefilename);
         
% ListStim: [3000×1 struct]
[cfg_game, data_passation, ListStim] = Convert_ACI_data_type(savefilename);
m = data_passation.expvar;
n_signal = data_passation.n_targets; % (old field data_passation.N_signal) 
                                      % 1 a target was presented, 2 a nontarget was presented
n_response = data_passation.n_responses;

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
        bLeo = flags.do_varnet2022a_JASA;
        bAlejandro = 0;
        bSFA = flags.do_varnet2022b_CFA;
        
        % Analyse noise in bands
        basef = 1000;
        BW_ERB = 1;
        undersampling = 100;
        if bLeo
            fcut = ERB2f(f2ERB(basef)+BW_ERB/2*[-1 1]);
        end
        if bAlejandro
            fcut = audtofreq(freqtoaud(basef)+BW_ERB/2*[-1 1]); % here there is a slight difference (ERB2f)
        end
        if bSFA
            fcut = [40 8000]; % Hz
        end

        Nchannel = length(fcut)-1; 
        
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
            bCalculate = 1;
        end
        if bAlejandro
            method = 'Gammatone_proc';
            bCalculate = 1;
        end
        if bSFA
            % method = 'Gammatone_proc2';
            method = 'Gammatone_proc'; warning('Temporarily changed to Gammatome_proc')
        end

        if bSFA
            switch method
                case 'Gammatone_proc'
                    method_here = 'gammatone';
                case 'Gammatone_proc2'
                    method_here = 'gammatone2'; % I added this flag in arg_fastACI_getACI.m
            end
            % glmfct = 'classic_revcorr';
            glmfct = 'l1glm';
            
            flags_here = {glmfct, ...
                'expvar_limits',[], ...
                'trialtype_analysis','total', ...
                'consistency_check',0, ... % new option to disable the consistency check
                method_here};
            
            switch glmfct
                case 'l1glm' % 'lassoglmslow'
                    N_lambda = 30;
                    Lambdas = logspace(-4, -1, N_lambda);
                    idx = find(Lambdas >= 10^-3);
                    Lambdas = Lambdas(idx);

                    flags_here{end+1} = 'lambda';
                    flags_here{end+1} = Lambdas;
            end
            if isfield(keyvals,'dir_out')
                if exist(keyvals.dir_out,'dir')
                    flags_here{end+1} = 'dir_out';
                    flags_here{end+1} = keyvals.dir_out;
                end
            end
            fnameACI = fastACI_getACI_fname(savefilename,flags_here{:});
            
            if exist(fnameACI,'file')
                bCalculate = 0;
            end
        end
        
        opts = [];
        opts.method = method;
        opts.dBFS = dBFS;
        opts.SNR  = SNR;
        opts.lvl  = lvl;
        opts.basef= basef;

        idx_order = cfg_game.stim_order; % data_passation.n_stim; % should be the same, remove n_stim
      
        if bCalculate
            noise_E = il_noisetone_converter(dir_target,dir_noise, ListStim, idx_order, ...
                    fcut, undersampling, fcut_noiseE, opts);
            [au,fsau] = audioread([fastACI_basepath 'Stimuli' filesep 'ready.wav']);
            sound(au,fsau);
            
            % Data_matrix = permute(noise_E, [3 1 2]);
            Data_matrix = permute(noise_E, [3 2 1]);
            % noise_E     = 1x360x3000
            % Data_matrix = 3000x1x360
            
            flags_here{end+1} = 'Data_matrix';
            flags_here{end+1} = Data_matrix;
        end

        % Analyse targets and compute ideal template
        idx_order = [1 2];
        [signal_E, fc, t] = il_noisetone_converter([],dir_target,ListSignals,idx_order,fcut,undersampling,fcut_noiseE,opts);
        ideal_template = signal_E(:,:,1)-signal_E(:,:,2);

        % Ideal template in noise: Deleted but you can find it back in the original Script3_AnalysisComplex.m
        %% Compute CI and target-present and target-absent sub-CI
        if bLeo || bAlejandro
            n_rand = 200;
            n_boot = 200;
            n_trials = length(n_response);

            %subselection of trials
            trials2analyse = 1:n_trials;
            is_norm = 0;
            [CI, CIrand, CIboot, ResponseMatrix, CI2, CI2rand, CI2boot, CI1, CI1rand, CI1boot] = computeCI(n_signal(trials2analyse),n_response(trials2analyse), noise_E(:,:,trials2analyse), n_rand, n_boot, is_norm);

            tE=t;%(1:size(CI,2))/(cfg_game.fs/undersampling);

            CIlim = [2.5 97.5];
            CIrand_ci = [prctile(CIrand,CIlim(2),3); prctile(CIrand,CIlim(1),3)];
            CI1rand_ci = [prctile(CI1rand,CIlim(2),3); prctile(CI1rand,CIlim(1),3)];
            CI2rand_ci = [prctile(CI2rand,CIlim(2),3); prctile(CI2rand,CIlim(1),3)];
            save([dir_subj 'CIt'],'CI','CI1','CI2','CIrand_ci','CI1rand_ci','CI2rand_ci','tE','ideal_template')
        end
         
         if bSFA
            [ACI,cfg_ACI,results] = fastACI_getACI(savefilename,flags_here{:});
            
            XT  = get(gca,'XTick');
            XTL = round(100*t(XT))/100;
            set(gca,'XTickLabel',XTL);
            xlabel('Time (s)');
            
            YT = get(gca,'YTick');
            YTL = round(fc(round(YT)));
            set(gca,'YTickLabel',YTL);
            ylabel('Frequency (Hz)');
            
            if flags.do_SA
                title('SA (S4 in varnet2022a, S-LV)')
            end
            if flags.do_SB
                title('SB (S1 in varnet2022a, S-AO)')
            end
            
            disp('')
        end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [undersmplE, fc, t] = il_noisetone_converter(dir_target,dir_noise,ListStim,n_stim,fcut,undersampling,fEcut,opts)
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
                undersmplE(i_channel, :, i_trial) = E;

                t=(1:size(E,2))/(fs/undersampling);
                fc = fcut;
            case 'Gammatone_proc'
                 binwidth = undersampling/fs; % 'undersampling' samples
                 flags_gamma = {'basef',basef,'flow',fcut(i_channel),'fhigh',fcut(i_channel+1), ...
                     'bwmul',.5,'dboffset',dBFS,'no_adt','binwidth',binwidth, ...
                    'no_outerear','no_middleear', 'hilbert'};
                 [E, fc, t, outs] = Gammatone_proc(S,fs,flags_gamma{:});
                 undersmplE(:, :, i_trial) = E;             
            case 'Gammatone_proc2'
                 binwidth = undersampling/fs; % 'undersampling' samples
                 flags_gamma = {'basef',basef,'flow',700,'fhigh',1200, ...%50,'fhigh',8000, ...
                     'bwmul',1,'dboffset',dBFS,'no_adt','binwidth',binwidth, ...
                     'no_outerear','no_middleear', 'hilbert'};
                 [E, fc, t, outs] = Gammatone_proc(S,fs,flags_gamma{:});
                 % Leo: I applied the transposition to get the same 'format'
                 %   as in 'Gammatone_proc' and 'butter', which is:
                 %   Dim_freq x Dim_time x Dim_trial. We need to swap the 
                 %   dimensions later to get Dim_trial x Dim_freq x Dim_trial,
                 %   compatibly with Data_matrix...
                 undersmplE(:, :, i_trial) = transpose(E); 
        end
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