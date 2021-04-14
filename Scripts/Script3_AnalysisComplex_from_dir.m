function data = Script3_AnalysisComplex_from_dir(dir_main,Subject_ID,opts)
% function data = Script3_AnalysisComplex_from_dir(dir_data,Subject_ID,opts)
%
% Abbreviations in this script (alphabetic order):
%   CR: Correct rejection
%   FA: False alarm
%   H : Hit
%   PC: Percentage correct
%   RT: Reaction time
%
% Notes for Leo on 4/02/2021: 
%   - histogram does not provide edges (please confirm I am correct) but only
%     the centre of the bins. I added a function called my_hist.m which is 
%     exactly as the MATLAB function hist.m, were as extra output I added 
%     the real edges.
%
%   - I added "greater or equal than" to the line 'idx_m = find(m>=m_edge(i_m) & m<m_edge(i_m+1));'
%     Otherwise some of the m values exactly on the edge were not being counted.
% 
% data_passation.resume_trial: Trial number at which a new session is started
%
% % Stand-alone example:
%       % Enter the full path to your data, on Unix, e.g.:
%       dir_data = '/home/alejandro/Downloads/S_LV/';
%       Script3_AnalysisComplex(dir_data);
%
%       % On Windows it could be something like:
%       dir_data = 'H:\Downloads\S_LV\';
%       Script3_AnalysisComplex(dir_data);
%
%       % If your data is the current directory:
%       dir_data = [pwd filesep];
%       Script3_AnalysisComplex(dir_data);
%
% Other scripts from where this processing is called:
%       g20210211_learning_revcorr;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data = [];

if nargin == 0
    close all

    dir_main = '/home/alejandro/Downloads/';
    Subject_ID = 'S_LV';
end

if nargin < 3
    opts = [];
end
opts = Ensure_field(opts,'do_behaviour',1);
opts = Ensure_field(opts,'do_analyse_noise',0);

dir_data = [dir_main Subject_ID filesep];
    
if ~exist(dir_data,'dir')
    help Script3_AnalysisComplex
    error('%s: Please specify (manually) a directory containing experimental data',upper(mfilename));
end

%%%
% [path,name,ext]=fileparts(which(mfilename)); % path will be the folder where this file is located...
% dir_main = [path filesep]; 
%%%

%%% Variables for storing figures:
% h = [];     % empty handles
% hname = []; % empty figure names 
% dir_out = dir_data;
%%%

% switch Subject_ID
%     case 'S_LV'
%         do_behaviour = 1; % Set to one to analyse/plot the behavioural results
%         do_analyse_noise = 0;
%         disp('For Subject S_LV I do not have the sound stimuli locally, running behavioural results only...');
%     case 'S_VR'
%         do_behaviour = input('do_behaviour (1=yes; 0=no): ');
%         do_analyse_noise =  input('do_analyse_noise (1=yes; 0=no): ');
%     otherwise
%         error('%s: Participant not validated yet, add manually to the script''s list and re-run...')
% end
% 
% dBFS = 93.6139; % previous knowledge

% -------------------------------------------------------------------------
% --- Load data:
files = Get_filenames(dir_data,'savegame_*.mat');
Show_cell(files);
if ~isempty(files)
    if length(files) == 1
        idx2proc = 1;
    else
        idx2proc = input(['Choose the file to process (enter 1-' num2str(length(files)) '): ']);
    end
else
    error('%s: No files to process were found. Please provide another ''dir_data'' folder and re-run this script',upper(mfilename))
end

file_savegame = [dir_data files{idx2proc}];

data = Script3_AnalysisComplex(file_savegame,opts);

% var = load(file_savegame);
% data_passation_old = var.data_passation;
% cfg_game_old   = var.cfg_game;
% ListStim_old   = var.ListStim; % used in do_analyse_noise
% 
% [cfg_game, data_passation, ListStim] = Convert_ACI_data_type(file_savegame);
% 
% n_responses = data_passation.n_responses;
% N_trials    = length(n_responses); % completed number of trials
% n_signal    = data_passation.n_targets(1:N_trials);
% % -------------------------------------------------------------------------
% if do_behaviour
%     if isfield(cfg_game,'N')
%         N_total = cfg_game.N; 
%     else
%         N_total = cfg_game.N_presentation*cfg_game.N_target; 
%     end
%     
%     if N_total ~= N_trials
%         disp('This participant did not complete the experiment...')
%         N_trials = floor(N_trials/100)*100;
%         fprintf('\tN_trials is being truncated to a multiple of 100 trials: %.0f\n',N_trials);
%     end
% 
%     expvar     = data_passation.expvar(1:N_trials);
%     is_correct = data_passation.is_correct;
%     RT         = data_passation.response_time(1:N_trials);
% 
%     resp_if_tar = 2; % 'if target' (AM)       - 2 modulated tone
%     resp_if_ref = 1; % 'if reference' (no AM) - 1 pure tone
% 
%     % -------------------------------------------------------------------------
%     % --- First figure: Behavioural results
%     n_window = 100;
%     minNformean = 50;
% 
%     m_bin = -15.5:1:-0.5;
%     % [N_m,m_edge] = histcounts(m, nbins); % Only in newer MATLAB versions
%     % N_m = hist(m,m_bin); m_edge = m_bin; % In older MATLAB versions
%     [N_m,m_edge] = my_hist(expvar,m_bin); % this really provides edges
% 
%     for i_m = 1:length(m_bin)
%         idx_m   = find(expvar>=m_edge(i_m) & expvar<m_edge(i_m+1));
%         idxs_m(i_m) = length(idx_m);
%         
%         n_is_TAR(i_m)  = sum(n_signal(idx_m)==resp_if_tar); 
%         n_is_REF(i_m)  = sum(n_signal(idx_m)==resp_if_ref);
%         
%         H(i_m)  = sum(n_signal(idx_m)==resp_if_tar & n_responses(idx_m)==resp_if_tar); % Hit
%         M(i_m)  = sum(n_signal(idx_m)==resp_if_tar & n_responses(idx_m)==resp_if_ref); % Miss
%         CR(i_m) = sum(n_signal(idx_m)==resp_if_ref & n_responses(idx_m)==resp_if_ref); % Correct rejection
%         FA(i_m) = sum(n_signal(idx_m)==resp_if_ref & n_responses(idx_m)==resp_if_tar); % False alarm
%     end
%     H_tot  = sum(H); 
%     M_tot  = sum(M);
%     CR_tot = sum(CR);
%     FA_tot = sum(FA);
%     tot_classified    = H_tot + CR_tot + M_tot + FA_tot;
%     counted_responses = sum(idxs_m);
% 
%     % Sanity check (to know that all responses were processed):
%     if counted_responses ~= N_total
%         error('%s: Not all responses were processed. Check whether the participant indeed completed the whole session...',upper(mfilename))
%     end
%     if tot_classified ~= N_total
%         error('%s: Not all responses were classified as H, M, CR, or FA. Check whether this is correct. If yes, convert this message into a warning only...',upper(mfilename))
%     end
% 
%     N_windows = length(n_responses)/n_window;
% 
%     % Memory allocation:
%     m_windowed       = nan(1,N_windows);
%     bias_windowed    = nan(1,N_windows);
%     PC_targetpresent = nan(1,N_windows);
%     PC_targetabsent  = nan(1,N_windows);
%     RT_windowed      = nan(1,N_windows);
% 
%     for i = 1:N_windows
%         idxs_here = (i-1)*n_window+1:i*n_window; % indexes of the trials within each window
%         response_windowed  = n_responses(idxs_here);
%         signal_windowed    = n_signal(idxs_here);
% 
%         m_windowed(i)      = mean(expvar(idxs_here));
%         bias_windowed(i)   = mean(response_windowed);
%         PC_targetpresent(i)= mean(response_windowed(signal_windowed==2))-1;
%         PC_targetabsent(i) = 2-mean(response_windowed(signal_windowed==1));
%         RT_windowed(i)     = mean(RT(idxs_here));
%     end
% 
%     % ---
%     m_presentation_nr     = 1:N_trials;
%     m_presentation_nr_win = 1:n_window:N_trials;
% 
%     figure('Position', [100 100 800 500]); 
%     subplot(2,2,1); 
%     plot(m_presentation_nr    , expvar         ,'g'); hold on; 
%     plot(m_presentation_nr_win, m_windowed,'k'); 
%     ylim([m_edge(1) m_edge(end)]); 
%     xlabel(' trial #'); ylabel('m (dB)'); 
%     xlim([1 length(expvar)]); ylimits=ylim;
% 
%     N_sessions = length(data_passation.resume_trial);
%     for i = 1:N_sessions
%         plot(data_passation.resume_trial(i)*[1 1],ylimits,'k:');
%     end
% 
%     % ---
%     subplot(2,2,3); 
%     plot(m_presentation_nr_win,PC_targetpresent); hold on; 
%     plot(m_presentation_nr_win,PC_targetabsent);  
%     xlim([1 length(expvar)]); xlabel(' trial #'); 
%     ylabel('correct response rate'); ylim([0 1]); hold on; 
%     plot([1 length(expvar)],[0.5 0.5],'k--'); 
%     ylimits=ylim;
% 
%     for i = 1:length(data_passation.resume_trial)
%         % Vertical dotted lines at the points where a new session was started:
%         plot(data_passation.resume_trial(i)*[1 1],ylimits,'k:','LineWidth',2);
%     end
% 
%     subplot(2,2,2); 
%     bar(m_bin, [H', M', CR', FA']); 
%     xlim([m_edge(1) m_edge(end)]); 
%     xlabel('m (dB)'); 
%     ylabel('Nb of trials'); hold on; 
%     plot([m_edge(1) m_edge(end)],[minNformean minNformean]/2,'k:');
%     legend({'H', 'M', 'CR', 'FA', 'Nmin'});
% 
%     subplot(2,2,4); 
%     bar(m_bin, [H'./(M'+H'), CR'./(CR'+FA')].*[M'+H'>minNformean,CR'+FA'>minNformean]);
%     xlim([m_edge(1) m_edge(end)]); 
%     xlabel('m (dB)'); 
%     ylabel('correct response rate'); 
%     hold on; 
%     plot([m_edge(1) m_edge(end)],[0.5 0.5],'k--'); 
%     legend({'target present', 'target absent','chance level'},'Location','southeast');
% 
%     h(end+1) = gcf;
%     hname{end+1} = [Subject_ID '-Behaviour'];
% 
%     %%%
%     H_rate = 100*H./n_is_TAR; % /N_tot_signal;
%     CR_rate = 100*CR./n_is_REF; % /N_tot_absent;
%     
%     figure;
%     plot(m_bin,H_rate,'bo-'); hold on;
%     plot(m_bin,CR_rate,'rs--','LineWidth',2);
%     ylim([-3 103]);
%     set(gca,'YTick',0:5:100); grid on
%     xlim([-17 1])
%     set(gca,'XTick',-16:0);
%     
%     legend('H rate','CR rate','Location','NorthWest');
%     ylabel(sprintf('Percentage correct\n(ref. # of presentations per m interval)'));
%     xlabel('Central bin of the tested modulation depths (dB)');
%     
%     h(end+1) = gcf;
%     hname{end+1} = [Subject_ID '-Hit-and-CR-rates-per-bin'];
%     
%     %%%
%     
%     trialnum = 1:n_window:length(expvar);
%     
%     data.trialnum = trialnum;
%     data.m_windowed = m_windowed;
%     data.PC_targetpresent = PC_targetpresent;
%     data.PC_targetabsent  = PC_targetabsent;
%     data.bias_windowed    = bias_windowed;
%     data.N_m = N_m;
%     data.m_edge = m_edge;
%     data.H = H;
%     data.M = M;
%     data.CR = CR;
%     data.FA = FA;
%     
%     data.h = h;
%     data.hname = hname;
%     
%     if nargout == 0
%         f2store = [dir_data 'Behaviour.mat'];
% 
%         if exist(f2store,'file')
%             bSave = input('File already exists, type 1 to overwrite, 0 to skip this step: ');
%         else
%             bSave = 1;
%         end
% 
%         save(f2store,'m_windowed','PC_targetpresent','PC_targetabsent','bias_windowed','trialnum','N_m','m_edge', 'H', 'M', 'CR', 'FA')
%     end
% end
% 
% bSave = 0;
% % -------------------------------------------------------------------------
% % --- Analyse noise in bands
% if do_analyse_noise
%     %%% My notes:
%     %    - noisetone_converter.m and noise_converter.m apply the same processing
%     %        but for the non-AM tone summed to the noises or for non-AM and AM tones,
%     %        respectively. It would be better to have the same function where you
%     %        use as inputs the waveform you want to extract the envelope from...
%     %    - The ideal template (ideal_template) is just the difference in envelopes
%     %        between the non-AM and AM tones.
%     
%     cfg = []; % structure to store parameters
%     %%% Plot options:
%     CI_lo =  2.5; % 'confidence intervals' lower limit
%     CI_up = 97.5; % 'confidence intervals' upper limit
%     %%% Experiment depentendent parameters:
%     lvl = 65;
%     SNR = -10;
%     fc  = 1000;
%     %%%
%     
%     BWerb = 2;
%     fcut = il_audtofreq(il_freqtoaud(fc)+[-BWerb/2 BWerb/2]);
%     Nchannel = length(fcut)-1; 
%     
%     undersampling = 100;% 10
%     fcut_noiseE   = 10; % 15;%30;%
%     foldername    = cfg_game.folder_name;
%     
%     dir_where = dir_data; warning('AO: Temporally using another folder...')
%     bLoad = 1;
%     
%     cfg.fs = cfg_game.fs;
%     cfg.save_undersmpl = undersampling;
%     
%     noise_E = [];
%     
%     % ---------------------------------------------------------------------
%     % Obtains the envelope of the tone+noise intervals (no AM in the tone) by:
%     %     bandpass filtering +/- 1 ERB_N around fc, second order Butterworth 
%     %     bandpass (i.e., IIR filter of order 4). The envelope is obtained from
%     %     lowpass filtering (fcut_noiseE=10 Hz, the defaul) the absolute value
%     %     of the filtered noise. The noise is finally downsampled.
%     if bLoad
%         warning('I skipped the use of noisetone_converter...I ran it once and I stored the noises...')
%         load('/home/alejandro/Downloads/S_LV/tmp_vars/noise_E.mat');
%     else
%         [noise_E] = noisetone_converter(dir_where, foldername, ListStim, data_passation.n_stim, fcut, undersampling, fcut_noiseE,lvl,SNR);
%         %[noise_E] = noise_converter(foldername, ListStim, data_passation.n_stim, fcut, undersampling, fcut_noiseE,lvl,SNR);
%         % % noise_E is 1 x 360 x 3000
%     end
%     % ---------------------------------------------------------------------
%     
%     % Analyse targets and compute ideal template
% 
%     foldername = 'TargetStims';
%     ListSignals(2).name = 'target.wav';
%     ListSignals(1).name = 'nontarget.wav';
%     
%     % ---------------------------------------------------------------------
%     % Sets the levels if needed (this won't be needed anymore with Alejandro's
%     % approach):
%     for i = 1:2
%         fname = ListSignals(i).name;
%         %signal = create_AM(cfg_game.fc, cfg_game.fm, 10^(mean(m)/10)*(i-1), cfg_game.stim_dur, cfg_game.fs)';
%         signal = audioread([dir_where foldername filesep fname]);
%         stim_normal = dBlvl(signal,cfg_game.SPL);
% 
%         lvl_signal = rmsdb(signal)+dBFS;
%         file2store = [dir_where foldername filesep fname];
%         if exist(file2store,'file')
%             fprintf('The signal already exists and its level is %.1f dB (should be %.1f dB), if incorrect, check your local files\n',lvl_signal,cfg_game.SPL);
%         else
%             audiowrite(file2store, stim_normal, cfg_game.fs);
%         end
%     end
% 
%     [signal_E] = noise_converter(dir_where,foldername,ListSignals,[1,2],fcut,undersampling,fcut_noiseE);
%     ideal_template = signal_E(:,:,2)-signal_E(:,:,1);
% 
%     %% Compute CI and target-present and target-absent sub-CI
% 
%     n_rand = 200;
%     n_boot = 200;
%     n_trials = length(n_responses);
% 
%     % Subselection of trials
%     trials2analyse = 1:n_trials; %find(m>=prctile(m,10) & m<=prctile(m,90));%%%%(iscorrect(i_trial)==0);%responsetime(i_trial)>=median(responsetime);%
%     is_normalised = 0; % 'no'
%     
%     %[CI, CIrand, CIboot]  = computeCI(n_signal(trials2analyze),n_response(trials2analyze), noise_E(:,:,trials2analyze), n_rand, n_boot, 'yes');
%     [CI , CIrand , CIboot , ResponseMatrix, ...
%      CI2, CI2rand, CI2boot, ...
%      CI1, CI1rand, CI1boot] = computeCI(n_signal(trials2analyse),n_responses(trials2analyse), noise_E(:,:,trials2analyse), n_rand, n_boot, is_normalised);
% 
%     tE=(1:size(CI,2))/(cfg.fs/cfg.save_undersmpl);
% 
%     dim = 3;
%     CIrand_ci  = [prctile(CIrand ,CI_up,dim); prctile(CIrand ,CI_lo,dim)];
%     CI1rand_ci = [prctile(CI1rand,CI_up,dim); prctile(CI1rand,CI_lo,dim)];
%     CI2rand_ci = [prctile(CI2rand,CI_up,dim); prctile(CI2rand,CI_lo,dim)];
% 
%     if bSave
%         save([dir_out 'CIt.mat'],'CI','CI1','CI2','CIrand_ci','CI1rand_ci','CI2rand_ci','tE','ideal_template')
%     end
% 
%     %%% Preparing further plots
%     Nchannels2plot      = 1:Nchannel;
%     
%     %%% Plot temporal revcorr
%     % 1. General kernel:
%     %    I do not understand the scaling by 'rms(ans(CI))', probably is for 
%     %    visual purposes only:
%     ideal_templ2plot = rms(abs(CI(:)))*ideal_template'/max(max(ideal_template)); 
%     
%     figure('Name', 'General kernel'); hold on
%     plot_channels(tE, ideal_templ2plot, Nchannels2plot, @(x,y)(plot(x,y,'r-')))
%     plot_channels(tE, zeros(size(CI')), Nchannels2plot, @(x,y)(plot(x,y,'k-')));
%     plot_channels(tE, CI'             , Nchannels2plot, @(x,y)(plot(x,y,'b-'))); 
%     
%     if ~isempty(CIrand) 
%         plot_channels(tE, prctile(CIrand,CI_lo,dim)',Nchannels2plot, @(x,y)(plot(x,y,'k:')))%, @plot, 1, size(undersmplE,1))
%         plot_channels(tE, prctile(CIrand,CI_up,dim)',Nchannels2plot, @(x,y)(plot(x,y,'k:')))%, @plot, 1, size(undersmplE,1))
%     end
%     if ~isempty(CIboot) 
%         plot_channels(tE, prctile(CIboot,CI_lo,dim)',Nchannels2plot, @(x,y)(plot(x,y,'b:')))%, @plot, 1, size(undersmplE,1))
%         plot_channels(tE, prctile(CIboot,CI_up,dim)',Nchannels2plot, @(x,y)(plot(x,y,'b:')))%, @plot, 1, size(undersmplE,1))
%     end
%     set(findobj(gcf,'type','axes'),'YLim',[-1 1]*1.1*max(abs([CI(:); CI1(:); CI2(:)])));
% 
%     if bSave
%         saveas(gcf,[dir_out 'CIt_all.png']);
%     end
% 
%     % 2. Target-absent kernel
%     ideal_templ2plot = rms(abs(CI1(:)))*ideal_template'/max(max(ideal_template));
%     figure('Name', 'Target-absent kernel');
%     plot_channels(tE, ideal_templ2plot , Nchannels2plot, @(x,y)(plot(x,y,'r-')))%, @plot, 1, size(undersmplE,1))
%     plot_channels(tE, CI1'             , Nchannels2plot, @(x,y)(plot(x,y,'b-')))%, @plot, 1, size(undersmplE,1))
%     plot_channels(tE, zeros(size(CI1')), Nchannels2plot, @(x,y)(plot(x,y,'k-')));
%     
%     if ~isempty(CIrand) 
%         plot_channels(tE, prctile(CI1rand,CI_up,dim)', Nchannels2plot, @(x,y)(plot(x,y,'k:')))%, @plot, 1, size(undersmplE,1))
%         plot_channels(tE, prctile(CI1rand,CI_lo,dim)', Nchannels2plot, @(x,y)(plot(x,y,'k:')))%, @plot, 1, size(undersmplE,1))
%     end
%     if ~isempty(CIboot) 
%         plot_channels(tE, prctile(CI1boot,CI_lo,dim)', Nchannels2plot, @(x,y)(plot(x,y,'b:')))%, @plot, 1, size(undersmplE,1))
%         plot_channels(tE, prctile(CI1boot,CI_up,dim)', Nchannels2plot, @(x,y)(plot(x,y,'b:')))%, @plot, 1, size(undersmplE,1))
%     end
%     set(findobj(gcf,'type','axes'),'YLim',[-1 1]*1.1*max(abs([CI(:); CI1(:); CI2(:)])));
%     
%     if bSave
%         saveas(gcf,[dir_out 'CIt_ta.png'])
%     end
% 
%     ideal_templ2plot = rms(abs(CI2(:)))*ideal_template'/max(max(ideal_template)); 
%     figure('Name', 'Target-present kernel');
%     plot_channels(tE, ideal_templ2plot , Nchannels2plot, @(x,y)(plot(x,y,'r-')))%, @plot, 1, size(undersmplE,1))
%     plot_channels(tE, zeros(size(CI2')), Nchannels2plot, @(x,y)(plot(x,y,'k-')));
%     plot_channels(tE, CI2'             , Nchannels2plot, @(x,y)(plot(x,y,'b-')))%, @plot, 1, size(undersmplE,1))
%     if ~isempty(CIrand)
%         plot_channels(tE, prctile(CI2rand,CI_up,dim)',Nchannels2plot, @(x,y)(plot(x,y,'k:')))%, @plot, 1, size(undersmplE,1))
%         plot_channels(tE, prctile(CI2rand,CI_lo,dim)',Nchannels2plot, @(x,y)(plot(x,y,'k:')))%, @plot, 1, size(undersmplE,1))
%     end
%     if ~isempty(CIboot)
%         plot_channels(tE, prctile(CI2boot,CI_lo,dim)',Nchannels2plot, @(x,y)(plot(x,y,'b:')))%, @plot, 1, size(undersmplE,1))
%         plot_channels(tE, prctile(CI2boot,CI_up,dim)',Nchannels2plot, @(x,y)(plot(x,y,'b:')))%, @plot, 1, size(undersmplE,1))
%     end
%     set(findobj(gcf,'type','axes'),'YLim',[-1 1]*1.1*max(abs([CI(:); CI1(:); CI2(:)])));
%     
%     if bSave
%         saveas(gcf,[dir_out 'CIt_tp.png'])
%     end
% 
%     fE=sqrt(fcut(1:end-1).*fcut(2:end));
% 
%     figure('Name', 'General kernel');
%     h=pcolor(tE,fcut,[CI; 1e-5*ones(1,360)]); set(h, 'edgecolor','none'); 
%     colorbar; 
%     C_axis = caxis; 
%     caxis([-1 1]*max(abs(C_axis))); hold on
%     plot(tE,50*ideal_template+fE','r'); 
%     
%     if bSave
%         saveas(gcf,[dir_out 'CItf_all.png'])
%     end
% 
%     figure('Name', 'Target-absent kernel');
%     h=pcolor(tE,fcut,[CI1; 1e-5*ones(1,360)]); set(h, 'edgecolor','none'); colorbar; caxis([-1 1]*max(abs(C_axis)))
%     hold on
%     plot(tE,50*ideal_template+fE','r');
%     
%     if bSave
%         saveas(gcf,[dir_out 'CItf_ta.png'])
%     end
% 
%     figure('Name', 'Target-present kernel');
%     h=pcolor(tE,fcut,[CI2; 1e-5*ones(1,360)]); set(h, 'edgecolor','none'); colorbar; caxis([-1 1]*max(abs(C_axis)))
%     hold on
%     plot(tE,50*ideal_template+fE','r');
%     
%     if bSave
%         saveas(gcf,[dir_out 'CItf_tp.png'])
%     end
% 
%     % ---------------------------------------------------------------------
%     %% complex fft of revcorr
%     bComplex_FFT = 1;
%     bComplex_FFT_plot = 0;
%     bExtract_metrics_using_Fourier = 0;
%     bExtract_metrics_using_xcorr   = 0; 
%     bComplex_FFT_revcorr = 0;
%     
%     if bComplex_FFT
%         clear CIfft CI1fft CI2fft ideal_templatefft ideal_templatecfft CIfftrand CI1fftrand CI2fftrand CIfftboot CI1fftboot CI2fftboot
%         Nfft = 512*100;
%         xl = [0 fcut_noiseE];
%         fE = (1:Nfft)/Nfft*(cfg.fs/cfg.save_undersmpl);
%         fidx2plot = fE>=xl(1) & fE<=xl(2);
%         CIfftrand = zeros(size(noise_E,1),Nfft,n_rand);
% 
%         for i_channel = 1:size(noise_E,1)
%             CIfft(i_channel,:) = (fft(CI(i_channel,:),Nfft));
%             CI1fft(i_channel,:) = (fft(CI1(i_channel,:),Nfft));
%             CI2fft(i_channel,:) = (fft(CI2(i_channel,:),Nfft));
%             ideal_templatecfft(i_channel,:) = (fft(ideal_template(i_channel,:),Nfft));
%             if n_rand > 0
%                 for i_rand = 1:n_rand
%                     CIfftrand(i_channel,:,i_rand) = (fft(CIrand(i_channel,:,i_rand),Nfft));%CIrand(i_channel,:,i_rand)+Enorm(i_channel,:,i_trial)*correlator*keeptrial/n_trials;
%                     CI1fftrand(i_channel,:,i_rand) = (fft(CI1rand(i_channel,:,i_rand),Nfft));%CIrand(i_channel,:,i_rand)+Enorm(i_channel,:,i_trial)*correlator*keeptrial/n_trials;
%                     CI2fftrand(i_channel,:,i_rand) = (fft(CI2rand(i_channel,:,i_rand),Nfft));%CIrand(i_channel,:,i_rand)+Enorm(i_channel,:,i_trial)*correlator*keeptrial/n_trials;
%                 end
%             end
%             if n_boot > 0
%                 for i_boot = 1:n_boot
%                     CIfftboot(i_channel,:,i_boot) = (fft(CIboot(i_channel,:,i_boot),Nfft));%CIrand(i_channel,:,i_rand)+Enorm(i_channel,:,i_trial)*correlator*keeptrial/n_trials;
%                     CI1fftboot(i_channel,:,i_boot) = (fft(CI1boot(i_channel,:,i_boot),Nfft));%CIrand(i_channel,:,i_rand)+Enorm(i_channel,:,i_trial)*correlator*keeptrial/n_trials;
%                     CI2fftboot(i_channel,:,i_boot) = (fft(CI2boot(i_channel,:,i_boot),Nfft));%CIrand(i_channel,:,i_rand)+Enorm(i_channel,:,i_trial)*correlator*keeptrial/n_trials;
%                 end
%             end
%         end
% 
%         CIfft = CIfft(:,fidx2plot,:);
%         CI1fft = CI1fft(:,fidx2plot,:);
%         CI2fft = CI2fft(:,fidx2plot,:);
%         ideal_templatecfft = ideal_templatecfft(:,fidx2plot,:);
%         CIfftboot = CIfftboot(:,fidx2plot,:);
%         CI1fftboot = CI1fftboot(:,fidx2plot,:);
%         CI2fftboot = CI2fftboot(:,fidx2plot,:);
%         CIfftrand = CIfftrand(:,fidx2plot,:);
%         CI1fftrand = CI1fftrand(:,fidx2plot,:);
%         CI2fftrand = CI2fftrand(:,fidx2plot,:);
%         fE = fE(fidx2plot);
% 
%         CIfftrand_ci  = [prctile(CIfftrand ,CI_up,dim); prctile(CIfftrand ,CI_lo,dim)];
%         CI1fftrand_ci = [prctile(CI1fftrand,CI_up,dim); prctile(CI1fftrand,CI_lo,dim)];
%         CI2fftrand_ci = [prctile(CI2fftrand,CI_up,dim); prctile(CI2fftrand,CI_lo,dim)];
% 
%         if bSave
%             save([dir_out 'CIf.mat'],'CIfft','CI1fft','CI2fft','CIfftrand_ci','CI1fftrand_ci','CI2fftrand_ci','fE','ideal_templatecfft')
%         end
%     end
%         if bComplex_FFT_plot
%         error('Not validated yet...')
%         %% plot complex fft of revcorr
%         figure('Name', 'General kernel, Fourier amplitude domain'); hold on
%         plot_channels(fE, abs(CIfft)', 1:1:size(noise_E,1), @(x,y)(plot(x,y,'b-')));%, @plot, 1, size(undersmplE,1)
%         plot_channels(fE, max(abs(CIfft(:)))*abs(ideal_templatecfft)'/max(max(abs(ideal_templatecfft))),1:1:size(noise_E,1), @(x,y)(plot(x,y,'r-')));%, @plot, 1, size(undersmplE,1)
%         plot_channels(fE, prctile(abs(CIfftrand),CI_lo(2),3)', 1:1:size(noise_E,1), @(x,y)(plot(x,y,'k:')));%, @plot, 1, size(undersmplE,1)
%         plot_channels(fE, prctile(abs(CIfftboot),CI_lo(1),3)', 1:1:size(noise_E,1), @(x,y)(plot(x,y,'b:')));%, @plot, 1, size(undersmplE,1)
%         plot_channels(fE, prctile(abs(CIfftboot),CI_lo(2),3)', 1:1:size(noise_E,1), @(x,y)(plot(x,y,'b:')));%, @plot, 1, size(undersmplE,1)
%         set(findobj(gcf,'type','axes'),'YLim',[-0.1 1]*1.5*max(abs([CIfft(:); CI1fft(:); CI2fft(:)])));
%         saveas(gcf,[dir_out 'CI_afft_all.png'])
% 
%         figure('Name', 'Target-absent kernel, Fourier amplitude domain'); hold on
%         plot_channels(fE, abs(CI1fft)', 1:1:size(noise_E,1), @(x,y)(plot(x,y,'b-')));%, @plot, 1, size(undersmplE,1)
%         plot_channels(fE, max(abs(CI1fft(:)))*abs(ideal_templatecfft)'/max(max(abs(ideal_templatecfft))),1:1:size(noise_E,1), @(x,y)(plot(x,y,'r-')));%, @plot, 1, size(undersmplE,1)
%         plot_channels(fE, prctile(abs(CI1fftrand),CI_lo(2),3)', 1:1:size(noise_E,1), @(x,y)(plot(x,y,'k:')));%, @plot, 1, size(undersmplE,1)
%         plot_channels(fE, prctile(abs(CI1fftboot),CI_lo(1),3)', 1:1:size(noise_E,1), @(x,y)(plot(x,y,'b:')));%, @plot, 1, size(undersmplE,1)
%         plot_channels(fE, prctile(abs(CI1fftboot),CI_lo(2),3)', 1:1:size(noise_E,1), @(x,y)(plot(x,y,'b:')));%, @plot, 1, size(undersmplE,1)
%         set(findobj(gcf,'type','axes'),'YLim',[-0.1 1]*1.5*max(abs([CIfft(:); CI1fft(:); CI2fft(:)])));
%         saveas(gcf,[dir_out 'CI_afft_ta.png'])
% 
%         figure('Name', 'Target-present kernel, Fourier amplitude domain'); hold on
%         plot_channels(fE, abs(CI2fft)', 1:1:size(noise_E,1), @(x,y)(plot(x,y,'b-')));%, @plot, 1, size(undersmplE,1)
%         plot_channels(fE, max(abs(CI2fft(:)))*abs(ideal_templatecfft)'/max(max(abs(ideal_templatecfft))),1:1:size(noise_E,1), @(x,y)(plot(x,y,'r-')));%, @plot, 1, size(undersmplE,1)
%         plot_channels(fE, prctile(abs(CI2fftrand),CI_lo(2),3)', 1:1:size(noise_E,1), @(x,y)(plot(x,y,'k:')));%, @plot, 1, size(undersmplE,1)
%         plot_channels(fE, prctile(abs(CI2fftboot),CI_lo(1),3)', 1:1:size(noise_E,1), @(x,y)(plot(x,y,'b:')));%, @plot, 1, size(undersmplE,1)
%         plot_channels(fE, prctile(abs(CI2fftboot),CI_lo(2),3)', 1:1:size(noise_E,1), @(x,y)(plot(x,y,'b:')));%, @plot, 1, size(undersmplE,1)
%         set(findobj(gcf,'type','axes'),'YLim',[-0.1 1]*1.5*max(abs([CIfft(:); CI1fft(:); CI2fft(:)])));
%         saveas(gcf,[dir_out 'CI_afft_tp.png'])
% 
%         figure('Name', 'General kernel, Fourier phase domain'); hold on
%         plot_channels(fE, angle(CIfft)', 1:1:size(noise_E,1), @(x,y)(plot(x,y,'b-')));%, @plot, 1, size(undersmplE,1)
%         %plot_channels(fE, max(abs(CIfft(:)))*abs(ideal_templatecfft)'/max(max(abs(ideal_templatecfft))),1:1:size(noise_E,1), @(x,y)(plot(x,y,'r-')));%, @plot, 1, size(undersmplE,1)
%         plot_channels(fE, prctile(angle(CIfftrand),CI_lo(1),3)', 1:1:size(noise_E,1), @(x,y)(plot(x,y,'k:')));%, @plot, 1, size(undersmplE,1)
%         plot_channels(fE, prctile(angle(CIfftrand),CI_lo(2),3)', 1:1:size(noise_E,1), @(x,y)(plot(x,y,'k:')));%, @plot, 1, size(undersmplE,1)
%         plot_channels(fE, prctile(angle(CIfftboot),CI_lo(1),3)', 1:1:size(noise_E,1), @(x,y)(plot(x,y,'b:')));%, @plot, 1, size(undersmplE,1)
%         plot_channels(fE, prctile(angle(CIfftboot),CI_lo(2),3)', 1:1:size(noise_E,1), @(x,y)(plot(x,y,'b:')));%, @plot, 1, size(undersmplE,1)
%         set(findobj(gcf,'type','axes'),'YLim',[-3.5 3.5]);
%         saveas(gcf,[dir_out 'CI_pfft_all.png'])
% 
%         figure('Name', 'Target-absent kernel, Fourier phase domain'); hold on
%         plot_channels(fE, angle(CI1fft)', 1:1:size(noise_E,1), @(x,y)(plot(x,y,'b-')));%, @plot, 1, size(undersmplE,1)
%         %plot_channels(fE, max(abs(CIfft(:)))*abs(ideal_templatecfft)'/max(max(abs(ideal_templatecfft))),1:1:size(noise_E,1), @(x,y)(plot(x,y,'r-')));%, @plot, 1, size(undersmplE,1)
%         plot_channels(fE, prctile(angle(CI1fftrand),CI_lo(1),3)', 1:1:size(noise_E,1), @(x,y)(plot(x,y,'k:')));%, @plot, 1, size(undersmplE,1)
%         plot_channels(fE, prctile(angle(CI1fftrand),CI_lo(2),3)', 1:1:size(noise_E,1), @(x,y)(plot(x,y,'k:')));%, @plot, 1, size(undersmplE,1)
%         plot_channels(fE, prctile(angle(CI1fftboot),CI_lo(1),3)', 1:1:size(noise_E,1), @(x,y)(plot(x,y,'b:')));%, @plot, 1, size(undersmplE,1)
%         plot_channels(fE, prctile(angle(CI1fftboot),CI_lo(2),3)', 1:1:size(noise_E,1), @(x,y)(plot(x,y,'b:')));%, @plot, 1, size(undersmplE,1)
%         set(findobj(gcf,'type','axes'),'YLim',[-3.5 3.5]);
%         saveas(gcf,[dir_out 'CI1_pfft_all.png'])
% 
%         figure('Name', 'Target-present kernel, Fourier phase domain'); hold on
%         plot_channels(fE, angle(CI2fft)', 1:1:size(noise_E,1), @(x,y)(plot(x,y,'b-')));%, @plot, 1, size(undersmplE,1)
%         %plot_channels(fE, max(abs(CIfft(:)))*abs(ideal_templatecfft)'/max(max(abs(ideal_templatecfft))),1:1:size(noise_E,1), @(x,y)(plot(x,y,'r-')));%, @plot, 1, size(undersmplE,1)
%         plot_channels(fE, prctile(angle(CI2fftrand),CI_lo(1),3)', 1:1:size(noise_E,1), @(x,y)(plot(x,y,'k:')));%, @plot, 1, size(undersmplE,1)
%         plot_channels(fE, prctile(angle(CI2fftrand),CI_lo(2),3)', 1:1:size(noise_E,1), @(x,y)(plot(x,y,'k:')));%, @plot, 1, size(undersmplE,1)
%         plot_channels(fE, prctile(angle(CI2fftboot),CI_lo(1),3)', 1:1:size(noise_E,1), @(x,y)(plot(x,y,'b:')));%, @plot, 1, size(undersmplE,1)
%         plot_channels(fE, prctile(angle(CI2fftboot),CI_lo(2),3)', 1:1:size(noise_E,1), @(x,y)(plot(x,y,'b:')));%, @plot, 1, size(undersmplE,1)
%         set(findobj(gcf,'type','axes'),'YLim',[-3.5 3.5]);
%         saveas(gcf,[dir_out 'CI2_pfft_all.png'])
%     end
%     
%     if bExtract_metrics_using_Fourier
%         error('Not validated yet...')
%         %% Extract metrics through fourier
%         clear maxA_CIfftboot idmaxA_CIfftboot freqmaxA_CIfftboot phasemaxA_CIfftboot Aftarget_CIfftboot phaseftarget_CIfftboot
%         clear maxA_CI1fftboot idmaxA_CI1fftboot freqmaxA_CI1fftboot phasemaxA_CI1fftboot Aftarget_CI1fftboot phaseftarget_CI1fftboot
%         clear maxA_CI2fftboot idmaxA_CI2fftboot freqmaxA_CI2fftboot phasemaxA_CI2fftboot Aftarget_CI2fftboot phaseftarget_CI2fftboot
%         clear maxA_CIfftrand idmaxA_CIfftrand freqmaxA_CIfftrand phasemaxA_CIfftrand Aftarget_CIfftrand phaseftarget_CIfftrand
%         clear maxA_CI1fftrand idmaxA_CI1fftrand freqmaxA_CI1fftrand phasemaxA_CI1fftrand Aftarget_CI1fftrand phaseftarget_CI1fftrand
%         clear maxA_CI2fftrand idmaxA_CI2fftrand freqmaxA_CI2fftrand phasemaxA_CI2fftrand Aftarget_CI2fftrand phaseftarget_CI2fftrand
% 
%         [~, idxtarget] = min(abs(fE-cfg_game.fm));
% 
%         for i_channel = 1:size(noise_E,1)
%             % amplitude and phase at target, max amplitude, frequency and phase at max amplitude
%             if n_boot > 0
%                 for i_boot = 1:n_boot
%                     [maxA_CIfftboot(i_channel,i_boot), idmaxA_CIfftboot(i_channel,i_boot)] = max(abs(CIfftboot(i_channel,:,i_boot)));
%                     [maxA_CI1fftboot(i_channel,i_boot), idmaxA_CI1fftboot(i_channel,i_boot)] = max(abs(CI1fftboot(i_channel,:,i_boot)));
%                     [maxA_CI2fftboot(i_channel,i_boot), idmaxA_CI2fftboot(i_channel,i_boot)] = max(abs(CI2fftboot(i_channel,:,i_boot)));
% 
%                     freqmaxA_CIfftboot(i_channel,i_boot) = fE(idmaxA_CIfftboot(i_channel,i_boot));
%                     phasemaxA_CIfftboot(i_channel,i_boot) = angle(CIfftboot(i_channel,idmaxA_CIfftboot(i_channel,i_boot),i_boot));
%                     freqmaxA_CI1fftboot(i_channel,i_boot) = fE(idmaxA_CI1fftboot(i_channel,i_boot));
%                     phasemaxA_CI1fftboot(i_channel,i_boot) = angle(CI1fftboot(i_channel,idmaxA_CI1fftboot(i_channel,i_boot),i_boot));
%                     freqmaxA_CI2fftboot(i_channel,i_boot) = fE(idmaxA_CI2fftboot(i_channel,i_boot));
%                     phasemaxA_CI2fftboot(i_channel,i_boot) = angle(CI2fftboot(i_channel,idmaxA_CI2fftboot(i_channel,i_boot),i_boot));
% 
%                     Aftarget_CIfftboot(i_channel,i_boot) = abs(CIfftboot(i_channel,idxtarget,i_boot));
%                     phaseftarget_CIfftboot(i_channel,i_boot) = angle(CIfftboot(i_channel,idxtarget,i_boot));
%                     Aftarget_CI1fftboot(i_channel,i_boot) = abs(CI1fftboot(i_channel,idxtarget,i_boot));
%                     phaseftarget_CI1fftboot(i_channel,i_boot) = angle(CI1fftboot(i_channel,idxtarget,i_boot));
%                     Aftarget_CI2fftboot(i_channel,i_boot) = abs(CI2fftboot(i_channel,idxtarget,i_boot));
%                     phaseftarget_CI2fftboot(i_channel,i_boot) = angle(CI2fftboot(i_channel,idxtarget,i_boot));
%                 end
%             end
%             if n_rand > 0
%                 for i_rand = 1:n_rand
%                     [maxA_CIfftrand(i_channel,i_rand), idmaxA_CIfftrand(i_channel,i_rand)] = max(abs(CIfftrand(i_channel,:,i_rand)));
%                     [maxA_CI1fftrand(i_channel,i_rand), idmaxA_CI1fftrand(i_channel,i_rand)] = max(abs(CI1fftrand(i_channel,:,i_rand)));
%                     [maxA_CI2fftrand(i_channel,i_rand), idmaxA_CI2fftrand(i_channel,i_rand)] = max(abs(CI2fftrand(i_channel,:,i_rand)));
% 
%                     freqmaxA_CIfftrand(i_channel,i_rand) = fE(idmaxA_CIfftrand(i_channel,i_rand));
%                     phasemaxA_CIfftrand(i_channel,i_rand) = angle(CIfftrand(i_channel,idmaxA_CIfftrand(i_channel,i_rand),i_rand));
%                     freqmaxA_CI1fftrand(i_channel,i_rand) = fE(idmaxA_CI1fftrand(i_channel,i_rand));
%                     phasemaxA_CI1fftrand(i_channel,i_rand) = angle(CI1fftrand(i_channel,idmaxA_CI1fftrand(i_channel,i_rand),i_rand));
%                     freqmaxA_CI2fftrand(i_channel,i_rand) = fE(idmaxA_CI2fftrand(i_channel,i_rand));
%                     phasemaxA_CI2fftrand(i_channel,i_rand) = angle(CI2fftrand(i_channel,idmaxA_CI2fftrand(i_channel,i_rand),i_rand));
% 
%                     Aftarget_CIfftrand(i_channel,i_rand) = abs(CIfftrand(i_channel,idxtarget,i_rand));
%                     phaseftarget_CIfftrand(i_channel,i_rand) = angle(CIfftrand(i_channel,idxtarget,i_rand));
%                     Aftarget_CI1fftrand(i_channel,i_rand) = abs(CI1fftrand(i_channel,idxtarget,i_rand));
%                     phaseftarget_CI1fftrand(i_channel,i_rand) = angle(CI1fftrand(i_channel,idxtarget,i_rand));
%                     Aftarget_CI2fftrand(i_channel,i_rand) = abs(CI2fftrand(i_channel,idxtarget,i_rand));
%                     phaseftarget_CI2fftrand(i_channel,i_rand) = angle(CI2fftrand(i_channel,idxtarget,i_rand));
%                 end
%             end
%         end
% 
%         phaseftarget_CIfftboot = unwrap(phaseftarget_CIfftboot);
%         phaseftarget_CI1fftboot = unwrap(phaseftarget_CI1fftboot);
%         phaseftarget_CI2fftboot = unwrap(phaseftarget_CI2fftboot);
%         phaseftarget_CIfftrand = unwrap(phaseftarget_CIfftrand);
%         phaseftarget_CI1fftrand = unwrap(phaseftarget_CI1fftrand);
%         phaseftarget_CI2fftrand = unwrap(phaseftarget_CI2fftrand);
% 
%         phasemaxA_CIfftboot = unwrap(phasemaxA_CIfftboot);
%         phasemaxA_CI1fftboot = unwrap(phasemaxA_CI1fftboot);
%         phasemaxA_CI2fftboot = unwrap(phasemaxA_CI2fftboot);
%         phasemaxA_CIfftrand = unwrap(phasemaxA_CIfftrand);
%         phasemaxA_CI1fftrand = unwrap(phasemaxA_CI1fftrand);
%         phasemaxA_CI2fftrand = unwrap(phasemaxA_CI2fftrand);
% 
% 
%         CIfftrand_ci = [prctile(CIfftrand,CI_lo(2),3); prctile(CIfftrand,CI_lo(1),3)];
%         CI1fftrand_ci = [prctile(CI1fftrand,CI_lo(2),3); prctile(CI1fftrand,CI_lo(1),3)];
%         CI2fftrand_ci = [prctile(CI2fftrand,CI_lo(2),3); prctile(CI2fftrand,CI_lo(1),3)];
% 
%         CIp = prctile(phaseftarget_CIfftboot(:),50);
%         CI1p = prctile(phaseftarget_CI1fftboot(:),50);
%         CI2p = prctile(phaseftarget_CI2fftboot(:),50);
% 
%         save('CIp','CIp','CI1p','CI2p')
% 
% 
%         figure('Name', 'Phase at target frequency'); 
%         polarplot(prctile(phaseftarget_CIfftboot(:),50),1,'bs',...
%             prctile(phaseftarget_CI1fftboot(:),50),1.1,'k*--',...
%             prctile(phaseftarget_CI2fftboot(:),50),1.2,'ko',...
%             [pi pi],[0 0.5],'r-',...
%             linspace(prctile(phaseftarget_CIfftboot(:),CI_lo(1)),prctile(phaseftarget_CIfftboot(:),CI_lo(2)),100),ones(1,100),'b-',...
%             linspace(prctile(phaseftarget_CI1fftboot(:),CI_lo(1)),prctile(phaseftarget_CI1fftboot(:),CI_lo(2)),100),1.1*ones(1,100),'k--',...
%             linspace(prctile(phaseftarget_CI2fftboot(:),CI_lo(1)),prctile(phaseftarget_CI2fftboot(:),CI_lo(2)),100),1.2*ones(1,100),'k-')
%         legend({'general','target-absent','target-present'},'Location','northeast')
%         set(gca, 'ThetaAxisUnits', 'radians');rticks([]);rlim([0 1.5])
%         saveas(gcf,'CIphase4Hz.png')
% 
%         figure('Name', 'Phase at peak frequency'); 
%         polarplot(prctile(phasemaxA_CIfftboot(:),50),1,'bs',...
%             prctile(phasemaxA_CI1fftboot(:),50),1.1,'k*--',...
%             prctile(phasemaxA_CI2fftboot(:),50),1.2,'ko',...
%             [pi pi],[0 0.5],'r-',...
%             linspace(prctile(phasemaxA_CIfftboot(:),CI_lo(1)),prctile(phasemaxA_CIfftboot(:),CI_lo(2)),100),ones(1,100),'b-',...
%             linspace(prctile(phasemaxA_CI1fftboot(:),CI_lo(1)),prctile(phasemaxA_CI1fftboot(:),CI_lo(2)),100),1.1*ones(1,100),'k--',...
%             linspace(prctile(phasemaxA_CI2fftboot(:),CI_lo(1)),prctile(phasemaxA_CI2fftboot(:),CI_lo(2)),100),1.2*ones(1,100),'k-')
%         legend({'general','target-absent','target-present'},'Location','northeast')
%         set(gca, 'ThetaAxisUnits', 'radians');rticks([]);rlim([0 1.5])
%         saveas(gcf,'CIphasemaxA.png')
% 
%         figure('Name', 'Amplitude at target frequency'); hold on
%         plot(1,prctile(Aftarget_CIfftboot(:),50),'bs');
%         plot([1,1],[prctile(Aftarget_CIfftboot(:),CI_lo(1)),prctile(Aftarget_CIfftboot(:),CI_lo(2))],'b');
%         plot([1.1,1.1],[prctile(Aftarget_CIfftrand(:),CI_lo(1)),prctile(Aftarget_CIfftrand(:),CI_lo(2))],'b:');
%         plot(2,prctile(Aftarget_CI1fftboot(:),50),'k*');
%         plot([2,2],[prctile(Aftarget_CI1fftboot(:),CI_lo(1)),prctile(Aftarget_CI1fftboot(:),CI_lo(2))],'k');
%         plot([2.1,2.1],[prctile(Aftarget_CI1fftrand(:),CI_lo(1)),prctile(Aftarget_CI1fftrand(:),CI_lo(2))],'k:');
%         plot(3,prctile(Aftarget_CI2fftboot(:),50),'ko');
%         plot([3,3],[prctile(Aftarget_CI2fftboot(:),CI_lo(1)),prctile(Aftarget_CI2fftboot(:),CI_lo(2))],'k');
%         plot([3.1,3.1],[prctile(Aftarget_CI2fftrand(:),CI_lo(1)),prctile(Aftarget_CI2fftrand(:),CI_lo(2))],'k:');
%         title('Amplitude at target frequency'); xlim([0.5 3.5])
%         set(gca, 'XTick', [1,2,3], 'XTickLabels', {'all', 'target-absent', 'target-present'})
%         ylabel('Amplitude at target frequency')
%         saveas(gcf,'CIAftarget.png')
% 
%         figure('Name', 'peak frequency'); hold on
%         plot(1,prctile(freqmaxA_CIfftboot(:),50),'bs');
%         plot([1,1],[prctile(freqmaxA_CIfftboot(:),CI_lo(1)),prctile(freqmaxA_CIfftboot(:),CI_lo(2))],'b');
%         plot([1.1,1.1],[prctile(freqmaxA_CIfftrand(:),CI_lo(1)),prctile(freqmaxA_CIfftrand(:),CI_lo(2))],'b:');
%         plot(2,prctile(freqmaxA_CI1fftboot(:),50),'k*');
%         plot([2,2],[prctile(freqmaxA_CI1fftboot(:),CI_lo(1)),prctile(freqmaxA_CI1fftboot(:),CI_lo(2))],'k');
%         plot([2.1,2.1],[prctile(freqmaxA_CI1fftrand(:),CI_lo(1)),prctile(freqmaxA_CI1fftrand(:),CI_lo(2))],'k:');
%         plot(3,prctile(freqmaxA_CI2fftboot(:),50),'ko');
%         plot([3,3],[prctile(freqmaxA_CI2fftboot(:),CI_lo(1)),prctile(freqmaxA_CI2fftboot(:),CI_lo(2))],'k');
%         plot([3.1,3.1],[prctile(freqmaxA_CI2fftrand(:),CI_lo(1)),prctile(freqmaxA_CI2fftrand(:),CI_lo(2))],'k:');
%         plot([0 4],[4 4],'r--')
%         title('peak frequency'); xlim([0.5 3.5])
%         set(gca, 'XTick', [1,2,3], 'XTickLabels', {'all', 'target-absent', 'target-present'})
%         ylabel('peak frequency (Hz)')
%         saveas(gcf,'CIfpeak.png')
%     end
% 
%     if bExtract_metrics_using_xcorr
%         %% Extract metrics through xcorr
% 
%         maxlag = cfg_game.fs/cfg_game.fm/undersampling/2;
% 
%         lags = (-maxlag:maxlag)/(cfg_game.fs/undersampling);%[-tE(end:-1:2) tE];
%         [xCI] = xcorr(CI,ideal_template,maxlag,'coeff'); [corrCI,phaseCI] = max(xCI); phaseCI=lags(phaseCI);
%         [xCI1] = xcorr(CI1,ideal_template,maxlag,'coeff'); [corrCI1,phaseCI1] = max(xCI1); phaseCI1=lags(phaseCI1);
%         [xCI2] = xcorr(CI2,ideal_template,maxlag,'coeff'); [corrCI2,phaseCI2] = max(xCI2); phaseCI2=lags(phaseCI2);
% 
%         % figure; hold on
%         % plot(lags, xCI);
%         % plot(lags, xCI1);
%         % plot(lags, xCI2);
%         % legend({'all', 'target-absent', 'target-present'})
% 
%         for i_rand = 1:n_rand
%             [xCIrand] = xcorr(CIrand(1,:,i_rand),ideal_template,maxlag,'coeff'); [corrCIrand(i_rand),phaseCIrand(i_rand)] = max(xCIrand); phaseCIrand(i_rand)=lags(phaseCIrand(i_rand));
%             [xCI1rand] = xcorr(CI1rand(1,:,i_rand),ideal_template,maxlag,'coeff'); [corrCI1rand(i_rand),phaseCI1rand(i_rand)] = max(xCI1rand); phaseCI1rand(i_rand)=lags(phaseCI1rand(i_rand));
%             [xCI2rand] = xcorr(CI2rand(1,:,i_rand),ideal_template,maxlag,'coeff'); [corrCI2rand(i_rand),phaseCI2rand(i_rand)] = max(xCI2rand); phaseCI2rand(i_rand)=lags(phaseCI2rand(i_rand));
%         end
% 
%         for i_boot = 1:n_boot
%             [xCIboot] = xcorr(CIboot(1,:,i_boot),ideal_template,maxlag,'coeff'); [corrCIboot(i_boot),phaseCIboot(i_boot)] = max(xCIboot); phaseCIboot(i_boot)=lags(phaseCIboot(i_boot));
%             [xCI1boot] = xcorr(CI1boot(1,:,i_boot),ideal_template,maxlag,'coeff'); [corrCI1boot(i_boot),phaseCI1boot(i_boot)] = max(xCI1boot); phaseCI1boot(i_boot)=lags(phaseCI1boot(i_boot));
%             [xCI2boot] = xcorr(CI2boot(1,:,i_boot),ideal_template,maxlag,'coeff'); [corrCI2boot(i_boot),phaseCI2boot(i_boot)] = max(xCI2boot); phaseCI2boot(i_boot)=lags(phaseCI2boot(i_boot));
%         end
% 
%         figure; hold on
%         plot(1,corrCI,'bo');plot([1,1],[prctile(corrCIboot(:),CI_lo(1)),prctile(corrCIboot(:),CI_lo(2))],'b');plot([1.1,1.1],[prctile(corrCIrand(:),CI_lo(1)),prctile(corrCIrand(:),CI_lo(2))],'k')
%         plot(2,corrCI1,'bo');plot([2,2],[prctile(corrCI1boot(:),CI_lo(1)),prctile(corrCI1boot(:),CI_lo(2))],'b');plot([2.1,2.1],[prctile(corrCI1rand(:),CI_lo(1)),prctile(corrCI1rand(:),CI_lo(2))],'k')
%         plot(3,corrCI2,'bo');plot([3,3],[prctile(corrCI2boot(:),CI_lo(1)),prctile(corrCI2boot(:),CI_lo(2))],'b');plot([3.1,3.1],[prctile(corrCI2rand(:),CI_lo(1)),prctile(corrCI2rand(:),CI_lo(2))],'k')
%         title('correlation with ideal template'); ylim([0 1]);xlim([0.5 3.5])
%         set(gca, 'XTick', [1,2,3], 'XTickLabels', {'all', 'target-absent', 'target-present'})
%         ylabel('correlation coefficient (r)')
%         saveas(gcf,'CIcorr.png')
% 
%         figure; hold on
%         plot(1,phaseCI,'bo');plot([1,1],[prctile(phaseCIboot(:),CI_lo(1)),prctile(phaseCIboot(:),CI_lo(2))],'b');plot([1.1,1.1],[prctile(phaseCIrand(:),CI_lo(1)),prctile(phaseCIrand(:),CI_lo(2))],'k')
%         plot(2,phaseCI1,'bo');plot([2,2],[prctile(phaseCI1boot(:),CI_lo(1)),prctile(phaseCI1boot(:),CI_lo(2))],'b');plot([2.1,2.1],[prctile(phaseCI1rand(:),CI_lo(1)),prctile(phaseCI1rand(:),CI_lo(2))],'k')
%         plot(3,phaseCI2,'bo');plot([3,3],[prctile(phaseCI2boot(:),CI_lo(1)),prctile(phaseCI2boot(:),CI_lo(2))],'b');plot([3.1,3.1],[prctile(phaseCI2rand(:),CI_lo(1)),prctile(phaseCI2rand(:),CI_lo(2))],'k')
%         title('phase shift with ideal template'); ylim([-1/cfg_game.fm/2-0.01 1/cfg_game.fm/2+0.01]);xlim([0.5 3.5])
%         set(gca, 'XTick', [1,2,3], 'XTickLabels', {'all', 'target-absent', 'target-present'})
%         ylabel('phase shift (s)')
%         saveas(gcf,'CIphaseshift.png')
%     end
%     
%     if bComplex_FFT_revcorr
%         %% Complex fft revcorr
%         error('Not validated yet...')
%         clear noise_E_fft noise_E_cfft ideal_templatefft ideal_templatecfft
% 
%         n_trials = length(n_responses);
%         Nfft = 512*100;
%         xl = [0 fcut_noiseE];
%         fE = (1:Nfft)/Nfft*(cfg.fs/cfg.save_undersmpl);
%         fidx2plot = fE>=xl(1) & fE<=xl(2);
% 
%         for i_channel = 1:size(noise_E,1)
%             for i_trial = 1:n_trials
%                 noise_E_wm = noise_E(i_channel,:,i_trial);% - mean(noise_E(i_channel,:,i_trial));
%                 noise_E_cfft(i_channel,:,i_trial) = fft(noise_E_wm,Nfft,2);
%             end
%             ideal_templatecfft(i_channel,:) = (fft(ideal_template(i_channel,:),Nfft));
% 
%             % for i_trial = 1:2*Ntemplate
%             %     noisysignal_E_wm = noisysignal_E(i_channel,:,i_trial);% - mean(noisysignal_E(i_channel,:,i_trial));
%             %     noisysignal_E_fft(i_channel,:,i_trial) = abs(fft(noisysignal_E_wm, Nfft, 2)); 
%             % end
%             % ideal_noisytemplatefft = mean(noisysignal_E_fft(:,:,Ntemplate+1:end),3)-mean(noisysignal_E_fft(:,:,1:Ntemplate),3);
%         end
%         noise_E_cfft = noise_E_cfft(:,fidx2plot,:);
%         ideal_templatecfft = ideal_templatecfft(:,fidx2plot,:);
%         %ideal_noisytemplatefft = ideal_noisytemplatefft(:,fidx2plot,:);
%         fE = fE(fidx2plot);
% 
%         %subselection of trials
%         trials2analyse = 1:n_trials;%m(i_trial)<=median(m);%(iscorrect(i_trial)==0);%responsetime(i_trial)>=median(responsetime);%
% 
%         [CI_F, CIrand_F, CIboot_F, ~, CI2_F, CI2rand_F, CI2boot_F, CI1_F, CI1rand_F, CI1boot_F] = computeCI(n_signal(trials2analyse), n_responses(trials2analyse), noise_E_cfft(:,:,trials2analyse), n_rand, n_boot, 'no');
% 
%         %% phase shift from complex fft
% 
%         for i_boot=1:n_boot
%             phaseshift_CIboot(i_boot) = (angle(CIboot_F(1,427,i_boot)));
%             phaseshift_CI1boot(i_boot) = (angle(CI1boot_F(1,427,i_boot)));
%             phaseshift_CI2boot(i_boot) = (angle(CI2boot_F(1,427,i_boot)));
%         end
%         for i_rand=1:n_rand
%             phaseshift_CIrand(i_rand) = (angle(CIrand_F(1,427,i_rand)));
%             phaseshift_CI1rand(i_rand) = (angle(CI1rand_F(1,427,i_rand)));
%             phaseshift_CI2rand(i_rand) = (angle(CI2rand_F(1,427,i_rand)));
%         end
% 
%         % phaseshift_CIboot = unwrap(phaseshift_CIboot-pi)*(1/4)/(2*pi);
%         % phaseshift_CI1boot = unwrap(phaseshift_CI1boot-pi)*(1/4)/(2*pi);
%         % phaseshift_CI2boot = unwrap(phaseshift_CI2boot-pi)*(1/4)/(2*pi);
%         % phaseshift_CIrand = (phaseshift_CIrand-pi)*(1/4)/(2*pi);
%         % phaseshift_CI1rand = (phaseshift_CI1rand-pi)*(1/4)/(2*pi);
%         % phaseshift_CI2rand = (phaseshift_CI2rand-pi)*(1/4)/(2*pi);
% 
%         figure;plot([1,2,3],[phaseshift_CIboot',phaseshift_CI1boot',phaseshift_CI2boot'],'o')
%         figure;plot([1,2,3],[phaseshift_CIrand',phaseshift_CI1rand',phaseshift_CI2rand'],'o')
% 
%         %% revcorr in the Fourier domain (phase-discarded)
% 
%         clear noise_E_fft noise_E_cfft ideal_templatefft ideal_templatecfft
% 
%         n_trials = length(n_responses);
%         Nfft = 512*100;
%         xl = [0 fcut_noiseE];
%         fE = (1:Nfft)/Nfft*(cfg.fs/cfg.save_undersmpl);
%         fidx2plot = fE>=xl(1) & fE<=xl(2);
% 
%         for i_channel = 1:size(noise_E,1)
%             for i_trial = 1:n_trials
%                 noise_E_wm = noise_E(i_channel,:,i_trial);% - mean(noise_E(i_channel,:,i_trial));
%                 noise_E_fft(i_channel,:,i_trial) = abs(fft(noise_E_wm, Nfft, 2)); 
%                 %noise_E_cfft(i_channel,:,i_trial) = fft(noise_E_wm,Nfft,2);
%             end
%             ideal_templatefft(i_channel,:) = abs(fft(ideal_template(i_channel,:),Nfft)); 
%             %ideal_templatecfft(i_channel,:) = (fft(ideal_template(i_channel,idx2fft),Nfft));
% 
%         %     for i_trial = 1:2*Ntemplate
%         %         noisysignal_E_wm = noisysignal_E(i_channel,:,i_trial);% - mean(noisysignal_E(i_channel,:,i_trial));
%         %         noisysignal_E_fft(i_channel,:,i_trial) = abs(fft(noisysignal_E_wm, Nfft, 2)); 
%         %     end
%         %     ideal_noisytemplatefft = mean(noisysignal_E_fft(:,:,Ntemplate+1:end),3)-mean(noisysignal_E_fft(:,:,1:Ntemplate),3);
%         end
%         noise_E_fft = noise_E_fft(:,fidx2plot,:);
%         ideal_templatefft = ideal_templatefft(:,fidx2plot,:);
%         %ideal_noisytemplatefft = ideal_noisytemplatefft(:,fidx2plot,:);
%         fE = fE(fidx2plot);
% 
%         %subselection of trials
%         trials2analyse = 1:n_trials;%m(i_trial)<=median(m);%(iscorrect(i_trial)==0);%responsetime(i_trial)>=median(responsetime);%
% 
%         [CI_F, CIrand_F, CIboot_F, ~, CI2_F, CI2rand_F, CI2boot_F, CI1_F, CI1rand_F, CI1boot_F] = computeCI(n_signal(trials2analyse), n_responses(trials2analyse), noise_E_fft(:,:,trials2analyse), n_rand, n_boot, 'no');
% 
%         figure('Name', 'General kernel');
%         plot_channels(fE, max(abs(CI_F(:)))*abs(ideal_templatefft')/max(max(abs(ideal_templatefft))),1:1:size(noise_E,1), @(x,y)(plot(x,y,'r-')));%, @plot, 1, size(undersmplE,1)
%         %plot_channels(fE, 0.3*abs(ideal_noisytemplatefft')/max(max(abs(ideal_noisytemplatefft))),1:1:size(noise_E,1), @(x,y)(plot(x,y,'g-')));%, @plot, 1, size(undersmplE,1)
%         plot_channels(fE, abs(CI_F)', 1:1:size(noise_E,1), @(x,y)(plot(x,y,'b-')));%, @plot, 1, size(undersmplE,1)
%         plot_channels(fE, prctile(abs(CIrand_F),CI_lo(2),3)', 1:1:size(noise_E,1), @(x,y)(plot(x,y,'k:')));%, @plot, 1, size(undersmplE,1)
%         plot_channels(fE, prctile(abs(CIboot_F),CI_lo(1),3)', 1:1:size(noise_E,1), @(x,y)(plot(x,y,'b:')));%, @plot, 1, size(undersmplE,1)
%         plot_channels(fE, prctile(abs(CIboot_F),CI_lo(2),3)', 1:1:size(noise_E,1), @(x,y)(plot(x,y,'b:')));%, @plot, 1, size(undersmplE,1)
%         set(findobj(gcf,'type','axes'),'YLim',[-0.1 1]*1.5*max(abs([CI_F(:); CI1_F(:); CI2_F(:)])));
%         saveas(gcf,'CIF_all.png')
% 
%         figure('Name', 'Target-absent');
%         plot_channels(fE, abs(CI1_F)', 1:1:size(noise_E,1), @(x,y)(plot(x,y,'b-')));%, @plot, 1, size(undersmplE,1)
%         plot_channels(fE, max(abs(CI1_F(:)))*abs(ideal_templatefft')/max(max(abs(ideal_templatefft))),1:1:size(noise_E,1), @(x,y)(plot(x,y,'r-')));%, @plot, 1, size(undersmplE,1)
%         %plot_channels(fE, 0.3*abs(ideal_noisytemplatefft')/max(max(abs(ideal_noisytemplatefft))),1:1:size(noise_E,1), @(x,y)(plot(x,y,'g-')));%, @plot, 1, size(undersmplE,1)
%         plot_channels(fE, prctile(abs(CI1rand_F),CI_lo(2),3)', 1:1:size(noise_E,1), @(x,y)(plot(x,y,'k:')));%, @plot, 1, size(undersmplE,1)
%         plot_channels(fE, prctile(abs(CI1boot_F),CI_lo(1),3)', 1:1:size(noise_E,1), @(x,y)(plot(x,y,'b:')));%, @plot, 1, size(undersmplE,1)
%         plot_channels(fE, prctile(abs(CI1boot_F),CI_lo(2),3)', 1:1:size(noise_E,1), @(x,y)(plot(x,y,'b:')));%, @plot, 1, size(undersmplE,1)
%         set(findobj(gcf,'type','axes'),'YLim',[-0.1 1]*1.5*max(abs([CI_F(:); CI1_F(:); CI2_F(:)])));
%         saveas(gcf,'CIF_ta.png')
% 
%         figure('Name', 'Target-present');
%         plot_channels(fE, abs(CI2_F)', 1:1:size(noise_E,1), @(x,y)(plot(x,y,'b-')));%, @plot, 1, size(undersmplE,1)
%         plot_channels(fE, max(abs(CI2_F(:)))*abs(ideal_templatefft')/max(max(abs(ideal_templatefft))),1:1:size(noise_E,1), @(x,y)(plot(x,y,'r-')));%, @plot, 1, size(undersmplE,1)
%         %plot_channels(fE, 0.3*abs(ideal_noisytemplatefft')/max(max(abs(ideal_noisytemplatefft))),1:1:size(noise_E,1), @(x,y)(plot(x,y,'g-')));%, @plot, 1, size(undersmplE,1)
%         plot_channels(fE, prctile(abs(CI2rand_F),CI_lo(2),3)', 1:1:size(noise_E,1), @(x,y)(plot(x,y,'k:')));%, @plot, 1, size(undersmplE,1)
%         plot_channels(fE, prctile(abs(CI2boot_F),CI_lo(1),3)', 1:1:size(noise_E,1), @(x,y)(plot(x,y,'b:')));%, @plot, 1, size(undersmplE,1)
%         plot_channels(fE, prctile(abs(CI2boot_F),CI_lo(2),3)', 1:1:size(noise_E,1), @(x,y)(plot(x,y,'b:')));%, @plot, 1, size(undersmplE,1)
%         set(findobj(gcf,'type','axes'),'YLim',[-0.1 1]*1.5*max(abs([CI_F(:); CI1_F(:); CI2_F(:)])));
%         saveas(gcf,'CIF_tp.png')
% 
%         % figure;
%         % plot_channels(fE, (angle(CI_F')), 1:1:size(noise_E,1), @(x,y)(plot(x,y,'b-')));%, @plot, 1, size(undersmplE,1)
%         % plot_channels(fE, (angle(ideal_templatecfft')),1:1:size(noise_E,1), @(x,y)(plot(x,y,'r-')));%, @plot, 1, size(undersmplE,1)
%         % plot_channels(fE, prctile(unwrap(angle(CIrand_F')),CIlim(2),3), 1:1:size(noise_E,1), @(x,y)(plot(x,y,'k:')));%, @plot, 1, size(undersmplE,1)
%         % plot_channels(fE, prctile(unwrap(angle(CIrand_F')),CIlim(1),3), 1:1:size(noise_E,1), @(x,y)(plot(x,y,'k:')));%, @plot, 1, size(undersmplE,1)
%         % plot_channels(fE, prctile(unwrap(angle(CIboot_F')),CIlim(2),3), 1:1:size(noise_E,1), @(x,y)(plot(x,y,'b:')));%, @plot, 1, size(undersmplE,1)
%         % plot_channels(fE, prctile(unwrap(angle(CIboot_F')),CIlim(1),3), 1:1:size(noise_E,1), @(x,y)(plot(x,y,'b:')));%, @plot, 1, size(undersmplE,1)
% 
%         % xlim([0 20])
% 
%         % figure;
%         % plot_channels(fE, abs(CI_F)', 1:1:size(noise_E,1), @(x,y)(plot(x,y,'b-')))%, @plot, 1, size(undersmplE,1))
%         % plot_channels(fE, ideal_templatefft'/max(max(ideal_templatefft)),1:1:size(noise_E,1), @(x,y)(plot(x,y,'r-')));%, @plot, 1, size(undersmplE,1)
%         % % for i_rand = 1:n_rand
%         % %     plot_channels((1:length(CI))/(cfg.fs/cfg.save_undersmpl), CIrand(:,:,i_rand)',1:1:size(noise_E,1))%, @plot, 1, size(undersmplE,1))
%         % % end
%         % plot_channels(fE, prctile(abs(CIrand_F),CIlim(2),3)',1:1:size(noise_E,1), @(x,y)(plot(x,y,'b:')))%, @plot, 1, size(undersmplE,1))
%         % plot_channels(fE, prctile(abs(CIrand_F),CIlim(1),3)',1:1:size(noise_E,1), @(x,y)(plot(x,y,'b:')))%, @plot, 1, size(undersmplE,1))
%         % 
%         % %ylim([-1 1])
%         % xlim([0 20])
% 
%         % %% Analyze noise in the internal representation domain
%         % 
%         % clear noise_E_fft noise_E_cfft ideal_templatefft ideal_templatecfft
%         % 
%         % foldername = 'ListeBruit';
%         % [noise_E] = noise_EPSM_converter(foldername, ListStim, data_passation.n_stim, 'yes');
%         % 
%         % foldername = 'ListeSignal';
%         % ListSignals(2).name = 'target.wav';ListSignals(1).name = 'nontarget.wav';
%         % 
%         % [signal_E, step] = noise_EPSM_converter(foldername,ListSignals,[1,2], 'yes');
%         % 
%         % ideal_template = signal_E(:,:,2)-signal_E(:,:,1); ideal_template = ideal_template/max(max(max(ideal_template)));
%         % 
%         % %% Compute EPSM CI and target-present and target-absent sub-CI
%         % 
%         % n_rand = 200;
%         % n_boot = 200;
%         % n_trials = size(noise_E,3);
%         % n_mchan = size(noise_E,2);
%         % n_fchan = size(noise_E,1);
%         % 
%         % noise_E_wm = zeros(size(noise_E));
%         % for i_channel = 1:size(noise_E,1)
%         %     for i_trial = 1:n_trials
%         %         noise_E_wm(i_channel,:,i_trial) = noise_E(i_channel,:,i_trial) - mean(noise_E(i_channel,:,i_trial));
%         %     end
%         % end
%         % 
%         % %subselection of trials
%         % trials2analyze = 1:n_trials;%m(i_trial)<=median(m);%(iscorrect(i_trial)==0);%responsetime(i_trial)>=median(responsetime);%
%         % 
%         % [CI_EPSM, CIrand_EPSM, CIboot_EPSM, ResponseMatrix, CI2_EPSM, CI2rand_EPSM, CI2boot_EPSM, CI1_EPSM, CI1rand_EPSM, CI1boot_EPSM] = computeCI(n_signal(trials2analyze),n_response(trials2analyze), noise_E(:,:,trials2analyze), n_rand, n_boot, 'yes');
%         % 
%         % %% plot EPSM revcorr
%         % 
%         % figure('Name', 'General kernel'); hold on
%         % plot_channels(step.fmc, CI_EPSM', 1:n_fchan, @(x,y)(semilogx(x,y,'b-')),1,n_fchan)%, @plot, 1, size(undersmplE,1))
%         % plot_channels(step.fmc, prctile(CIboot_EPSM,CIlim(1),3)',1:n_fchan, @(x,y)(semilogx(x,y,'b:')),1,n_fchan)%, @plot, 1, size(undersmplE,1))
%         % %plot_channels(step.fmc, 60*mean(Template)',1:n_mchan, @(x,y)(semilogx(x,y,'r-')),1,5)%, @plot, 1, size(undersmplE,1))
%         % plot_channels(step.fmc, 0.15*ideal_template',1:n_fchan, @(x,y)(semilogx(x,y,'r-')),1,n_fchan)%, @plot, 1, size(undersmplE,1))
%         % plot_channels(step.fmc, prctile(CIrand_EPSM,CIlim(1),3)',1:n_fchan, @(x,y)(semilogx(x,y,'k:')),1,n_fchan)%, @plot, 1, size(undersmplE,1))
%         % plot_channels(step.fmc, zeros(size(CI_EPSM')), 1:n_fchan, @(x,y)(semilogx(x,y,'k-')), 1,n_fchan, num2str(step.fc'));
%         % plot_channels(step.fmc, prctile(CIrand_EPSM,CIlim(2),3)',1:n_fchan, @(x,y)(semilogx(x,y,'k:')),1,n_fchan)%, @plot, 1, size(undersmplE,1))
%         % plot_channels(step.fmc, prctile(CIboot_EPSM,CIlim(2),3)',1:n_fchan, @(x,y)(semilogx(x,y,'b:')),1,n_fchan)%, @plot, 1, size(undersmplE,1))
%         % xlabel('modulation frequency (Hz)')
%         % 
%         % figure('Name', 'Target-absent kernel');
%         % plot_channels(step.fmc, CI1_EPSM', 1:n_fchan, @(x,y)(semilogx(x,y,'b-')))%, @plot, 1, size(undersmplE,1))
%         % plot_channels(step.fmc, prctile(CI1boot_EPSM,CIlim(1),3)',1:n_fchan, @(x,y)(semilogx(x,y,'b:')))%, @plot, 1, size(undersmplE,1))
%         % plot_channels(step.fmc, zeros(size(CI1_EPSM')), 1:n_fchan, @(x,y)(semilogx(x,y,'k-')));
%         % plot_channels(step.fmc, 0.15*ideal_template',1:n_fchan, @(x,y)(semilogx(x,y,'r-')))
%         % %plot_channels(step.fmc, 60*squeeze(Template)',1:n_mchan, @(x,y)(semilogx(x,y,'r-')))
%         % plot_channels(step.fmc, prctile(CI1rand_EPSM,CIlim(2),3)',1:n_fchan, @(x,y)(semilogx(x,y,'k:')))%, @plot, 1, size(undersmplE,1))
%         % plot_channels(step.fmc, prctile(CI1rand_EPSM,CIlim(1),3)',1:n_fchan, @(x,y)(semilogx(x,y,'k:')), [], [], num2str(step.fc'))%, @plot, 1, size(undersmplE,1))%, @plot, 1, size(undersmplE,1))
%         % plot_channels(step.fmc, prctile(CI1boot_EPSM,CIlim(2),3)',1:n_fchan, @(x,y)(semilogx(x,y,'b:')))%, @plot, 1, size(undersmplE,1))
%         % xlabel('modulation frequency (Hz)')
%         % 
%         % figure('Name', 'Target-present kernel');
%         % plot_channels(step.fmc, CI2_EPSM', 1:n_fchan, @(x,y)(semilogx(x,y,'b-')))%, @plot, 1, size(undersmplE,1))
%         % plot_channels(step.fmc, prctile(CI2boot_EPSM,CIlim(1),3)',1:n_fchan, @(x,y)(semilogx(x,y,'b:')))%, @plot, 1, size(undersmplE,1))
%         % plot_channels(step.fmc, zeros(size(CI2_EPSM')), 1:n_fchan, @(x,y)(semilogx(x,y,'k-')));
%         % plot_channels(step.fmc, 0.15*ideal_template',1:n_fchan, @(x,y)(semilogx(x,y,'r-')))
%         % %plot_channels(step.fmc, 60*squeeze(Template)',1:n_mchan, @(x,y)(semilogx(x,y,'r-')))%, @plot, 1, size(undersmplE,1))
%         % plot_channels(step.fmc, prctile(CI2rand_EPSM,CIlim(2),3)',1:n_fchan, @(x,y)(semilogx(x,y,'k:')))%, @plot, 1, size(undersmplE,1))
%         % plot_channels(step.fmc, prctile(CI2rand_EPSM,CIlim(1),3)',1:n_fchan, @(x,y)(semilogx(x,y,'k:')), [], [], num2str(step.fc'))%, @plot, 1, size(undersmplE,1))
%         % plot_channels(step.fmc, prctile(CI2boot_EPSM,CIlim(2),3)',1:n_fchan, @(x,y)(semilogx(x,y,'b:')))%, @plot, 1, size(undersmplE,1))
%         % xlabel('modulation frequency (Hz)')
% 
%         %% TF CI
% 
%         foldername = cfg_game.folder_name;
%         [noise_spec,f_spec,t_spec] = noisetone_tf_converter(foldername, ListStim, data_passation.n_stim, [700 1300]);
%         %[noise_spec,f_spec,t_spec] = noise_tf_converter(foldername, ListStim, data_passation.n_stim, [700 1300]);
% 
%         [CI, ~, ~, ResponseMatrix, CI2, ~, ~, CI1, ~, ~] = computeCI(n_signal,n_responses, noise_spec, 0, 0, 'no');
%         %%
%         %fE=sqrt(fcut(1:end-1).*fcut(2:end));
% 
%         figure('Name', 'General kernel');
%         h=pcolor(t_spec,f_spec,[CI]); set(h, 'edgecolor','none'); colorbar; C_axis = caxis; caxis([-1 1]*max(abs(C_axis)))
%         % hold on
%         % plot(t_spec,50*ideal_template+fE','r')
%         saveas(gcf,'CItf_all.png')
% 
%         figure('Name', 'Target-absent kernel');
%         h=pcolor(t_spec,f_spec,[CI1]); set(h, 'edgecolor','none'); colorbar; C_axis = caxis; caxis([-1 1]*max(abs(C_axis)))
%         % hold on
%         % plot(tE,50*ideal_template+fE','r')
%         saveas(gcf,'CItf_ta.png')
% 
%         figure('Name', 'Target-present kernel');
%         h=pcolor(t_spec,f_spec,[CI2]); set(h, 'edgecolor','none'); colorbar; C_axis = caxis; caxis([-1 1]*max(abs(C_axis)))
%         % hold on
%         % plot(tE,50*ideal_template+fE','r')
%         saveas(gcf,'CItf_tp.png')
%     end
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if bSave
%     for i = 1:length(h);
%         opts = [];
%         opts.format = 'epsc'; % vectorial format for figures
%         Saveas(h(i),[dir_out hname{i}],opts);
% 
%         opts.format = 'png'; % vectorial format for figures
%         Saveas(h(i),[dir_out hname{i}],opts);
%     end
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [f,bUsingAMT] = il_audtofreq(erb)
% 
% try
%     f = audtofreq(erb);
%     bUsingAMT = 1;
% catch
%     f = ERB2f(erb);
%     bUsingAMT = 0;
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [erb,bUsingAMT] = il_freqtoaud(f)
% 
% try
%     erb = freqtoaud(f);
%     bUsingAMT = 1;
% catch
%     erb = f2ERB(f);
%     bUsingAMT = 0;
% end