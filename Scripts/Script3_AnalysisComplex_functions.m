function data = Script3_AnalysisComplex_functions(cfg_game,data_passation,fct)
% function data = Script3_AnalysisComplex_functions(file_savegame,data_passation,fct)
%
% List of function processing in this script:
%
%   prctile: preselects the trials having an data_passation.expvar values 
%       between percentiles 5 and 95.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    clc
    help Script3_AnalysisComplex_functions;
    return;
end

h = [];
hname = [];

if isfield(cfg_game,'Subject_ID')
    Subject_ID = cfg_game.Subject_ID;
else
    Subject_ID = '';
end

switch lower(fct)
    case 'prctile'
        data = il_get_percentile(data_passation);
        return;
        
    case {'histogram-leo','histogram'}
        
        %%%
        % if isfield(cfg_game,'N')
        %     N_total = cfg_game.N; 
        % else
        %     N_total = cfg_game.N_presentation*cfg_game.N_target; 
        % end
        %%%
        
        expvar = data_passation.expvar;
        is_correct = data_passation.is_correct;
        n_responses = data_passation.n_responses;
        RT = data_passation.response_time; % (1:N_trials);
        n_targets   = data_passation.n_targets; %(1:N_trials);
        
        if isfield(cfg_game,'expvar_description')
            unit   = cfg_game.expvar_description;
        else
            unit = 'm (dB)';
            warning('Default x-axis variable label');
        end
        
        resp_if_2 = 2; % 'if target' (AM)       - 2 modulated tone
        resp_if_1 = 1; % 'if reference' (no AM) - 1 pure tone
    
        bPrint = 0;
        [data,PercL,PercU,expvar,outs] = il_get_percentile(data_passation,cfg_game,bPrint);
        expvar_L = data.expvar_L;
        expvar_U = data.expvar_U;
        if outs.bTransform == 0
            % Regular case
            bin_step = 1;
        else
            % If bump noises
            bin_step = 0.25;
        end
        bin_centres   = expvar_L: bin_step : expvar_U;  % linspace(expvar_L,expvar_U,9);
        
        switch lower(fct)
            case 'histogram-leo'
                % bin_centres = linspace(prctile(expvar,SNRprctile(1)),prctile(expvar,SNRprctile(2)),10);
                bin_edges   = mean([bin_centres(1:end-1);bin_centres(2:end)]); 
                bin_edges   = [bin_edges(1)-bin_step bin_edges bin_edges(end)+bin_step];
            case 'histogram'
                % bin_centres = -15.5:1:-0.5;

                % [N_m,m_edge] = histcounts(m, nbins); % Only in newer MATLAB versions
                % N_m = hist(m,m_bin); m_edge = m_bin; % In older MATLAB versions
                [N_m_all,bin_edges] = my_hist(expvar,bin_centres); % this really provides edges
                bin_edges(1)   = bin_centres(1)-bin_step/2;
                bin_edges(end) = bin_centres(end)+bin_step/2;
                N_m = histc(expvar,bin_edges); % this really provides edges
        end
        
        for i_m = 1:length(bin_centres)
            
            idx_m   = find(expvar>=bin_edges(i_m) & expvar<bin_edges(i_m+1));
            N_hist(i_m) = length(idx_m);

            N_correct(i_m) = sum(is_correct(idx_m));
            
            N_if_2(i_m)  = sum(n_targets(idx_m)==resp_if_2); 
            N_if_1(i_m)  = sum(n_targets(idx_m)==resp_if_1);

            N_resp_1(i_m) = sum(n_responses(idx_m)==resp_if_1);
            N_resp_2(i_m) = sum(n_responses(idx_m)==resp_if_2);
            
            H(i_m)  = sum(n_targets(idx_m)==resp_if_2 & n_responses(idx_m)==resp_if_2); % Hit
            M(i_m)  = sum(n_targets(idx_m)==resp_if_2 & n_responses(idx_m)==resp_if_1); % Miss
            CR(i_m) = sum(n_targets(idx_m)==resp_if_1 & n_responses(idx_m)==resp_if_1); % Correct rejection
            FA(i_m) = sum(n_targets(idx_m)==resp_if_1 & n_responses(idx_m)==resp_if_2); % False alarm
        end
        H_tot  = sum(H); 
        M_tot  = sum(M);
        CR_tot = sum(CR);
        FA_tot = sum(FA);
        tot_classified    = H_tot + CR_tot + M_tot + FA_tot;
        responses_counted = sum(N_hist);
        
        % % Sanity check (to know that all responses were processed):
        % if responses_counted ~= N_total
        %     error('%s: Not all responses were processed. Check whether the participant indeed completed the whole session...',upper(mfilename))
        % end
        % if tot_classified ~= N_total
        %     error('%s: Not all responses were classified as H, M, CR, or FA. Check whether this is correct. If yes, convert this message into a warning only...',upper(mfilename))
        % end
        
        idxs_excluded = find(expvar < bin_edges(1) | expvar >= bin_edges(end));
        responses_excluded = length(idxs_excluded);
        
        prop_correct      = N_correct./N_hist;
        bias              = N_resp_1./N_hist;
        bias_description  = 'bias towards response 1';
        
        figure; 
        yyaxis left
        area(bin_centres,N_hist,'FaceAlpha',0.2,'EdgeColor','none');
        
        Location = 'NorthWest';
        if strfind(lower(cfg_game.expvar_description),'bump')
            str4label = '(log2[number of bumps])';
            Location = 'SouthWest';
        elseif strfind(lower(cfg_game.expvar_description),'snr')
            str4label = '(SNR in dB)';
        else
            str4label = '';
        end
        
        xlabel(['expvar ' str4label]);
        
        ylabel('N of trials')
        yyaxis right
        
        hprop = plot(bin_centres,prop_correct,'-o');
        ylabel('proportion'); 
        ylim([0 1]); hold on
        
        hbias = plot(bin_centres,bias, '--o');
        ylabel('prop'); 
        ylim([0 1])
        grid on
        legend([hprop, hbias],{'correct','bias (resp. was 1)'},'Location',Location)
        
        h(end+1) = gcf;
        hname{end+1} = [Subject_ID '-Confusion'];
        
        % yyaxis right
        % hprop = plot(SNRbin_centres,prop_correct,'-o');ylabel('prop'); ylim([0 1]); hold on
        % hbias = plot(SNRbin_centres,bias, '--o');      ylabel('prop'); ylim([0 1])
        % grid on
        % legend([hprop, hbias],{'correct','bias'},'Location','northwest')
        
        %%%
        % -----------------------------------------------------------------
        % --- First figure: Behavioural results
        n_window = 100;
        minNformean = 50;
        
        % Memory allocation:
        N_windows = length(n_responses)/n_window;

        m_windowed       = nan(1,N_windows);
        bias_windowed    = nan(1,N_windows);
        PC_targetpresent = nan(1,N_windows);
        PC_targetabsent  = nan(1,N_windows);
        RT_windowed      = nan(1,N_windows);
        
        for i = 1:N_windows
            idxs_here = (i-1)*n_window+1:i*n_window; % indexes of the trials within each window
            response_windowed  = n_responses(idxs_here);
            signal_windowed    = n_targets(idxs_here);

            m_windowed(i)      = mean(expvar(idxs_here));
            bias_windowed(i)   = mean(response_windowed);
            PC_targetpresent(i)= mean(response_windowed(signal_windowed==2))-1;
            PC_targetabsent(i) = 2-mean(response_windowed(signal_windowed==1));
            RT_windowed(i)     = mean(RT(idxs_here));
        end
        
        N_trials = length(expvar);
        % ---
        m_presentation_nr     = 1:N_trials;
        m_presentation_nr_win = 1:n_window:N_trials;

        figure('Position', [100 100 800 500]); 
        subplot(2,2,1); 
        plot(m_presentation_nr    , expvar         ,'g'); hold on; 
        plot(m_presentation_nr_win, m_windowed,'k'); 
        ylim([bin_edges(1) bin_edges(end)]); 
        xlabel(' trial #'); ylabel(unit); 
        xlim([1 length(expvar)]); ylimits=ylim;

        N_sessions = length(data_passation.resume_trial);
        for i = 1:N_sessions
            plot(data_passation.resume_trial(i)*[1 1],ylimits,'k:');
        end
        title(sprintf('expvar per trial (of %.0f trials)',length(expvar)))

        % ---
        subplot(2,2,3); 
        plot(m_presentation_nr_win,PC_targetpresent); hold on; 
        plot(m_presentation_nr_win,PC_targetabsent);  
        xlim([1 length(expvar)]); xlabel(' trial #'); 
        ylabel('correct response rate'); ylim([0 1]); hold on; 
        plot([1 length(expvar)],[0.5 0.5],'k--'); 
        ylimits=ylim;
        title(sprintf('CR rate (bins of %.0f trials)',n_window))

        for i = 1:length(data_passation.resume_trial)
            % Vertical dotted lines at the points where a new session was started:
            plot(data_passation.resume_trial(i)*[1 1],ylimits,'k:','LineWidth',2);
        end

        subplot(2,2,2); 
        var2plot = [H', M', CR', FA'];
        bar(bin_centres, var2plot); 
        xlim([bin_edges(1) bin_edges(end)]); 
        xlabel(unit); 
        ylabel('Number of trials'); hold on; 
        plot([bin_edges(1) bin_edges(end)],[minNformean minNformean]/2,'k:');
        legend({'H', 'M', 'CR', 'FA', 'Nmin'},'Location','best');
        title(sprintf('Histogram (of %.0f trials)',sum(sum(var2plot))))
        
        subplot(2,2,4); 
        bar(bin_centres, [H'./(M'+H'), CR'./(CR'+FA')].*[M'+H'>minNformean,CR'+FA'>minNformean]);
        xlim([bin_edges(1) bin_edges(end)]); 
        xlabel(unit); 
        ylabel('correct response rate'); 
        hold on; 
        plot([bin_edges(1) bin_edges(end)],[0.5 0.5],'k--'); 
        hl = legend({sprintf('target %.0f present',resp_if_2), sprintf('target %.0f absent',resp_if_2),'chance level'},'Location','southeast');
        set(hl,'FontSize',8);
        title(sprintf('CRs (of %.0f trials)',tot_classified))
        
        h(end+1) = gcf;
        hname{end+1} = [Subject_ID '-Behaviour'];

        %%%
        H_rate  = 100*H./N_if_2; % /N_tot_signal;
        CR_rate = 100*CR./N_if_1; % /N_tot_absent;

        figure;
        plot(bin_centres,H_rate,'bo-'); hold on;
        plot(bin_centres,CR_rate,'rs--','LineWidth',2);
        ylim([-3 103]);
        set(gca,'YTick',0:5:100); grid on
        
        % xlim([-17 1])
        % set(gca,'XTick',-16:0);

        legend('H rate','CR rate','Location','best');
        ylabel(sprintf('Percentage correct\n(ref. # of presentations per m interval)'));
        xlabel(['Central bin of the tested ' unit]);

        h(end+1) = gcf;
        hname{end+1} = [Subject_ID '-Hit-and-CR-rates-per-bin'];        
        
        data.trialnum = m_presentation_nr_win;
        data.m_windowed = m_windowed;
        
        data.bin_edges = bin_edges;
        data.H = H;
        data.M = M;
        data.CR = CR;
        data.FA = FA;
        
        data.h = h;
        data.hname = hname;
        return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data,PercL,PercU,expvar,outs] = il_get_percentile(data_passation,cfg_game,bPrint)

outs = []; 
if nargin < 3
    bPrint = 1;
end
expvar = data_passation.expvar;

bTransform = 0;
if isfield(cfg_game,'Condition')
    switch cfg_game.Condition
        case 'bump'
            bTransform = 1;
        otherwise
            bTransform = 0;
    end         
end

if bTransform
    expvar = -expvar;
    
    outs.exp2eval_transformation = 'expvar = log10(expvar)/log10(2);';
    eval(outs.exp2eval_transformation);
    
    idxs = find(expvar<0);
    expvar(idxs) = 0;
end

outs.bTransform = bTransform;

PercL = 5;
PercU = 95;
expvar_L = prctile(expvar,PercL);
expvar_U = prctile(expvar,PercU);

idxs_select_boolean = expvar>=expvar_L & expvar<= expvar_U; 
idxs_select =    find(expvar>=expvar_L & expvar<= expvar_U);

data.expvar_L = expvar_L;
data.expvar_U = expvar_U;
data.idxs_select_boolean = idxs_select_boolean;
data.idxs_select_boolean_description = '0 or 1 if the corresponding expvar is not selected or selected, respectively';
data.idxs_select = idxs_select;
data.idxs_select_description = 'Indexes between percentiles 5 and 95 of expvar';

if bPrint == 1
    fprintf('prctile preproc: Selected idxs=%.0f of %.0f trials\n',length(idxs_select),length(expvar));
end