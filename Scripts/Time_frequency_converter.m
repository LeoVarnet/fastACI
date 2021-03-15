function [P_dB,f_spec,t_spec] = Time_frequency_converter(dir_where, fname, N, opts)
% function [P_dB,f_spec,t_spec,E_tones,outs] = Time_frequency_converter(foldername,ListStim,n_stim, fcut)
%
% This script expects that you input all sounds you want to convert to a T-F
% representation.

% [E,f_spec,t_spec,E_tones,outs] = 

% outs = [];
% 
% SNR = -10;
% lvl = 65;
% dBFS = 100;
% % ATTENTION BRICOLAGE %
% [tone,   fs] = audioread(['.' filesep 'TargetStims' filesep 'nontarget.wav']);
% [target, fs] = audioread(['.' filesep 'TargetStims' filesep 'target.wav']);
% 
% tone   = setdbspl(tone,lvl+SNR,dBFS);
% target = setdbspl(target,lvl+SNR,dBFS);
% 
% outs.tone = tone;
% outs.target = target;

P_dB = [];

if nargin < 4
    opts = [];
end

if nargin < 3
    N = length(fname);
end

opts = Ensure_field(opts,'type','tf');
switch opts.type
    case 'tf'
        opts = Ensure_field(opts,'window','hamming');
    case 'lyon'
        path = fileparts(which('data_load'));
        path = fileparts(path); % one level up
        path = [path filesep 'tb_AuditoryToolbox' filesep];
        addpath(path);
        
        subpath = [path 'src' filesep];
        addpath(subpath);
        
        decimation = 350; % Parametres de calcul du cochleogramme
        earQ = 8;         %
        stepfactor = 0.4; %
end

L = 512;
L_overlap = 0;
L_f = L;
for i=1:N
    fprintf(['stim #' num2str(i) '\n'])
    file2load = [dir_where fname{i}];
    [noise, fs] = audioread(file2load);
    
    % S = tone + noise;
    S = noise;
    
    switch opts.type
        case 'tf'
            switch opts.window
                case 'hamming'
                    win = hamming(L);
            end
            [~,f_spec,t_spec,P] = spectrogram(S,win,L_overlap,L_f,fs);
            % P is the power spectral density
            % P=P(find(f_spec>=fcut(1) & f_spec<=fcut(end)),:,:);
            P_dB(:,:,i) = 10*log10(abs(P));
            
        case 'lyon'
            Tcoch = 1/(fs/decimation);
            t_spec=(Tcoch:Tcoch:(length(noise)/fs))-Tcoch/2;
            [~, f_spec] = DesignLyonFilters(fs,earQ,stepfactor);
            f_spec=f_spec(end:-1:1);
            
            noise_coch = LyonPassiveEar(noise,fs,decimation,earQ,stepfactor,1,1);
            noise_coch=noise_coch(end:-1:1,:);
            % noise_coch=noise_coch(cfg.freq_analysis_index,cfg.time_analysis_index);
            P_dB(:,:,i) = noise_coch;
    end
    
    %S=S/rms(S);
end


% case 'lyon'
% 

% 
%         Tcoch = 1/(fs/cfg.decimation);
%         t=(Tcoch:Tcoch:(length(bruit)/fs))-Tcoch/2;
%         [~, f] = DesignLyonFilters(fs,cfg.earQ,cfg.stepfactor);
%         f=f(end:-1:1);
%         cfg.freq_analysis_index = find(f>=cfg.freq_analysis(1) & f<=cfg.freq_analysis(2));
%         f = f(cfg.freq_analysis_index);
%         cfg.time_analysis_index = find(t>=cfg.time_analysis(1) & t<=cfg.time_analysis(2));
%         t = t(cfg.time_analysis_index);
%         size_CI=[length(f),length(t)];
%         if isempty(ROI)
%             Data_matrix = zeros(cfg.N_TrialsLoad, size_CI(1), size_CI(2));
%         else
%             Data_matrix = zeros(cfg.N_TrialsLoad, size(ROI,3));
%         end
%         tf.f=f;
%         tf.t=t;
        
switch opts.type
    case 'lyon'
        rmpath(path);
        rmpath(subpath);
end