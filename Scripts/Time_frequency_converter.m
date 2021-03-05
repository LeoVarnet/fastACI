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

opts = Ensure_field(opts,'window','hamming');

L = 512;
L_overlap = 0;
L_f = L;
for i=1:N
    fprintf(['stim #' num2str(i) '\n'])
    file2load = [dir_where fname{i}];
    [noise, fs] = audioread(file2load);
    
    % S = tone + noise;
    S = noise;
    
    switch opts.window
        case 'hamming'
            win = hamming(L);
    end
    
    [~,f_spec,t_spec,P] = spectrogram(S,win,L_overlap,L_f,fs);
    % P is the power spectral density
    
    % P=P(find(f_spec>=fcut(1) & f_spec<=fcut(end)),:,:);
    P_dB(:,:,i) = 10*log10(abs(P));
        
    %S=S/rms(S);
end

% [~,f_spec,t_spec,P] = spectrogram(tone,512,0,512,fs);
% P=P(find(f_spec>=fcut(1) & f_spec<=fcut(end)),:,:);
% pdB=10*log10(abs(P));
% E_tones(:,:,1)=P;
% 
% [~,f_spec,t_spec,P] = spectrogram(target,512,0,512,fs);
% P=P(find(f_spec>=fcut(1) & f_spec<=fcut(end)),:,:);
% pdB=10*log10(abs(P));
% E_tones(:,:,2)=P;
% 
% f_spec = f_spec(find(f_spec>=fcut(1) & f_spec<=fcut(end)));