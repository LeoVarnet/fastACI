function data = publ_osses2019a_JASA_figs(varargin)
% function data = publ_osses2019a_JASA_figs(varargin)
%
% 1. Description: Add a descripotion
%
% % To display Fig. 3ab of Osses, Kohlrausch and Chaigne, (2019, JASA) use :::
%     publ_osses2019a_JASA_figs('fig3ab'); % Piano P1 (anechoic) + ICRA noise representation 
%
% % To display Fig. 3de of Osses, Kohlrausch and Chaigne, (2019, JASA) use :::
%     publ_osses2019a_JASA_figs('fig3de'); % Piano P3 (anechoic) + ICRA noise representation 
%
% Author: Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% close all, clc

if nargin == 0
    help publ_osses2022b_JASA_figs;
    return
end

data = [];
h = [];
hname = [];
 
definput.flags.type={'missingflag', ...
    'fig3ab', ... % old name: fig4_simu
    'fig3de'}; % different lambdas
definput.flags.plot={'plot','no_plot'};
% definput.keyvals.dir_out=[];
 
[flags,keyvals]  = ltfatarghelper({},definput,varargin);

if flags.do_fig3ab
    Piano_ID = '1';
end
if flags.do_fig3de
    Piano_ID = '3';
end
if flags.do_fig3ab || flags.do_fig3de
    version_nr = 6.0; % version A is the same as version 6.0
end

if flags.do_fig3ab || flags.do_fig3de
    dir_stim = [fastACI_basepath 'Stimuli' filesep 'osses2019' filesep];
    
    files = Get_filenames(dir_stim,['dataset-1-' Piano_ID '-*.wav']);
    dBFS = 100; % prior knowledge about the stored waveforms
    
    [insig,fs] = audioread([dir_stim files{1}]);
    insig = From_dB(dBFS-94) * insig; % converting to pressure in Pa
    dBFS = 94; % this is the new dB FS convention
    
    switch version_nr
        case 6.0
            Noise_ID = 'A';
        case 8.6 
            Noise_ID = 'B'; % This is the recommended ICRA noise version
    end
    noiseA = Generate_osses2019_icra_noise4piano(insig,fs,version_nr);
    
    fc = 20; % Hz, for LPF applied to the Hilbert envelope
    insig_env = il_Get_envelope_piano(insig,fs,fc);
    
    noiseA_env = il_Get_envelope_piano(noiseA,fs,fc);
    
    t = (1:length(insig))/fs;
    
    re = 2e-5; % 20 uPa, reference
    insig_dB = 20*log10(abs(insig)/re); % in dB SPL
    insig_env_dB = 20*log10(abs(insig_env)/re);
    
    figure;
    subplot(1,2,1)
    plot(t,insig_dB,'r'); hold on;
    plot(t,insig_env_dB,'k-','LineWidth',2);
    grid on;
    ylim([30 90])
    xlim([0 max(t)]);
    
    XT = 0.1:.3:max(t);
    set(gca,'XTick',XT);
    YT = 40:10:80;
    set(gca,'YTick',YT);
    
    ylabel('Amplitude (dB SPL)');
    xlabel('Time (s)');
    
    title(sprintf('Piano P%s',Piano_ID));
    %%% Now we plot the noise:
    noiseA_dB     = 20*log10(abs(noiseA)    /re);
    noiseA_env_dB = 20*log10(abs(noiseA_env)/re);
    
    subplot(1,2,2)
    plot(t,noiseA_dB,'r'); hold on;
    plot(t,noiseA_env_dB,'k-','LineWidth',2);
    grid on;
    ylim([30 90])
    xlim([0 max(t)]);
    
    set(gca,'XTick',XT);
    set(gca,'YTick',YT);
    
    ylabel('Amplitude (dB SPL)');
    xlabel('Time (s)');
    
    title(sprintf('Noise N%s (version %s)',Piano_ID,Noise_ID));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function outsig = il_Get_envelope_piano(insig,fs,fc)
% function outsig = il_Get_envelope_piano(insig,fs,fc)
%
% 1. Description:
%       Applies a 4th-order Butterworth filter centred at 20-Hz.
% 
% 2. Stand-alone example:
%       outsig = Get_envelope_piano(insig,fs,20)
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Created on    : 24/11/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    fc = 20;
end

yin    = abs(hilbert( insig ));
[b, a] = butter(4,fc/(fs/2),'low');

yenv   = filtfilt(b,a,yin); % linear phase implementation

outsig = yenv; % yin + yout;
