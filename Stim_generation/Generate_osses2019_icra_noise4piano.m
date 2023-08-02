function [outsig,outs] = Generate_osses2019_icra_noise4piano(insig,fs,version_nr,bDebug)
% function [outsig,outs] = Generate_osses2019_icra_noise4piano(insig,fs,version_nr,bDebug)
%
% 1. Description:
%       Imitate the making of the icra5 noise (Dreschler et al, 2005, Int J Aud)
%       The processing introduced here is similar to that of icra5_noise4piano
%       using method = 3.
%       This function requires as input any input signal 'insig' from which
%       a temporally and spectrally shaped waveform is generated.
%       The recommended version_nr is 0 = 8.6 (version B in Osses et al 2019).
%       The other published version is 6.0 or 6.5 which correspond to version A
%       in Osses et al 2019, however, version A of the algorithm has an
%       undersiderd spectral tilt towards high frequencies (of about +10 dB at most).
% 
%       Existing versions of our ICRA noise (versions 0,9,8,7,6,5 available in this script):
%%% icra_noise4piano:
%       version_nr 0  - 17/04/2017 (same as version_nr 8.6)
%       version_nr 9  - 17/04/2017 (testing Hohmann implementation)
%       version_nr 8.6- 17/04/2017 (introducing delay compensation for all freq. bands)
%       version_nr 8.5- 17/04/2017 (introducing delay compensation for one frequency band)
%       version_nr 8  - 06/04/2017 (similar to version 7, but adjusting band levels according to the stored rms values)
%       version_nr 7.9- 07/04/2017 (similar to version 8, but still keeping the phase randomisation)
%       version_nr 7  - 05/04/2017 - Adjusts band level after re-filtering using the theoretical BW
%       version_nr 6.5- 23/04/2017 - same as 6 but without Phase randomisation + basef to 554 Hz
%       version_nr 6  - 14/04/2016 - as used in the experiments
%       version_nr 5  - 06/04/2016
%       version_nr 4  - 02/04/2016
%%% icra5_noise4piano:
%       version_nr 3 - 18/01/2016
%       version_nr 2 - 26/12/2015 - as sent to Antoine (band levels adjusted)
%       version_nr 1 - 16/12/2015
% 
% 2. Stand-alone example:
%       file       = [Get_TUe_data_paths('db_speechmaterials') 'Spanish' delim 'Matrix' delim '00131.wav'];
%       [insig,fs] = audioread(file);
%       outsig     = icra_noise4piano_debug(insig,fs);
%       lvl        = rmsdb(insig);
%       outsig = setdbspl(outsig,lvl+100);
%       sound(insig,fs);
%       pause(length(insig)/fs*1.2);
%       sound(outsig,fs);
% 
%       file = [Get_TUe_data_paths('piano') '01-Chabassier' delim 'SONS' delim 'Cd5' delim 'pressionexpe.wav'];
%       [insig,fs] = audioread(file);
%       outsig = icra_noise4piano_debug(insig,fs);
%       lvl = rmsdb(insig);
%       outsig = setdbspl(outsig,lvl+100);
%       sound(insig,fs);
%       pause(length(insig)/fs*1.2);
%       sound(outsig,fs);
% 
%       dur = 4;
%       fs = 44100;
%       finf = 0;
%       fsup = fs/2;
%       SPL = 70;
%       % insig = AM_random_noise(finf,fsup,SPL,dur,fs); % default, no modulation, this can be useful to have a look at the auditory filter response
% 
%       insig = [zeros(4095,1); 1; zeros(4096,1)];
%       bDebug = 1;
%       version_nr = 8;
%       outsig0 = icra_noise4piano_debug(insig,fs,bDebug,version_nr);
%       version_nr = 7;
%       outsig7 = icra_noise4piano_debug(insig,fs,bDebug,version_nr); 
%       version_nr = 6;
%       outsig6 = icra_noise4piano_debug(insig,fs,bDebug,version_nr); 
% 
% 3. See also Alejandro's local scripts:
%       r20160316_update_Antoine.m, icra5_noise4piano.m
% 
% Author: Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2017
% Created on    : 21/03/2016
% Created on    : 02/08/2023 (integrated in fastACI)
% Original name : icra_noise4piano_debug.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 4
    bDebug = 0;
end

if nargin < 3
    version_nr = 0;
end

dBFS = 100; % dB FS convention used inside this script: amplitude 1 is equal to 100 dB SPL
RMS_insig = rmsdb(insig) +dBFS; % level of the input signal
switch version_nr
    case 7.9
        bwmul = input('Enter the target BW of the ERB filters in ERB-N (suggested 1, 0.5 ERB-N): '); % the overlap (above the -3 dB points) of the filters will be approximately 1/bwmul;
    otherwise
        bwmul = 1;
end

switch version_nr
    case {0,8,7.9,7,6,5}
        basef = []; % 554; % [];
        % warning('Temporal...')
        
        % auditory filterbank with 1 ERB spacing (independent of bwmul):
        if bwmul == 1 
            % [outsig, fc] = auditoryfilterbank_filtfilt(insig,fs,'bwmul',bwmul); % default AMT design
            [outsig, fc] = auditoryfilterbank(insig,fs,'bwmul',bwmul,'basef',basef); % default AMT design
        else
            % [outsig, fc] = auditoryfilterbank_filtfilt(insig,fs,'bwmul',bwmul,'betamul',bwmul); % design with narrower/wider bands
            [outsig, fc] = auditoryfilterbank(insig,fs,'bwmul',bwmul,'betamul',bwmul,'basef',basef); % design with narrower/wider bands
        end
    case 6.5
        
        basef = 554;
        if bwmul == 1 
            % [outsig, fc] = auditoryfilterbank_filtfilt(insig,fs,'bwmul',bwmul); % default AMT design
            [outsig, fc] = auditoryfilterbank(insig,fs,'bwmul',bwmul,'basef',basef); % default AMT design
        else
            % [outsig, fc] = auditoryfilterbank_filtfilt(insig,fs,'bwmul',bwmul,'betamul',bwmul); % design with narrower/wider bands
            [outsig, fc] = auditoryfilterbank(insig,fs,'bwmul',bwmul,'betamul',bwmul,'basef',basef); % design with narrower/wider bands
        end
        
    case 9
        
        error('Not tested for fastACI yet')
        basef = []; % 1000; %[];
        delay_time    = 0.036; % ms
        delay_samples = round(2*delay_time*fs);
        insig = [insig; il_Gen_silence(2*delay_time,fs)];
        
        % Construct new analyzer object;
        
        flow = 80;
        fhigh = 8000;
        analyser = gfb_analyzer_new(fs,flow,basef,fhigh,1/bwmul,4,bwmul);
        synthe   = gfb_synthesizer_new(analyser, delay_time); 
        % Filter signal;
        [outsig, analyser] = gfb_analyzer_process(analyser, insig);
        
        outsig = transpose(real(outsig));
        fc     = transpose(analyser.center_frequencies_hz);
        
        % [outsig,fc] = auditoryfilterbank(insig,fs,'bwmul',bwmul);
    case {8.6,8.5}
        
        % auditory filterbank with 1 ERB spacing (independent of bwmul):
        basef = 554;
        fc = erbspacebw(80,8000,bwmul,basef); % default fc, without basef
        [b,a,delay,z,p,k] = gammatone(fc,fs,'complex','allpole');
        fc2comp = basef-1; % input('Enter the lowest frequency you want to compensate for: ');
        idx = find(fc >= fc2comp,1,'first');
        if bDebug
            fprintf('Compensating delay for band %d, centred at %.1f Hz\n',idx,fc(idx));
        end
        delay_time = ones(size(fc));
        
        if version_nr == 8.5
            delay_time = delay(idx)*delay_time; % ms
        elseif version_nr == 8.6
            delay_time = delay;
            delay_time(1:idx-1) = delay_time(idx);
        end
        delay_samples = round(2*delay_time*fs);
        insig = [insig; il_Gen_silence(2*max(delay_time),fs)];
        
        if bwmul == 1 
            % [outsig, fc] = auditoryfilterbank_filtfilt(insig,fs,'bwmul',bwmul); % default AMT design
            [outsig, fc] = auditoryfilterbank(insig,fs,'bwmul',bwmul,'basef',basef); % default AMT design
        else
            % [outsig, fc] = auditoryfilterbank_filtfilt(insig,fs,'bwmul',bwmul,'betamul',bwmul); % design with narrower/wider bands
            [outsig, fc] = auditoryfilterbank(insig,fs,'bwmul',bwmul,'betamul',bwmul,'basef',basef); % design with narrower/wider bands
        end
        
end

Nbands  = length(fc);

switch version_nr
    case {0,9,8.6,8.5,8,7.9,7}
        BL = rmsdb(outsig)+dBFS; % Band level
        % SL = BL-10*log10(f2-f1); % Spectrum level
    case {6.5,6,5}
        % Nothing is additionally done
end

%%%
if bDebug
    fprintf('%.2f dB (Level of the input signal)\n',RMS_insig);
    interimRMS0   = rmsdb(outsig)+dBFS; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Schroeder filter: included in all versions:
for i = 1:Nbands
    outsig(:,i)  = il_schroeder(outsig(:,i));
end

%%%
if bDebug
    interimRMS1   = rmsdb(outsig)+dBFS;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Re-filtering:
switch version_nr
    case {0,8.6,8.5,8,7.9,7,6.5,6,5}
        for i = 1:Nbands
            outsig(:,i) = auditoryfilterbank(outsig(:,i),fs,fc(i),'flow',fc(i),'fhigh',fc(i),'basef',fc(i)); % compensates for the 'considerabe' change of level
        end   
    case 9
        for i = 1:Nbands
            outsig(:,i) = auditoryfilterbank(outsig(:,i),fs,fc(i),'flow',fc(i),'fhigh',fc(i),'basef',fc(i)); % compensates for the 'considerabe' change of level
        end 
        
        outsig = gfb_delay_process(synthe.delay, transpose(outsig)); % first compensation
        % outsig = transpose(outsig.*repmat(synthe.mixer.gains',1,size(outsig,2)));
        outsig = transpose(outsig);
end
BL_filt = rmsdb(outsig)+dBFS; % Band level of the re-filtered signals

switch version_nr
    case 7
        fcut_low  = audtofreq(freqtoaud(fc)-0.5); % Hz
        fcut_high = audtofreq(freqtoaud(fc)+0.5); % Hz
        BL_BW     = 10*log10(fcut_high-fcut_low); % band level of a noise assuming flat amplitude
                                                  % (which is the case after il_schroeder).
        BL_full   = 10*log10(fs/2-20); % audible range
        BL_corr   = BL_full-BL_BW; 
    otherwise
        % Nothing is additionally done at this stage
end

for i = 1:Nbands
    switch version_nr
        case {0,9,8.6,8.5,8,7.9}
            outsig(:,i) = From_dB( BL(i)-BL_filt(i) )*outsig(:,i);
            % The previous line is equivalent to: outsig(:,i) = setdbspl(outsig(:,i),BL(i));
        case 7
            outsig(:,i) = From_dB( BL_corr(i) )*outsig(:,i); % compensates for the 'considerabe' change of level
        case {6.5,6,5}
            % Nothing is additionally done
    end
end

if bDebug
    interimRMS2   = rmsdb(outsig)+dBFS;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adding signals together: included in all versions:
switch version_nr
    case 9
        % synthesizer.mixer.gains = ones(size(synthesizer.mixer.gains));
        % Resynthesize filtered impulse response from above.
        % outsig = gfb_delay_process(synthe.delay, transpose(outsig));
        % outsig = transpose(outsig);
        % synthe.mixer.gains = ones(size(synthe.mixer.gains));
        outsig = gfb_synthesizer_process(synthe, transpose(outsig));
        outsig = transpose(outsig); 
        outsig = outsig(delay_samples+1:end,1);
        outsig = From_dB(3)*outsig; % 3 dB, from empirical observation
        
        % outsig = setdbspl(outsig,interimRMStot);
        % figure; plot(insig,'r'); hold on; plot(outsig); 
        % outsig3 = sum(outsig,2); % figure; plot(insig,'r'); hold on; plot(outsig3); 
    case 8.5
        outsig = outsig(delay_samples+1:end,:);
        outsig = sum(outsig,2);
    
    case 8.6
        M = size(outsig,1);
        for i = 1:Nbands
            outsig(1:M-delay_samples(i),i) = outsig(delay_samples(i)+1:end,i);
        end
        outsig(end-max(delay_samples)+1:end,:) = [];
        outsig = sum(outsig,2);
        
    otherwise
        outsig = sum(outsig,2); % Just to check some levels
end

if bwmul < 1
    disp([mfilename ': the centre frequencies of the filter bank are spaced '])
    disp(           ' in less than 1 ERB while the filters, while the bandwidth is approximately 1 ERB.')
    disp(           ' the output signal may require a level scaling (not applied here).')
end
if bwmul > 1
    disp([mfilename ': the centre frequencies of the filter bank are spaced in more'])
    disp(           ' than 1 ERB while the bandwidth of the filters is approximately 1 ERB.')
    disp(           ' the output signal may require a level scaling (not applied here).')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch version_nr
    case {0,9,8.5,8,7}
        % Nothing is additionally done

    case 7.9
        N = input('Enter the length N of the overlap-add filter (suggested values of 512, 4096, 8192): ');
        % bDoRandomisation = 1; % default
        % bDoRandomisation = 0; % no randomisation, only overlap-add
        bDoRandomisation = 1; % randomisation without overlap-add
        outsig    = il_randomise_phase(outsig,N,bDoRandomisation);
        if bDoRandomisation ~= 1
            display('Phase randomisation will be done in an ALTERNATIVE way. Press any button to continue, press Ctrl+C if this is by mistake...')
            pause()
        end
        outsig    = setdbspl( outsig,RMS_insig ); % this should be a small adjustment only
    case 6.5
        % No phase randomisation, no IIR filter (compared to v. 6)
        outsig    = setdbspl( outsig,RMS_insig ); 
    case 6
        outsig    = il_randomise_phase(outsig,512); % it can (slightly) increase or decrease the level
        outsig    = il_apply_IIR_Butter(outsig,fs,audtofreq(33.5),'low',8);
        outsig    = scaletodbspl( outsig,RMS_insig,'dboffset',dBFS); % this should be a small adjustment only
    case 5
        outsig    = il_randomise_phase(outsig,512); % it can (slightly) increase or decrease the level
        outsig    = il_apply_IIR_Butter(outsig,fs,audtofreq( 2.5),'high',4);
        outsig    = il_apply_IIR_Butter(outsig,fs,audtofreq(33.5), 'low',8);
        outsig    = scaletodbspl( outsig,RMS_insig,'dboffset',dBFS); % this should be a small adjustment only
end

if bDebug
    outsig_tmp   = auditoryfilterbank(outsig,fs,'bwmul',bwmul,'basef',basef);
    interimRMS3  = rmsdb(outsig_tmp)+dBFS;
end

%%%
if bDebug
    fprintf('Levels after Gammatone filterbank (col 1)\n')
    fprintf('Levels after Schroeder randomisaton (col 2)\n')
    switch version_nr
        case {0,8,7}
            fprintf('Levels after re-filtering, with band level adjustment (col 3)\n')
        otherwise
            fprintf('Levels after re-filtering (col 3)\n')
    end
    fprintf('Levels after phase randomisaton (converted back to band levels, col 4)\n')
    disp( [interimRMS0' interimRMS1' interimRMS2' interimRMS3'] );
end

outs.fc = fc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function outsig=il_schroeder(insig)
outsig=insig.*sign(rand(length(insig),1)-0.5);

function outsig = il_randomise_phase(insig,N,bDoRandomisation)
% function outsig = il_randomise_phase(insig,N,bDoRandomisation)
%
% Taken from Dreschler et al 2001:
% The phases are modified in a 512-point FFT procedure. This was accomplished
% by randomising the phase and overlapping/adding the segments after an inverse FFT.

if nargin < 3
    bDoRandomisation = 1; % this is the default
end

hop    = N/8;
window = hamming(N);
outsig = zeros(length(insig),1);

switch bDoRandomisation
    case {0,1}
        for i=0:round(length(insig)/hop)-1-N/hop

            from = i*hop+1;
            to   = i*hop+N;
            stuk0 = insig(from:to);
            stuk1 = fft(stuk0.*window,N);
            if bDoRandomisation
                % Version as it was...
                stuk2 = stuk1.*exp(j*2*pi*rand(N,1)); % filter with unit amplitudes 
                                                      % and uniformly distributed phases from 0 to 2*pi
                stuk3 = real(ifft(stuk2));
                disp('')
            else
                % If no randomisation, stuk3 should be equal to stuk 0
                stuk2 = stuk1;
                stuk3 = real(ifft(stuk2)); % stuk3 is a windowed version of stuk0

                disp(''); 
                % figure; plot(stuk0-stuk3)
                % figure; plot(stuk0); hold on; plot(stuk3,'r')
            end

            outsig(from:to)=outsig(from:to)+stuk3;
        end
    case -1
        stuk1 = fft(insig);
        stuk2 = stuk1.*exp(j*2*pi*rand(length(insig),1)); % filter with unit amplitudes 
                                              % and uniformly distributed phases from 0 to 2*pi
        stuk3 = real(ifft(stuk2));
        outsig= stuk3;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function outsig = il_apply_IIR_Butter(insig,fs,fc,type,order)
% function outsig = il_apply_IIR_Butter(insig,fs,type,fc,order)
%
% 1. Description:
%       Designs and apply an IIR Butterworth filter to a signal 'insig'.
%       type should be 'low' or 'high'

if nargin < 5
    order = 4; % slope order*6 dB/Oct
end

wc = fc/(fs/2); % normalised frequency (fs/2 = 1)
[b, a] = butter(order,wc,type);

outsig = filtfilt(b,a,insig);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function silence = il_Gen_silence(duration, fs)
% function silence = il_Gen_silence(duration, fs)
%
% 1. Description: 
%       Generates a silence with a specific duration at a specific sample
%       rate (fs).

N = round(duration * fs); % Number of time samples.
silence = zeros(N,1);     % Column vector.