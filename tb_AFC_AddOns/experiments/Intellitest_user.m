function IntellitestConstant_user
% function IntellitestConstant_user
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global def
global work
global setup

MaterialPath = intellitest_getpaths('material_location'); % Get_UGent_data_paths('VlMatrix'); 

if ( def.interleaved && def.OLSA_SplitTestlistAcrossTracks )
	% take every other sentence in interleaved tracks with split testlist
	currentSentence = work.pvind + (work.stepnum{work.pvind}(end)-1)*def.interleavenum;
else
	currentSentence = work.stepnum{work.pvind}(end);
end

try
    wavName = [MaterialPath setup.testlistStrings{work.pvind}{currentSentence} '.wav'];
    [sentence,fs] = audioread(wavName);
catch
    disp('')
end

% check samplerate
if ( fs ~= def.samplerate )
    error('wav files have wrong samplerate')
end

% DEBUG: print current rms of OLSA sentence 
%       20*log10(rms(sentence))
% DO NOT NORMALIZE INDIVIDUALLY, sentences have to be assumed to have -24 dB rms re full scale.

% we might want to change the track mapping
%def.trackMap =[0 0 0 1];
if def.continuous == 1
    afc_sound('bgloopwav_restart');
    bAddNoise = 0;
else
    bAddNoise = 1;
    dur_onset = 500e-3;
    dur_onset_samples = dur_onset*def.samplerate;
    
    N     = length(sentence)+2*dur_onset_samples;
    noise = il_Randomise_insig(work.noise,N);
    
    Nstart = 1+dur_onset_samples;
    Nend   = Nstart+size(sentence,1)-1;

    Nwin                = length(setup.win_up);
    noise(1:Nwin)       = setup.win_up.*noise(1:Nwin); % ramp up
    noise(end-Nwin:end) = setup.win_dn.*noise(end-Nwin:end); % ramp down
end

if work.varspeech % scale the signal and take left channel to make it mono
    sentence = sentence(:,1) * 10^(work.expvaract/20);
end
if work.varnoise
    noise    = noise(:,1) * 10^(-work.expvaract/20);
end

% pre-, post- and pausesignals (all zeros)

presig = zeros(def.presiglen,def.outputChannels);
postsig = zeros(def.postsiglen,def.outputChannels);
pausesig = zeros(def.pauselen,def.outputChannels);

% make required fields in work
if bAddNoise == 0
    % as in bglooped or speech-in-silence
    work.signal = [ sentence sentence ];	% left = right (diotic) first two columns holds the test signal (left right)
elseif bAddNoise == 1
    work.signal = [noise noise];
    work.signal(Nstart:Nend,:) = work.signal(Nstart:Nend,:)+[sentence sentence];
end
work.presig = presig;					% must contain the presignal
work.postsig = postsig;					% must contain the postsignal
work.pausesig = pausesig;				% must contain the pausesignal

% eof

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function outsig = il_Randomise_insig(insig,N)
% function outsig = il_Randomise_insig(insig,N)
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Created on    : 12/08/2015
% Last update on: 21/12/2016 
% Last use on   : 21/12/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Nstart = round( length(insig)*random('unif',[0 1]) );
Nstart = round( length(insig)*random('unif',0,1) );
Nstart = max(1,Nstart);

outsig = [insig(Nstart:end,:); insig(1:Nstart-1,:)];

if nargin >= 2
    outsig = outsig(1:N);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end