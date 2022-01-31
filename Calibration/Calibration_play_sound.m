function Calibration_play_sound(fname,dur_playback)
% function Calibration_play_sound(fname,dur_playback)

if nargin == 0
    fname = [fastACI_basepath 'Calibration' filesep 'calibration-18dBFS.wav'];
    fprintf('%s: The default calibration sound (%s) will be played\n',upper(mfilename),fname);
end

if nargin < 2
    dur_playback = 60; % s
end

[insig,fs] = audioread(fname);

dur_signal = length(insig)/fs;

N_loops = ceil(dur_playback/dur_signal);

insig = repmat(insig,N_loops,1);

fprintf('During the reproduction, type: clear insig, if you want to stop the sound reproduction...\n')
sound(insig,fs);

