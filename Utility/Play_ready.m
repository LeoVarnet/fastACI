function Play_ready
% function Play_ready

fname = [fastACI_basepath 'Stimuli' filesep 'ready.wav'];

[insig,fs] = audioread(fname);

sound(insig,fs);

