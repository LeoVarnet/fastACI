function [stim,struct_extra] = generate_stim( signal, bruit, SNR, fadein_ech, noise_type)
%GENERATE_STIM(signal, bruit, SNR, fadein_ech) cree un stim a partir des
% vecteurs signal et bruit, avec un certain RSB, en ajoutant un fadein/out
% gaussien de fadein_ech (en nbr d'echantillons)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bruit = bruit/std(bruit);

if fadein_ech>0
    fadein = gausswin(2*fadein_ech);
    fadein = fadein(1:floor((length(fadein)/2)));
    if strcmp(noise_type, 'white')
        bruit_fadein  = randn(size(fadein));
        bruit_fadeout = randn(size(fadein));
    elseif strcmp(noise_type, 'pink')
        error('%s: Not validated yet...',upper(mfilename))
        bruit_fadein  = pinknoise(size(fadein));
        bruit_fadeout = pinknoise(size(fadein));
    end
    bruit_fadein = randn(size(fadein)); bruit_fadein = bruit_fadein/rms(bruit_fadein);
    bruit_fadein = bruit_fadein.*fadein; 
    bruit_fadeout = randn(size(fadein)); bruit_fadeout = bruit_fadeout/rms(bruit_fadeout);
    bruit_fadeout = bruit_fadeout.*fadein(end:-1:1); 
else
    bruit_fadein=[];
    bruit_fadeout=[];
end

% Noise level is kept constant, signal is adjusted to get the desired SNR:
[stimRSB, A, extra] = Addition_RSB(signal,bruit,SNR); 

% extra contains extra.N (noise after the SNR cal) and extra.S (signal after the SNR cal)

stim_N = [bruit_fadein; extra.N_scaled; bruit_fadeout];
sil = zeros(size(bruit_fadein));
stim_S = [sil; extra.S_scaled; sil];
struct_extra.stim_N = stim_N;
struct_extra.lvl_N_dBFS  = 20*log10(rms(extra.N_scaled));

struct_extra.stim_S = stim_S;
struct_extra.lvl_S_dBFS  = 20*log10(rms(extra.S_scaled));

stim = [bruit_fadein; stimRSB; bruit_fadeout];