function [ stim ] = generate_stim( signal, bruit, RSB, fadein_ech, color_noise)
%GENERATE_STIM(signal, bruit, RSB, fadein_ech) cree un stim � partir des
%vecteurs signal et bruit, avec un certain RSB, en ajoutant un fadein/out
%gaussien de fadein_ech (en nbr d'echantillons)
if fadein_ech>0
    fadein = gausswin(2*fadein_ech);
    fadein = fadein(1:floor((length(fadein)/2)));
            if strcmp(color_noise, 'white')
            bruit_fadein = randn(size(fadein));
            bruit_fadeout = randn(size(fadein));
        elseif strcmp(color_noise, 'pink')
            bruit_fadein = pinknoise(size(fadein));
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
[stimRSB, A] = Addition_RSB(signal,bruit,RSB);
stim = [bruit_fadein; stimRSB; bruit_fadeout];

end

