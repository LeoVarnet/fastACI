function [b, a] = vandorpschuitman2013a_outmiddlefilter(fs,order)
%vandorpschuitman2013a_outmiddlefilter simulates the outer- and middle ear transfer function from van Dorp et. al. 2013
%
%   Usage: [b, a] = vandorpschuitman2013a_outmiddlefilter(insig,fs,order)
%
%   Input parameters:
%        fs     : sampling rate.
%        order  : order of the high- and low-pass filter. The final order 
%                 of the combined outfmiddlefilter will be 2*order.
%
%   `vandorpschuitman2013_outmiddlefilter(fs,order)` returns the IIR filter
%   coefficients of a filter with a transfer function as the combined outer
%   and middle ear using a sampling frequency *fs* Hz. The first row of b 
%   and a contains the coefficients for the high-pass filter (fcut=1000 Hz)
%   and the second row contains the coefficients for the low-pass filter
%   (fcut=4000 Hz).
%
%   References: breebaart2001a vandorp2013a osses2017a
%

%   AUTHOR: Alejandro Osses

freqs=[1000 4000]; % Hz, fixed cut-off frequencies

if nargin < 2
    order = 1; % butterworth: slope=20*order/decade or 6 dB/oct
end

[B1,A1]=butter(order, freqs(1)/fs*2,'high');
[B2,A2]=butter(order, freqs(2)/fs*2,'low');

[H1,F]=freqz(B1,A1,[],fs);
[H2  ]=freqz(B2,A2,[],fs);

Ht_max = max(abs(H1.*H2));
B1 = sqrt(1/Ht_max)*B1;
B2 = sqrt(1/Ht_max)*B2;

if nargout == 0
    H1=freqz(B1,A1,[],fs);
    H2=freqz(B2,A2,[],fs);
    figure
    semilogx(F,20*log10([abs(H1) abs(H2)])); grid on
    title(sprintf('fs=%.0f [Hz]',fs))
end

b = [B1; B2];
a = [A1; A2];

end
