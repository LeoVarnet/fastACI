function h = Get_raster_plot(insig,fs)
% function h = Get_raster_plot(insig,fs)
%
% h is the figure handle
%
% Author: Alejandro Osses
% Date: 4/04/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numTot    = size(insig,2);
N_samples = size(insig,1);

if nargin < 2
    t = 1:N_samples;
    xlab = 'Time (samples)';
else
    t = (1:N_samples)/fs;
    xlab = 'Time (s)';
end

for i=1:numTot
    stairs(t,.8*insig(:,i)+(i-1)*1,'k'); hold on; grid on
    % the amplitude of 0.8 is arbitrary assuming that the pulse trains are 
    %   arranged in a matrix 'insig' of 0s or 1s only.
end
xlabel(xlab);
ylabel('Spike trains per simulated neuron');
    
h = gcf;