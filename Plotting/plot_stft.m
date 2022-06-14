function plot_stft(t,f,insig_stft,opts)
% function plot_stft(t,f,insig_stft,opts)

if nargin < 4
    opts = [];
end
opts = Ensure_field(opts,'bColourBar',1);

% % See stft_example for explanations about the following lines...
% K       = sum(win)/wlen;
% stft    = abs(stft)/wlen/K;
% 
% if rem(nfft, 2) % odd nfft excludes Nyquist point
%     stft(2:end, :)      = stft(2:end, :).*2;
% else            % even nfft includes Nyquist point
%     stft(2:end-1, :)    = stft(2:end-1, :).*2;
% end
% 
% stft = 20*log10(stft + 1e-6);

% if isfield(opts,'dB_min')
%     idx = find(stft < info.dB_min); stft(idx) = opts.dB_min;
% end

normalise_time_factor = 1;

% in = transpose(insig_stft);

imagesc(t/normalise_time_factor, f, insig_stft);

opts = Ensure_field(opts,'Colour','Gray');

switch opts.Colour
    case 'Gray'
        cmap = colormap('Gray'); 
        cmap = flipud(cmap); % lighter for lower levels
    case 'RdGy'
        mymap = Get_colourmap_rgb('RdGy');
        cmap = colormap(mymap);
end

colormap(cmap);
 
set(gca,'YDir','normal'); % set(gca,'YDir','reverse');

% set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
NfrequencyTicks = 9;
if ~isfield(opts,'YTickLabel')
    YTickLabel = round(f(round(linspace(1,end,NfrequencyTicks))));
else
    YTickLabel = opts.YTickLabel;
end

if ~isempty(NfrequencyTicks)
    set(gca,'YTick',linspace(f(1),f(end),NfrequencyTicks))
    set(gca,'YTickLabel',YTickLabel)
end

opts = Ensure_field(opts,'XLabel','Time (s)');

xlabel(opts.XLabel)
ylabel(sprintf('Frequency (Hz)'))
% info = Ensure_field(info,'txtTitle', 'Amplitude spectrogram of the signal');
% 
% if exist('stftR','var')
%     info.txtTitle = [info.txtTitle ' (only L-channel plotted)'];
% end
% 
% title(info.txtTitle)

if opts.bColourBar
    opts = Ensure_field(opts,'label','Magnitude (dB)');
    handl = colorbar;
    % set(handl, 'FontName', 'Times New Roman', 'FontSize', 14)
    ylabel(handl, opts.label)
end
 
disp('')
% if isfield(info,'scale_dB')
%     YTick = get(handl,'YTick');
%     idx = find(YTick+info.scale_dB<0);
%     YTick(idx) = [];
%     YTickLabel = YTick + info.scale_dB;
%     set(handl,'YTick',YTick);
%     set(handl,'YTickLabel',YTickLabel);
% end