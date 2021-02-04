function [ ] = plot_channels( X, Y, chan2plot, plotfun, msubplot, nsubplot, chantitles )
%PLOT_CHANNELS Conveniently plots a certain number of channels in a matrix
% of filterbank responses. To be used as plot_channels( X, Y, chan2plot)
%
% plot_channels( X, Y, chan2plot, plotfun, msubplot, nsubplot, chanlabels) uses the
% plotting function in plotfun (default: plot) and creates a
% msubplot-by-nsubplot array of plots with titles in chantitles
%
% Leo Varnet 2016 - last modified 2019

Nplot = length(chan2plot);
if nargin<4
    plotfun = @plot;
end
if nargin<6 | (isempty(msubplot) & isempty(nsubplot))
    nsubplot=ceil(sqrt(Nplot));
    if nsubplot*(nsubplot-2)>=Nplot
        msubplot=nsubplot-2;
    elseif nsubplot*(nsubplot-1)>=Nplot
        msubplot=nsubplot-1;
    else
        msubplot=nsubplot;
    end
end
if nargin<7
    chantitles = [];
end
if nargin==5
    error('You must specify both msubplot and nsubplot: plot_channels( X, Y, chan2plot, plotfun, msubplot, nsubplot )');
end

%figure;

X = X(:);
if ndims(Y)==2
    maxval = max(max(abs(Y(:,chan2plot))));
elseif ndims(Y)==3
    maxval = max(max(max(abs(Y(:,chan2plot,:)))));
end

isubplot = 1;
for ichan=fliplr(chan2plot)
    subplot(nsubplot,msubplot,isubplot);
    if ndims(Y)==2
        plotfun(X,Y(:,ichan));
        maxval = max(max(abs(Y(:,chan2plot))));
    elseif ndims(Y)==3
        plotfun(X,squeeze(Y(:,ichan,:)));
    end
    hold on
    if isempty(chantitles)
        title(['Channel #' num2str((ichan))]);
    else
        title(chantitles(ichan,:));
    end
    
    xlim([X(1) X(end)]); 
    if maxval>0
        ylim([-maxval maxval]);
    end
    isubplot = isubplot+1;
end
xlabel('time (s)');

end

