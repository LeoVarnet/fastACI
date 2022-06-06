function outs = affichage_tf( data, type, varargin )
% function AFFICHAGE_TF( data, type, t, f)
% 
% 1. Description:
%       It plots the information available in 'data'. The data will be 
%       displayed using the format specified by 'type'.
%       t, f indices de temps et frequence de la matrice
%
% Input arguments:
%   - data
%   - type: ('CI', 'CInorm','pow', 'prob', 'zscore', 'zscorealpha', 'zscore7', 'formants', 'tvalue')
% varargin (Optional input parameters):
%   - 'time_affiche' et 'freq_affiche' intervalles de temps (s) et frequence (Hz) � afficher
%   - 'caxis' [min max] specifies the minimum and maximum value that will be
%         used in the colourbar.
%   - 'title' le(s) nom(s) de la (des) CI(s)
%   - 'suptitle' titre general pour le subplot
%   - 'colorbar' 'yes' pour afficher la colorbar
%   - 'NameResponse' le nom des reponses (reponse positive en deuxieme); 
%          ex: {'signal#1', 'signal#2'}
%   - 'alphamap' pour sp�cifier une map de transparence
%   - 'binaryalphamap' pour sp�cifier une map de transparence binaire
%   - 'alphalevel' pour le niveau de probabilite represente en orange
%   - 'formant' associe a un tableau de cellules {t_F, F} ou {t_F, F, style_F}
%         affiche les formants indiqu�s dans la matrice F aux temps t_F,
%         eventuellement avec les styles contenus dans style_F.
%   - 'formantload' nom d'un fichier .mat � charger et afficher automatiquement
%         (le fichier doit contenir des donn�es de formants structurees comme
%         indique precedemment)
%   - 'convertform' passe les formants de l'�chelle lineaire a l'�chelle
%         logarithmique, avec d�calage de 0.05s
%   - 'ROI' pour afficher un ensemble de ROIs rectangulaires [tminROI1 tmaxROI1
%         fminROI1 fmaxROI1] definies dans un tableau de cellules
%   - 'numROI' pour afficher le numero de la ROI au centre de celle-ci
%   - 'meanf' pour afficher la moyenne (en valeur absolue) des poids en 
%         frequence (seulement pour CI unique) 
%   - 'meant' pour afficher la moyenne (en valeur absolue) des poids en 
%         temps (seulement pour CI unique) 
%   - 'NlinesSubplot' nombre de lignes pour l'affichage subplot
%   - 'NfrequencyTicks' nombre de ticks en fr�quence
%   - 'NcolorTicks' nombre de ticks sur la colorbar
%   - 'FormantWidth' largeur de trait des formants
%   - 'displayaxes' pour afficher les axes ('on' par defaut)
%
% It is possible to call affichage_tf( data, type, 'cfg', cfg, varargin)
% to automatically load the adequate data according to the struct 'cfg'.
%
% The input data can aksi be a CI-matrix (dimension 3), auquel cas la
% fonction affiche des subplots, eventuellement ordonnes par le parametre
% 'order' 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

outs = [];

if iscell(data)
    N_CI=length(data);
    datatemp = zeros(size(data{1},1),size(data{1},2),N_CI);
    for i = 1:N_CI
        data_temp(:,:,i) = data{i};
    end
    data=data_temp;
    clear data_temp
end

time_affiche = [];
freq_affiche = [];
title_CI = [];
suptitle_CI = [];
c_axis = [];
NameResponse={};
alphamap =[];
binaryalphamap =[];
alphalevel = 0.05;
grid_on = 0;
F = [];
t_F = [];
ROI = [];
numROI = 0;
order = 1:size(data,3);
convertform='no';
meant = 'no';
meanf='no';
NlinesSubplot = [];
NfrequencyTicks = 9;
NtimeTicks = [];
NcolorTicks = [];
displaycolorbar = 'yes';
FormantWidth = 2;
displayaxes = 'on';

while ~isempty(varargin)
    if ~isstr(varargin{1})
        t=varargin{1};
        f=varargin{2};
    else
        if strcmp(varargin{1}, 'cfg')
            cfg = varargin{2};
            if isfield(cfg,'time_analysis')
                time_affiche = cfg.time_analysis;
            end
            if isfield(cfg,'freq_analysis')
                freq_affiche = cfg.freq_analysis;
            end
            if isfield(cfg,'time_affiche')
                time_affiche = cfg.time_affiche;
            end
            if isfield(cfg,'freq_affiche')
                freq_affiche = cfg.freq_affiche;
            end
            if isfield(cfg,'NameCond')
                title_CI = ['CI for condition ' cfg.NameCond];
            end
            if isfield(cfg,'NameResponse')
                NameResponse = cfg.NameResponse;
            end
            if isfield(cfg,'ROI')
                if ~isempty(cfg.ROI) && isnumeric(cfg.ROI)
                    ROI={cfg.ROI};
                else
                    ROI = cfg.ROI;
                end
            end
            if isfield(cfg,'t')
                t = cfg.t;
            end
            if isfield(cfg,'f')
                f = cfg.f;
            end
            grid_on = 1;
        elseif strcmp(varargin{1}, 'time_affiche')
            time_affiche = varargin{2};
        elseif strcmp(varargin{1}, 'freq_affiche')
            freq_affiche = varargin{2};
        elseif strcmp(varargin{1}, 'title')
            title_CI = varargin{2};
        elseif strcmp(varargin{1}, 'suptitle')
            suptitle_CI = varargin{2};        
        elseif strcmp(varargin{1}, 'colorbar')
            displaycolorbar = varargin{2};
        elseif strcmp(varargin{1}, 'caxis')
            c_axis = varargin{2};
        elseif strcmp(varargin{1}, 'NameResponse')
            NameResponse = varargin{2};
        elseif strcmp(varargin{1}, 'alphamap')
            alphamap = varargin{2};
        elseif strcmp(varargin{1}, 'binaryalphamap')
            binaryalphamap = varargin{2};
        elseif strcmp(varargin{1}, 'NfrequencyTicks')
            NfrequencyTicks = varargin{2};        
        elseif strcmp(varargin{1}, 'NtimeTicks')
            NtimeTicks = varargin{2};        
        elseif strcmp(varargin{1}, 'NcolorTicks')
            NcolorTicks = varargin{2};
        elseif strcmp(varargin{1}, 'alphalevel')
            alphalevel = varargin{2};
        elseif strcmp(varargin{1}, 'ROI')
            ROI = varargin{2};
            if ~isempty(ROI) && isnumeric(ROI)
                ROI={ROI};
            end
        elseif strcmp(varargin{1}, 'numROI')
            numROI = isyes(varargin{2});
        elseif strcmp(varargin{1}, 'grid')
            grid_on = isyes(varargin{2});
        elseif strcmp(varargin{1}, 'formant')
            var = varargin{2};
            if length(var)<=2
                var{3}={};
            end
            if isempty(F)
                F = {var{2}};
                t_F = {var{1}};
                style_F ={var{3}};
            else
                F = {F{:} var{2}};
                t_F = {t_F{:} var{1}};
                style_F = {style_F{:} var{3}};
            end
        elseif strcmp(varargin{1}, 'order')
            order = varargin{2};
        elseif strcmp(varargin{1}, 'convertform')
            convertform = varargin{2};
        elseif strcmp(varargin{1}, 'formantload')
            if exist(varargin{2})
                load(varargin{2});
            end
        elseif strcmp(varargin{1}, 'meant')
            meant = varargin{2};
        elseif strcmp(varargin{1}, 'meanf')
            meanf = varargin{2};
        elseif strcmp(varargin{1}, 'NlinesSubplot')
            NlinesSubplot = varargin{2};
        elseif strcmp(varargin{1}, 'FormantWidth')
            FormantWidth = varargin{2};
        elseif strcmp(varargin{1}, 'displayaxes')
            displayaxes = varargin{2};
        else
            error(['error : unknown argument ' varargin{1}]);
        end
    end
    varargin=varargin(3:end);
end

if strcmp(type,'tvalue')
    [~,~,~,stats]=ttest(data,[],[],[],3);
    data=stats.tstat;
end

if ~exist('t','var')
    t=1:size(data,2);
end
if ~exist('f','var')
    f=1:size(data,1);
end
if isempty(time_affiche)
    time_affiche = [t(1) t(end)];
end
if isempty(freq_affiche)
    freq_affiche = [f(1) f(end)];
end
time_affiche_Index=find(t>=time_affiche(1) & t<=time_affiche(2));
freq_affiche_Index=find(f>=freq_affiche(1) & f<=freq_affiche(2));

% calcule le nombre de subplot pour une image multiCI
N_CI = size(data,3);
if ~isempty(NlinesSubplot)
    msubplot=NlinesSubplot;
    nsubplot=ceil(N_CI/msubplot);
else
    nsubplot=ceil(sqrt(N_CI));
    if nsubplot*(nsubplot-2)>=N_CI
        msubplot=nsubplot-2;
    elseif nsubplot*(nsubplot-1)>=N_CI
        msubplot=nsubplot-1;
    else
        msubplot=nsubplot;
    end
end

for num_CI = 1:N_CI
    i_CI=order(num_CI);
    if N_CI==1
        if isyes(meant) || isyes(meanf)
            axes('Position',[0.2,0.23,0.75,0.7])
        end
    else
        subplot(msubplot,nsubplot,num_CI)
    end
    try
        maxabsdata = max(max(max(abs(data(freq_affiche_Index,time_affiche_Index,:)))));
    catch
        disp('')
    end
    datatoplot=data(freq_affiche_Index,time_affiche_Index,i_CI);
    
    if ~isempty(binaryalphamap)
        if iscell(binaryalphamap)
            if length(binaryalphamap)==1
                datatoplot=datatoplot.*binaryalphamap{1};
            else
                datatoplot=datatoplot.*binaryalphamap{num_CI};
            end
        else
            if ndims(binaryalphamap)==2
                datatoplot=datatoplot.*binaryalphamap;
            else
                datatoplot=datatoplot.*binaryalphamap(:,:,num_CI);
            end
        end
    end
    if strcmp(type, 'CInorm')
        imagesc(t(time_affiche_Index),f(freq_affiche_Index),datatoplot/max(abs((datatoplot(:)))));
    elseif length(type)>5 && strcmp(type(1:6), 'zscore')
        imagesc(t(time_affiche_Index),f(freq_affiche_Index),(datatoplot-mean(datatoplot(:)))/std(datatoplot(:)));
    else
        imagesc(t(time_affiche_Index),f(freq_affiche_Index),datatoplot);
    end
    if strcmp(type, 'pow')
        Audition = [0,0,0;0.00603318260982633,0.0129713425412774,0.0422322787344456;0.0120663652196527,0.0259426850825548,0.0844645574688911;0.0180995482951403,0.0389140285551548,0.126696839928627;0.0241327304393053,0.0518853701651096,0.168929114937782;0.0301659144461155,0.0648567155003548,0.211161404848099;0.0361990965902805,0.0778280571103096,0.253393679857254;0.0422322787344456,0.0907993987202644,0.295625954866409;0.0482654608786106,0.103770740330219,0.337858229875565;0.0542986430227757,0.116742081940174,0.380090504884720;0.0603318288922310,0.129713431000710,0.422322809696198;0.0663650110363960,0.142684772610664,0.464555084705353;0.0723981931805611,0.155656114220619,0.506787359714508;0.0784313753247261,0.168627455830574,0.549019634723663;0.121895432472229,0.154575169086456,0.553268015384674;0.165359482169151,0.140522882342339,0.557516336441040;0.208823531866074,0.126470595598221,0.561764717102051;0.252287596464157,0.112418301403522,0.566013097763062;0.295751631259918,0.0983660146594048,0.570261478424072;0.339215695858002,0.0843137279152870,0.574509859085083;0.382679760456085,0.0702614411711693,0.578758180141449;0.426143795251846,0.0562091507017612,0.583006560802460;0.469607859849930,0.0421568639576435,0.587254941463471;0.513071894645691,0.0281045753508806,0.591503262519836;0.556535959243774,0.0140522876754403,0.595751643180847;0.600000023841858,0,0.600000023841858;0.622222244739533,0,0.566666662693024;0.644444465637207,0,0.533333361148834;0.666666686534882,0,0.500000000000000;0.688888907432556,0,0.466666698455811;0.711111128330231,0,0.433333337306976;0.733333349227905,0,0.400000005960465;0.755555570125580,0,0.366666674613953;0.777777791023254,0,0.333333343267441;0.800000011920929,0,0.300000011920929;0.822222232818604,0,0.266666680574417;0.844444453716278,0,0.233333349227905;0.866666674613953,0,0.200000002980232;0.888888895511627,0,0.166666671633720;0.911111116409302,0,0.133333340287209;0.933333337306976,0,0.100000001490116;0.955555558204651,0,0.0666666701436043;0.977777779102325,0,0.0333333350718021;1,0,0;1,0.0500000007450581,0;1,0.100000001490116,0;1,0.150000005960464,0;1,0.200000002980232,0;1,0.250000000000000,0;1,0.300000011920929,0;1,0.349999994039536,0;1,0.400000005960465,0;1,0.449999988079071,0;1,0.500000000000000,0;1,0.550000011920929,0;1,0.600000023841858,0;1,0.649999976158142,0;1,0.699999988079071,0;1,0.750000000000000,0;1,0.800000011920929,0;1,0.850000023841858,0;1,0.899999976158142,0;1,0.949999988079071,0;1,1,0;];
        colormap(Audition);
    elseif strcmp(type(1:2), 'CI') || strcmp(type,'tvalue')
        colormap(jet);
        try
            caxis([-maxabsdata maxabsdata]);
        catch
            warning('maxansdata seems to be a NaN, setting to -1 and 1')
            caxis([-1 1]);
        end
    elseif strcmp(type, 'ACI')
        % New by Alejandro
        mymap = Get_colourmap_rgb('RdGy');
        cmap = colormap(mymap);
        caxis([-maxabsdata maxabsdata])    
    elseif strcmp(type, 'formants')
        colormap([1 1 1])
    elseif length(type)>5 && strcmp(type(1:6), 'zscore')
        colormap(jet);c=caxis; caxis([-max(abs(c)) max(abs(c))])
    elseif strcmp(type, 'prob')
        map=[[zeros(1,66) 0:1/34:1];[zeros(1,33) 0:1/33:1 ones(1,34)];[0:1/33:1 ones(1,67)]]';
        if ~isempty(alphalevel)
            map(1:alphalevel*100,:)=ones(1,alphalevel*100)'*[1 0.5 0];
        end
        colormap(map)
    else
        error(['error : unknown display option : ' type]);
    end
    
    if ~isempty(c_axis)
        caxis(c_axis)
    elseif strcmp(type, 'CInorm')
        caxis([-1 1]);
    elseif strcmp(type, 'prob')
        caxis([0 1]);
    elseif strcmp(type, 'zscore7')
        caxis([-7 7]);
    elseif strcmp(type, 'zscore5')
        caxis([-5 5]);
    elseif strcmp(type, 'zscore3')
        caxis([-3 3]);
    end
    c_axis=caxis;
    
    % colorbar
    if isyes(displaycolorbar)
        if N_CI==1
            tcolorbar=colorbar;
            if ~isempty(NameResponse)
                title(tcolorbar,NameResponse{1});
                xlabel(tcolorbar,NameResponse{2});
            end
            if ~isempty(NcolorTicks)
                set(tcolorbar,'YTick',linspace(c_axis(1),c_axis(2),NcolorTicks))
            end
        end
        if N_CI==1 && length(type)>5 && strcmp(type(1:6), 'zscore')
            ylabel(tcolorbar, 'Z-score')
        elseif N_CI==1 && length(type)>5 && strcmp(type(1:6), 'tvalue')
            ylabel(tcolorbar, 't-value')
        elseif strcmp(type, 'prob')
            ylabel(tcolorbar, 'p-value')
            set(tcolorbar,'YTick',[alphalevel 1])
        end
        outs.tcolorbar = tcolorbar;
        outs.tcolorbar_description = 'Colour bar handle';
    end
   
    axis xy;

    if ~isempty(alphamap)
        if iscell(alphamap)
            if length(alphamap)==1
                alpha(alphamap{1});
            else
                alpha(alphamap{num_CI});
            end
        else
            if ndims(alphamap)==2
                alpha(alphamap);
            else
                alpha(alphamap(:,:,num_CI));
            end
        end
    elseif strcmp(type, 'zscorealpha')
        zscore_thres=abs(norminv(alphalevel/2));
        alpha(0.4+0.6*double(abs((datatoplot-mean(datatoplot(:)))/std(datatoplot(:)))>=zscore_thres));
    end
    
    view(0,90);
    
    lastline=(num_CI>(msubplot*nsubplot)-nsubplot);
    firstcolumn=(mod(num_CI-1,nsubplot)==0);
    
    % afficher les axes (if requested)

    if isyes(displayaxes) && lastline
        xlabel('Time (s)');
        if ~isempty(NtimeTicks)
            XT = linspace(round(10*t(1))/10,round(10*t(end))/10,NtimeTicks);
            set(gca,'XTick',XT);
            
            outs.XTick = XT;
        end
    else
        set(gca,'XTick',[])
    end
    if isyes(displayaxes) && firstcolumn
        ylabel('Frequency (Hz)');
        if ~isempty(NfrequencyTicks)
            N_f = length(f);
            YT  = linspace(f(1),f(N_f),NfrequencyTicks);
            YTL = round(f(round(linspace(1,N_f,NfrequencyTicks))));
            set(gca,'YTick',YT);
            set(gca,'YTickLabel',YTL);
            
            outs.YTick = YT;
            outs.YTickLabel = YTL;
        end
    else
        set(gca,'YTick',[])
    end
    set(gca,'TickDir','out');
    if ~isyes(displayaxes)
        axis off
    end
    
    if ~isempty(title_CI)
        if iscell(title_CI)
            if length(title_CI)==N_CI
                title(strrep(title_CI{i_CI},'_',' '));
            else
                title(strrep(title_CI,'_',' '));
            end
        elseif size(title_CI,1)>1 && size(title_CI,2)>1
            title(strrep(title_CI(i_CI,:),'_',' '));
        else
            title(strrep(title_CI,'_',' '));
        end
    end
    
    % afficher une grille (if requested)
    if grid_on
        grid on
    end
    
    
    % afficher les trac�s des formants (if requested)
    Ftemp=F;
    t_Ftemp=t_F;
        minf=min(f);
        maxf=max(f);
    for i=1:length(F)
        hold on
        % convertir le trac� dans l'�chelle appropri�e (if requested)
        Ftemp{i}(find(Ftemp{i}<minf))=NaN;
        Ftemp{i}(find(Ftemp{i}>maxf))=NaN;
        
        if isyes(convertform)
            param = polyfit(1:length(f),f,10);
            fprecis=polyval(param,1:(1/10):length(f));
            fbis=linspace(f(1),f(end),10*length(f));
            Ftemp{i}(find(~isnan(Ftemp{i})))=arrayfun(@(x) fbis(nearest(fprecis,x)),Ftemp{i}(find(~isnan(Ftemp{i}))));
            t_Ftemp{i}=t_Ftemp{i}+0.05;
        end
        style=style_F{i};
        plot(t_Ftemp{i},Ftemp{i},'--','LineWidth',FormantWidth,style{:})
        clear style
    end
    
    % afficher une ROI (if requested)
    for i=1:length(ROI)
        hold on
        if ROI{i}(1)>ROI{i}(2)
            temp=ROI{i}(2);ROI{i}(2)=ROI{i}(1);ROI{i}(1)=temp;
        end
        if ROI{i}(3)>ROI{i}(4)
            temp=ROI{i}(4);ROI{i}(4)=ROI{i}(3);ROI{i}(3)=temp;
        end
        rectangle('Position',[ROI{i}(1) ROI{i}(3) ROI{i}(2)-ROI{i}(1) ROI{i}(4)-ROI{i}(3)], 'LineWidth', 2);
        if numROI
            f_text = ROI{i}(3)+(ROI{i}(4)-ROI{i}(3))/2;
            t_text = ROI{i}(1)+(ROI{i}(2)-ROI{i}(1))/2;
            text(t_text,f_text,num2str(i),'HorizontalAlignment','center','Color',[1 1 1], 'FontWeight', 'bold');
        end
    end
end

% afficher un super-titre (if requested)
if ~isempty(suptitle_CI)
    mtit(strrep(suptitle_CI,'_',' '),'xoff', 0, 'yoff',.03,'fontsize',14,'FontWeight', 'bold')
end

% calcul et affichage de la moyenne temporelle et/ou fr�quentielle de l'image (if requested)
if N_CI==1
    if isyes(meant)
        xlabel('')
        set(gca,'XTickLabel',' ')
    end
    if isyes(meanf)
        ylabel('')
        set(gca,'YTickLabel',' ')
        posit = get(gca,'position');
    end
    if isyes(meant)
        CI_f = sum((data(freq_affiche_Index,time_affiche_Index)).^2,2);
        axes('Position',[0.1,posit(2),posit(1)-0.1-0.05,posit(4)])
        plot(1:length(freq_affiche_Index),CI_f);view([90 90])
        xlabel('Frequency (Hz)');
        
        set(gca,'XTick',linspace(1, length(freq_affiche_Index),9))
        set(gca,'XTickLabel',round(f(round(linspace(1,end,9)))))
        set(gca,'XDir','reverse')
        xlim([1 length(freq_affiche_Index)])
    end
    if isyes(meanf)
        CI_t = sum((data(freq_affiche_Index,time_affiche_Index)).^2,1);
        axes('Position',[posit(1),0.1,posit(3),posit(2)-0.1-0.05])
        plot(t(time_affiche_Index),CI_t)
        xlabel('Time (s)');
        xlim([t(time_affiche_Index(1)) t(time_affiche_Index(end))])
    end
end

drawnow
end