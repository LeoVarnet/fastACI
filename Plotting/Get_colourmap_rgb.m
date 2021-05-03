function cmp = Get_colourmap_rgb(nam); % ,n,keepAlpha,pyCmd)
% cmp = Get_colourmap_rgb(nam [, n, keepAlpha, pyCmd])
%
%
% ::INPUT::
%
% nam:      Colormap name from matplotlib library (case sensitve!). See
%           below. Alternatively, '!GetNames' returns a cellstring of
%           available colormaps.
% n:        Number of Colors; defaults to 128
% keepAlpha: Switch to keep the alpha channel of the colormap (4th colum);
%           defaults to false. If true, a Nx4 matrix is returned in cmp
%           (instead of Nx3).
% pyCmd:    python command; defaults to 'python'
%
% 
% ::OUTPUT::
%
% cmp       A Nx3 (Nx4 if keepAlpha is true) double array of RGB(A) values.
%
%
% Colormap name can be one of the following:
%
%     Accent      gist_earth        Oranges           RdYlBu
%     Accent_r    gist_earth_r      Oranges_r         RdYlBu_r
%     afmhot      gist_gray         OrRd              RdYlGn
%     afmhot_r    gist_gray_r       OrRd_r            RdYlGn_r
%     autumn      gist_heat         Paired            Reds
%     autumn_r    gist_heat_r       Paired_r          Reds_r
%     binary      gist_ncar         Pastel1           seismic
%     binary_r    gist_ncar_r       Pastel1_r         seismic_r
%     Blues       gist_rainbow      Pastel2           Set1
%     Blues_r     gist_rainbow_r    Pastel2_r         Set1_r
%     bone        gist_stern        pink              Set2
%     bone_r      gist_stern_r      pink_r            Set2_r
%     BrBG        gist_yarg         PiYG              Set3
%     BrBG_r      gist_yarg_r       PiYG_r            Set3_r
%     brg         GnBu              PRGn              Spectral
%     brg_r       GnBu_r            PRGn_r            spectral
%     BuGn        gnuplot           prism             spectral_r
%     BuGn_r      gnuplot_r         prism_r           Spectral_r
%     BuPu        gnuplot2          PuBu              spring
%     BuPu_r      gnuplot2_r        PuBu_r            spring_r
%     bwr         gray              PuBuGn            summer
%     bwr_r       gray_r            PuBuGn_r          summer_r
%     CMRmap      Greens            PuOr              terrain
%     CMRmap_r    Greens_r          PuOr_r            terrain_r
%     cool        Greys             PuRd              winter
%     cool_r      Greys_r           PuRd_r            winter_r
%     coolwarm    hot               Purples           Wistia
%     coolwarm_r  hot_r             Purples_r         Wistia_r
%     copper      hsv               rainbow           YlGn
%     copper_r    hsv_r             rainbow_r         YlGn_r
%     cubehelix   jet               RdBu              YlGnBu
%     cubehelix_r jet_r             RdBu_r            YlGnBu_r
%     Dark2       nipy_spectral     RdGy              YlOrBr
%     Dark2_r     nipy_spectral_r   RdGy_r            YlOrBr_r
%     flag        ocean             RdPu              YlOrRd
%     flag_r      ocean_r           RdPu_r            YlOrRd_r
% 
% V 1.3; Konrad Schumacher, 07.2019

if nargin == 0
    nam = 'RdGy';
end
% if strcmpi(nam,'!GetNames')
%     % translate switch to retrieve colormap names into python-switch:
%     nam = 'listCMapNames';
% end
% 
% 
% % defaults:
% if ~exist('n','var') || isempty(n)
%     n = 128;
% end
% if ~exist('keepAlpha','var') || isempty(keepAlpha)
%     keepAlpha = 0;
% end
% if ~exist('pyCmd','var') || isempty(pyCmd)
%     pyCmd = 'python';
% end
% 
% 
% % check if python script is present
% pyScript = which('pyplotCMap2txt.py');
% assert(~isempty(pyScript), 'getPyPlot_cMap:PyScriptNotFound', ...
%     'Could not find python script (%s).','pyplotCMap2txt.py');
% 
% tmpf = tempname;
% 
% % call python script
% comd = sprintf('%s "%s" %s -o "%s" -n %d', pyCmd, pyScript, nam, tmpf, n);
% [s,m] = system(comd);
% 
% % check if system command ran w/o error
% assert(s==0, 'getPyPlot_cMap:SystemCMDFailed', ...
%         'There was an error executing the command\n\t%s\nSystem returned:\n\t%s', ...
%         comd, m);
% 
% 
% if strcmp(nam,'listCMapNames')
%     % cMap names retrieved; read text file
%     fid = fopen(tmpf,'r');
%     cmp = textscan(fid,'%s');
%     fclose(fid);
%     cmp = cmp{1};
%     
% else
%     % load cMap data from text file
%     cmp = load(tmpf,'-ascii');
% 
%     if keepAlpha
%     else % remove 4th column of alpha values
%         cmp = cmp(:,1:3);
%     end
% end
% 
% 
% % delete temp-file
% delete(tmpf);

cmp = il_get_colour(nam);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cMap = il_get_colour(colour)

switch colour
    case 'RdGy'
cMap = [ ...
    0.4039         0    0.1216; 0.4271    0.0074    0.1253; 0.4502    0.0148    0.1290; 0.4734    0.0222    0.1327; ...
    0.4966    0.0296    0.1364; 0.5197    0.0371    0.1401; 0.5429    0.0445    0.1438; 0.5660    0.0519    0.1475; ...
    0.5892    0.0593    0.1512; 0.6124    0.0667    0.1549; 0.6355    0.0741    0.1586; 0.6587    0.0815    0.1623; ...
    0.6818    0.0889    0.1660; 0.7014    0.1008    0.1718; 0.7125    0.1230    0.1823; 0.7236    0.1453    0.1928; ...
    0.7347    0.1675    0.2033; 0.7458    0.1897    0.2138; 0.7570    0.2119    0.2243; 0.7681    0.2342    0.2348; ...
    0.7792    0.2564    0.2453; 0.7903    0.2786    0.2558; 0.8014    0.3009    0.2663; 0.8125    0.3231    0.2768; ...
    0.8237    0.3453    0.2873; 0.8348    0.3676    0.2978; 0.8448    0.3893    0.3118; 0.8540    0.4106    0.3281; ...
    0.8633    0.4319    0.3445; 0.8726    0.4532    0.3609; 0.8818    0.4745    0.3772; 0.8911    0.4958    0.3936; ...
    0.9004    0.5171    0.4100; 0.9096    0.5384    0.4263; 0.9189    0.5597    0.4427; 0.9281    0.5810    0.4591; ...
    0.9374    0.6023    0.4754; 0.9467    0.6236    0.4918; 0.9559    0.6449    0.5082; 0.9594    0.6621    0.5290; ...
    0.9621    0.6787    0.5503; 0.9649    0.6954    0.5716; 0.9677    0.7121    0.5929; 0.9705    0.7288    0.6142; ...
    0.9733    0.7454    0.6355; 0.9760    0.7621    0.6568; 0.9788    0.7788    0.6781; 0.9816    0.7955    0.6994; ...
    0.9844    0.8121    0.7207; 0.9872    0.8288    0.7420; 0.9899    0.8455    0.7633; 0.9923    0.8610    0.7839; ...
    0.9929    0.8722    0.8011; 0.9935    0.8833    0.8184; 0.9941    0.8944    0.8357; 0.9948    0.9055    0.8530; ...
    0.9954    0.9166    0.8703; 0.9960    0.9277    0.8876; 0.9966    0.9389    0.9049; 0.9972    0.9500    0.9222; ...
    0.9978    0.9611    0.9395; 0.9985    0.9722    0.9568; 0.9991    0.9833    0.9741; 0.9997    0.9944    0.9914; ...
    0.9952    0.9952    0.9952; 0.9856    0.9856    0.9856; 0.9761    0.9761    0.9761; 0.9665    0.9665    0.9665; ...
    0.9569    0.9569    0.9569; 0.9474    0.9474    0.9474; 0.9378    0.9378    0.9378; 0.9282    0.9282    0.9282; ...
    0.9186    0.9186    0.9186; 0.9091    0.9091    0.9091; 0.8995    0.8995    0.8995; 0.8899    0.8899    0.8899
    0.8803    0.8803    0.8803
    0.8690    0.8690    0.8690
    0.8573    0.8573    0.8573
    0.8456    0.8456    0.8456
    0.8338    0.8338    0.8338
    0.8221    0.8221    0.8221
    0.8104    0.8104    0.8104
    0.7986    0.7986    0.7986
    0.7869    0.7869    0.7869
    0.7752    0.7752    0.7752
    0.7634    0.7634    0.7634
    0.7517    0.7517    0.7517
    0.7400    0.7400    0.7400
    0.7278    0.7278    0.7278
    0.7121    0.7121    0.7121
    0.6963    0.6963    0.6963
    0.6806    0.6806    0.6806
    0.6648    0.6648    0.6648
    0.6491    0.6491    0.6491
    0.6333    0.6333    0.6333
    0.6176    0.6176    0.6176
    0.6019    0.6019    0.6019
    0.5861    0.5861    0.5861
    0.5704    0.5704    0.5704
    0.5546    0.5546    0.5546
    0.5389    0.5389    0.5389
    0.5222    0.5222    0.5222
    0.5043    0.5043    0.5043
    0.4864    0.4864    0.4864
    0.4685    0.4685    0.4685
    0.4506    0.4506    0.4506
    0.4327    0.4327    0.4327
    0.4148    0.4148    0.4148
    0.3969    0.3969    0.3969
    0.3790    0.3790    0.3790
    0.3611    0.3611    0.3611
    0.3432    0.3432    0.3432
    0.3252    0.3252    0.3252
    0.3073    0.3073    0.3073
    0.2909    0.2909    0.2909
    0.2752    0.2752    0.2752
    0.2594    0.2594    0.2594
    0.2437    0.2437    0.2437
    0.2279    0.2279    0.2279
    0.2122    0.2122    0.2122
    0.1964    0.1964    0.1964
    0.1807    0.1807    0.1807
    0.1650    0.1650    0.1650
    0.1492    0.1492    0.1492
    0.1335    0.1335    0.1335
    0.1177    0.1177    0.1177
    0.1020    0.1020    0.1020];
    otherwise
        error('Colour not added yet...')
end

disp('')