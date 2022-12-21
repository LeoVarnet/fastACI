function [cmp,max_colour,min_colour] = Get_colourmap_rgb(nam); % ,n,keepAlpha,pyCmd)
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

% 'DG_red_blue' % Double gradient: from red (max) to blue (min)

if nargin == 0
    nam = 'RdGy';
end

cmp = il_get_colour(nam);
if nargout >= 2
    max_colour = cmp(end,:);
end
if nargout >= 3
    min_colour = cmp(1,:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cMap = il_get_colour(colour)

switch colour
    case {'RdGy','rdgy'}
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
    0.9186    0.9186    0.9186; 0.9091    0.9091    0.9091; 0.8995    0.8995    0.8995; 0.8899    0.8899    0.8899; ...
    0.8803    0.8803    0.8803; 0.8690    0.8690    0.8690; 0.8573    0.8573    0.8573; 0.8456    0.8456    0.8456; ...
    0.8338    0.8338    0.8338; 0.8221    0.8221    0.8221; 0.8104    0.8104    0.8104; 0.7986    0.7986    0.7986; ...
    0.7869    0.7869    0.7869; 0.7752    0.7752    0.7752; 0.7634    0.7634    0.7634; 0.7517    0.7517    0.7517; ...
    0.7400    0.7400    0.7400; 0.7278    0.7278    0.7278; 0.7121    0.7121    0.7121; 0.6963    0.6963    0.6963; ...
    0.6806    0.6806    0.6806; 0.6648    0.6648    0.6648; 0.6491    0.6491    0.6491; 0.6333    0.6333    0.6333; ...
    0.6176    0.6176    0.6176; 0.6019    0.6019    0.6019; 0.5861    0.5861    0.5861; 0.5704    0.5704    0.5704; ...
    0.5546    0.5546    0.5546; 0.5389    0.5389    0.5389; 0.5222    0.5222    0.5222; 0.5043    0.5043    0.5043; ...
    0.4864    0.4864    0.4864; 0.4685    0.4685    0.4685; 0.4506    0.4506    0.4506; 0.4327    0.4327    0.4327; ...
    0.4148    0.4148    0.4148; 0.3969    0.3969    0.3969; 0.3790    0.3790    0.3790; 0.3611    0.3611    0.3611; ...
    0.3432    0.3432    0.3432; 0.3252    0.3252    0.3252; 0.3073    0.3073    0.3073; 0.2909    0.2909    0.2909; ...
    0.2752    0.2752    0.2752; 0.2594    0.2594    0.2594; 0.2437    0.2437    0.2437; 0.2279    0.2279    0.2279; ...
    0.2122    0.2122    0.2122; 0.1964    0.1964    0.1964; 0.1807    0.1807    0.1807; 0.1650    0.1650    0.1650; ...
    0.1492    0.1492    0.1492; 0.1335    0.1335    0.1335; 0.1177    0.1177    0.1177; 0.1020    0.1020    0.1020];
    case {'Audition','audition'}
cMap = [ ...
    0 0 0; ...
    0.00603318260982633,0.0129713425412774,0.0422322787344456; ...
    0.0120663652196527, 0.0259426850825548,0.0844645574688911; ...
    0.0180995482951403, 0.0389140285551548,0.126696839928627; ...
    0.0241327304393053, 0.0518853701651096,0.168929114937782; ...
    0.0301659144461155, 0.0648567155003548,0.211161404848099; ...
    0.0361990965902805, 0.0778280571103096,0.253393679857254; ...
    0.0422322787344456, 0.0907993987202644,0.295625954866409; ...
    0.0482654608786106, 0.103770740330219, 0.337858229875565; ...
    0.0542986430227757, 0.116742081940174, 0.380090504884720; ...
    0.0603318288922310, 0.129713431000710, 0.422322809696198; ...
    0.0663650110363960, 0.142684772610664, 0.464555084705353; ...
    0.0723981931805611,0.155656114220619,0.506787359714508;0.0784313753247261,0.168627455830574,0.549019634723663;0.121895432472229,0.154575169086456,0.553268015384674;0.165359482169151,0.140522882342339,0.557516336441040;0.208823531866074,0.126470595598221,0.561764717102051;0.252287596464157,0.112418301403522,0.566013097763062;0.295751631259918,0.0983660146594048,0.570261478424072;0.339215695858002,0.0843137279152870,0.574509859085083;0.382679760456085,0.0702614411711693,0.578758180141449;0.426143795251846,0.0562091507017612,0.583006560802460;0.469607859849930,0.0421568639576435,0.587254941463471;0.513071894645691,0.0281045753508806,0.591503262519836;0.556535959243774,0.0140522876754403,0.595751643180847;0.600000023841858,0,0.600000023841858;0.622222244739533,0,0.566666662693024;0.644444465637207,0,0.533333361148834;0.666666686534882,0,0.500000000000000;0.688888907432556,0,0.466666698455811;0.711111128330231,0,0.433333337306976;0.733333349227905,0,0.400000005960465;0.755555570125580,0,0.366666674613953;0.777777791023254,0,0.333333343267441;0.800000011920929,0,0.300000011920929;0.822222232818604,0,0.266666680574417;0.844444453716278,0,0.233333349227905;0.866666674613953,0,0.200000002980232;0.888888895511627,0,0.166666671633720;0.911111116409302,0,0.133333340287209;0.933333337306976,0,0.100000001490116;0.955555558204651,0,0.0666666701436043;0.977777779102325,0,0.0333333350718021;1,0,0;1,0.0500000007450581,0;1,0.100000001490116,0;1,0.150000005960464,0;1,0.200000002980232,0;1,0.250000000000000,0;1,0.300000011920929,0;1,0.349999994039536,0;1,0.400000005960465,0;1,0.449999988079071,0;1,0.500000000000000,0;1,0.550000011920929,0;1,0.600000023841858,0;1,0.649999976158142,0;1,0.699999988079071,0;1,0.750000000000000,0;1,0.800000011920929,0;1,0.850000023841858,0;1,0.899999976158142,0;1,0.949999988079071,0;1,1,0;];

    case {'DG_red_blue'}
    max_colour_val = 1; % Yves: 1
    N_steps = 128;
%%% As suggested by Yves (Slack, 30/11/2022):    
%     cMap = [ [ linspace(0,max_colour_val,N_steps)' ; ones(N_steps,1)] ,...
%              [ linspace(0,max_colour_val,N_steps)' ; linspace(max_colour_val,0,N_steps)' ] ,...
%              [ ones(100,max_colour_val)            ; linspace(max_colour_val,0,N_steps)' ] ];
    cMap = [ [ linspace(0,1,N_steps)' ; ones(N_steps,1)] ,...
             [ linspace(0,1,N_steps)' ; linspace(1,0,N_steps)' ] ,...
             [ ones(N_steps,1)        ; linspace(1,0,N_steps)' ] ];
    cMap = max_colour_val*cMap;     
    
    case {'DG_jet'}
    max_colour_val = 0.5; 
    N_steps = 128;
    
    cMap = [ [ linspace(0,1,N_steps)' ; linspace(1,max_colour_val,N_steps)'] ,...
             [ linspace(0,1,N_steps)' ; linspace(1,0,N_steps)' ] ,...
             [ linspace(max_colour_val,1,N_steps)' ; linspace(1,0,N_steps)' ] ]; % col 3
    
    case {'SQ_red'}
    max_colour_val = 1;
    N_steps = 128;
    cMap = [ 1+0*linspace(max_colour_val,0,N_steps)' ,...
             linspace(max_colour_val,0,N_steps)' ,...
             linspace(max_colour_val,0,N_steps)' ];

    otherwise
        error('Colour not added yet...')
end

disp('')