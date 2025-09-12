function publ_osses2025_Frontiers_figs(varargin)
% function publ_osses2025_Frontiers_figs(varargin)
%
% 1. Description: Generates the figures
%
% % To display Fig. 3 of Osses, Le Bagousse and Varnet, (2025, Frontiers) use :::
%     publ_osses2025_Frontiers_figs('fig3'); 
%
% % To display Fig. 4 of Osses, Le Bagousse and Varnet, (2025, Frontiers) use :::
%     publ_osses2025_Frontiers_figs('fig4'); 
%
% % To display Fig. 5 of Osses, Le Bagousse and Varnet, (2025, Frontiers) use :::
%     publ_osses2025_Frontiers_figs('fig5'); 
%
% % To display Fig. 6 of Osses, Le Bagousse and Varnet, (2025, Frontiers) use :::
%     publ_osses2025_Frontiers_figs('fig6'); 
%
% % To display Fig. 7 of Osses, Le Bagousse and Varnet, (2025, Frontiers) use :::
%     publ_osses2025_Frontiers_figs('fig7'); 
%
% % To display Fig. 8 of Osses, Le Bagousse and Varnet, (2025, Frontiers) use :::
%     publ_osses2025_Frontiers_figs('fig8'); 
%
% % To display Fig. 9 of Osses, Le Bagousse and Varnet, (2025, Frontiers) use :::
%     publ_osses2025_Frontiers_figs('fig9'); 
%
% Author: Leo Varnet
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% close all, clc

if nargin == 0
    help publ_osses2025_Frontiers_figs;
    return
end

definput.flags.type={'missingflag', ...
    'fig3', ...  
    'fig4', ...
    'fig5', ... 
    'fig6', ...  
    'fig7', ...
    'fig8', ...  
    'fig9'}; 
definput.flags.plot={'plot','no_plot'};
definput.flags.local={'local','zenodo'};

definput.keyvals.dir_zenodo=[];
definput.keyvals.dir_out=[];

[flags,keyvals]  = ltfatarghelper({},definput,varargin);

if flags.do_fig3 
    publ_osses2025_fig3;
end

if flags.do_fig4 
    publ_osses2025_fig4;
end

if flags.do_fig5 
    publ_osses2025_fig5;
end

if flags.do_fig6 
    publ_osses2025_fig6;
end

if flags.do_fig7 
    publ_osses2025_fig789;
end

if flags.do_fig8 
    publ_osses2025_fig789;
end

if flags.do_fig9 
    publ_osses2025_fig789;
end
