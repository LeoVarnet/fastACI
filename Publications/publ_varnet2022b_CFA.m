function [h,hname] = publ_varnet2022b_CFA(varargin)
% function [h,hname] = publ_varnet2022b_CFA(varargin)
%
% Generates the figures
%
% % To display Fig. 2 of Varnet, Lorenzi, Osses (2022, CFA) use :::
%     publ_varnet2022b_CFA('fig2');
%
% Author: Alejandro Osses and Leo Varnet
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    help publ_varnet2022b_CFA;
    return
end

h = [];
hname = [];

definput.flags.type={'missingflag','fig2_mod22'};
definput.flags.subject = {'SA','SB'};
definput.flags.publ = {'varnet2022b_CFA','varnet2022a_JASA'}; % default is varnet2022b_CFA
definput.keyvals.models=[];

[flags,keyvals]  = ltfatarghelper({},definput,varargin);

dir_data = [fastACI_paths('dir_data') 'modulationACI' filesep];

if flags.do_fig2_mod22
    if flags.do_SA
        files = {'S4'}; % S_LV2
    end
    if flags.do_SB
        files = {'S1'}; % S_AO
    end
end

N_subjects = length(files);

data.N_subjects = N_subjects;

for i_subject = N_subjects:-1:1 % First participant read at last...
    dir_local = [dir_data files{i_subject} filesep];
    
    if flags.do_fig2_mod22
        publ_varnet2022a_utils(files{i_subject},'Get_CIt',flags);
    end
end

%%% End of the script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
