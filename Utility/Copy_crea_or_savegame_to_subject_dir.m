function bCopied = Copy_crea_or_savegame_to_subject_dir(crea_or_savegame,fname_src,type)
% function bCopied = Copy_crea_or_savegame_to_subject_dir(crea_or_savegame,fname_src,type)
%
% Copies from the stored cfg_crea or savegame files in the fastACI repository
%   to the participants' local folder.
%
% Programmed by Alejandro Osses, ENS, France 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dir_src = [fileparts(fname_src) filesep];
file_src = fname_src(length(dir_src)+1:end);
    
dir_res = Get_fastACI_subject_dir(crea_or_savegame.experiment_full,crea_or_savegame.Subject_ID);
dir_res = [dir_res 'Results' filesep]; 
if ~exist(dir_res,'dir')
    mkdir(dir_res);
end

idxf = strfind(file_src,crea_or_savegame.experiment)-1;
idxi = strfind(file_src(1:idxf-1),'_'); idxi = idxi(end);
file_dst = [dir_res file_src(1:idxi) crea_or_savegame.Subject_ID file_src(idxf:end)];

bCopied = 0;
if ~exist(file_dst,'file')
    
    switch type
        case 'cfg_crea'
            cfg_crea = crea_or_savegame;
            save(file_dst,'cfg_crea');
            
        case 'cfg_game'
            cfg_game = crea_or_savegame;
            load([dir_src file_src],'data_passation');
            save(file_dst,'cfg_game','data_passation');
    end
    
    % copyfile(fname_src,file_dst);
    bCopied = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end