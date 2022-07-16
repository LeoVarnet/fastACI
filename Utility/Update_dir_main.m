function cfg_inout = Update_dir_main(dir_new,cfg_inout)
% function Check_if_dir(dir2check,number_lvl_up)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isfield(cfg_inout,'dir_main')
    cfg_inout.dir_main_old = cfg_inout.dir_main;
end
cfg_inout.dir_main = dir_new;

%%%
if isfield(cfg_inout,'dir_speech')
    if strcmp(cfg_inout.dir_speech(end),'/') || strcmp(cfg_inout.dir_speech(end),'\')
        dir_check = fileparts(fileparts(fileparts( cfg_inout.dir_speech )));
        dir_check = [dir_check filesep];
    else
        %%% Nothing to do, but an error will raise
    end
    if strcmp(dir_check,dir_main_old)
        L = length(dir_check);
        cfg_inout.dir_speech = [dir_new cfg_inout.dir_speech(L+1:end-1) filesep];
        disp('dir_speech updated...')
    end
end

%%%
if isfield(cfg_inout,'dir_noise')
    if strcmp(cfg_inout.dir_noise(end),'/') || strcmp(cfg_inout.dir_noise(end),'\')
        dir_check = fileparts(fileparts(fileparts( cfg_inout.dir_noise )));
        dir_check = [dir_check filesep];
    else
        %%% Nothing to do, but an error will raise
    end
    if ~isfield(cfg_inout,'dir_main_old') || strcmp(dir_check,cfg_inout.dir_main_old)
        L = length(dir_check);
        cfg_inout.dir_noise = [dir_new cfg_inout.dir_noise(L+1:end-1) filesep];
        disp('dir_noise updated...')
    end
    
end

disp('')