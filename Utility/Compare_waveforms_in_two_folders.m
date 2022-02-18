function [bSame_sounds, diffe] = Compare_waveforms_in_two_folders(dir1,dir2)
% function [load_name_full,load_name] = Compare_waveforms_in_two_folders(dir1,dir2)
%
% Programmed by Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

files1 = Get_filenames(dir1,'*.wav');
files2 = Get_filenames(dir2,'*.wav');

N = [length(files1) length(files2)];
[xx,idx_ref]  = min(N);
[xx,idx_test] = max(N);

switch idx_ref
    case 1
        files = files1;
        dir_ref = dir1;
        dir_test = dir2;
    case 2
        files = files2;
        dir_ref = dir2;
        dir_test = dir1;
end

for i = 1:length(files)
    
    insig_ref  = audioread([dir_ref  files{i}]);
    insig_test = audioread([dir_test files{i}]);
    
    if length(insig_ref) ~= length(insig_test)
        error('The sounds in dir1 and dir2 don''t seem to be comparable');
    end
    diffe(i) = sum(sum(insig_ref - insig_test)); % second sum, in case the sounds are stereo
    
end

diffe_total = sum(diffe);
if diffe_total == 0
    bSame_sounds = 1;
else
    bSame_sounds = 0;
end
