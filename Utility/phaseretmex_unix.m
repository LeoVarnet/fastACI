function phaseretmex_unix(dir_phaseret)
%PHASERETMEX_UNIX(dir_phaseret)
%
%AUTHOR: Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 0
    if ~exist('fastACI_dir_phaseret.m','file')
        fastACI_set_phaseret;
    end
    dir_phaseret = fastACI_dir_phaseret;
end

if isunix
    currdir = pwd;
    
    cd([dir_phaseret 'libltfat' filesep]);
    [status,res] = system('make','-echo'); % makes sure the libraries are compiled
    
    cd([dir_phaseret 'libphaseret' filesep]);
    [status,res] = system('make','-echo'); % makes sure the libraries are compiled
    
    % We now compile the mex file:
    dir_where = [dir_phaseret 'mex' filesep];
    
    cd(dir_where)

    % -v is verbose (to print on screen)
    % -I is to include a folder or library
    % ostools.mk defines some useful operating-system dependent commands
    % if the following command gives an error install the library libfftw3: sudo apt-get install libfftw3-dev
    mex -v comp_rtpghiupdate.c -I../libltfat/include/ -I../libphaseret/include/ -I../libltfat/thirdparty/ -I../libphaseret/ostools.mk  ../libphaseret/build/libphaseretd.a ../libltfat/build/libltfatd.a -lfftw3
end

cd(currdir)

%%% Some trial an error experiences during the compilation:
%.$(EXT):  Makefile_unix config.h
%	$(CC) $(CFLAGS) $(SEARCHPATHS) -o $@ -lc -lm $(MEXLINKFLAGS) $(TARGETUP) $(FFTWLIB)
% 
% mex -v comp_rtpghiupdate.c
% 
%%% The following was the output of running the mex file, but with errors:
% /usr/bin/gcc -c -DMATLAB_DEFAULT_RELEASE=R2017b  -DUSE_MEX_CMD   -D_GNU_SOURCE -DMATLAB_MEX_FILE  -I"/usr/local/MATLAB/R2020b/extern/include" -I"/usr/local/MATLAB/R2020b/simulink/include" -fexceptions -fPIC -fno-omit-frame-pointer -pthread -O2 -fwrapv -DNDEBUG "/home/alejandro/Documents/Databases/Toolbox/phaseret-0.2.1/mex/comp_rtpghiupdate.c" -o /tmp/mex_45882243441094_13353/comp_rtpghiupdate.o
% 
% $(CC) $(CFLAGS) $(MEXCOMPFLAGS) $< -o $@ -Wl,--dll -L./  $(MEXLINKFLAGS) $(TARGETUP) $(FFTWLIB) -static-libgcc 
% gcc   -shared -s -Wall -std=c99 -I../libphaseret/include -I../libltfat/include -I../libltfat/thirdparty   -L"$(MATLABROOT)\bin\$(ARCH)" -lmex -lmx    ../libphaseret/% build/libphaseretd.a ../libltfat/build/libltfatd.a                -static-libgcc
