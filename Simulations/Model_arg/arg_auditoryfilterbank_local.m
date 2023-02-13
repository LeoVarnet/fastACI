function definput=arg_auditoryfilterbank_local(definput)
% ARG_AUDITORYFILTERBANK_LOCAL
%
%   #License: GPL
%   #Author: Peter Soendergaard (2011): Initial version
%   #Author: Alejandro Osses (2020-2021): Extensions
%   #Author: Piotr Majdak (2021): Adapted to AMT 1.0

%% General
definput.keyvals.dboffset = dbspl(1); % dB Full scale convention
definput.flags.outerear = {'no_outerear','outerear'};
definput.flags.middleear= {'no_middleear','middleear','jepsen2008'};

%% Auditory filterbank
definput.keyvals.flow=80;
definput.keyvals.fhigh=8000;
definput.keyvals.basef=[];
definput.keyvals.bwmul=1;
definput.keyvals.betamul=[]; % if empty, betamul will be re-assessed in gammatone.m

definput.keyvals.fs_up  = [];
definput.flags.internalnoise= {'no_internalnoise', 'internalnoise'};

%% Groups
definput.groups.afb_dau1997 = {'dboffset',100,'basef',1000};
definput.groups.afb_bruce2018_debug = {'no_outerear','middleear','bwmul',[], ...
         'flow',125,'fhigh',16000,'fs_up',100e3};
definput.groups.afb_zilany2014     = {'dboffset',94,'fs_up',100e3,'flow',125};
definput.groups.afb_zilany2014_erb = {'fs_up',100e3,'flow',125,'fhigh',8000, ...
    'bwmul',.6};
definput.groups.afb_zilany2014_aabba = {'fs_up',100e3,'flow',125,'fhigh',8000, ...
    'bwmul',.8};

% relanoiborra2019: basef=8000 Hz, to match this frequency
%                   bwmul=0.5 because their filter design uses 60 bands
%                      spaced at 0.5 ERBN between 100 and 8000 Hz
definput.groups.drnl_relanoiborra2019 = {'outerear','middleear','basef',8000, ...
    'hearing_profile','NH','internalnoise','bwmul',0.5};

definput.groups.afb_osses2021 = {'outerear','middleear','basef',[]}; 
definput.groups.afb_osses2017 = {'outerear','middleear','flow',168,'fhigh',1840,'basef',[]};
definput.groups.afb_vandorpschuitman2013 = {'outerear','middleear','flow',168,'fhigh',1840,'basef',168};
