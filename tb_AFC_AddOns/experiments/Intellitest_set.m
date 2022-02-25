function Intellitest_set
% function Intellitest_set - setup function of experiment 'exampleVlMatrix'
%
% See also help exampleVlMatrix_cfg, exampleVlMatrix_user, afc_main

global def
global work
global setup
global msg

global global_vars

% make condition dependend entries in structure set

Conditions = work.condition;

Cond_str = strsplit(Conditions,'_');

Cond1 = Cond_str{1};
if length(Cond_str) > 1
    Cond2 = Cond_str{2};
else
    Cond2 = '';
end

switch Cond1
case {'varnoise16'} % 'condition' is so far unused
   work.sentences_nr = 16;
   work.varnoise = 1;
case {'varspeech16'} % 'condition' is so far unused
   work.sentences_nr = 16;
   work.varnoise = 0;
otherwise
   error('Condition not recognised, this should indicate the number of sentences in the lists');
end
work.varspeech = ~work.varnoise;

switch Cond2
    case ''
        % Nothing to do
    otherwise
        idx = strfind(Cond2,'mod');
        if ~isempty(idx)
            fmod = str2double(Cond2(idx+3:end));
            setup.fmod = fmod;
            setup.fmod_unit = 'Hz';
        end
end

%%% Loading some locations:
MaterialPath = intellitest_getpaths('material_location'); 
testListPath = intellitest_getpaths('TestLists'); 

%%% Loading the noise:
work.noise_filename = [MaterialPath 'LTASS_noise_Intellitest.wav'];

[noise,fs_noise] = audioread(work.noise_filename);
if fs_noise ~= def.samplerate
    error('The noise to be used has a different sampling frequency than the speech samples');
end
%%%

if isfield(setup,'fmod')
    start_phase = -pi/2; % Use -pi/2 to begin in minimum
    m = 100;
    [noise,env] = ch_am(noise,setup.fmod,m,'m',fs_noise,start_phase);
end
work.noise = noise;

if work.sentences_nr ~= def.sentences_nr
    error('Condition specified in afc_main mismatches the number of lists specified in cfg.m file (changed manually)...')
end

% select 20 or 30 testlists or demo (go to exampleVlMatrix_cfg)
switch work.sentences_nr
    case 16
        testListFolderFilePart = 'intellitest16.';
        
    otherwise
        error('Not added yet')
end
%testListFolderFilePart = 'demo/demo15.'; % \ changed by / by AO

words_ref = msg.buttonString(:);

for idx = 1:def.interleavenum
    testListFile = [ testListPath testListFolderFilePart num2str(work.int_exppar1{idx}) ];
     
    tmp = strvcat(textread( testListFile ,'%s','headerlines',4));
    count = 1;
    for i = length(tmp):-1:1
        str = tmp(i,:);
        if isempty(strfind(str,':'));
            tmp(i,:) = [];
        else
            idx_loc = strfind(str,':');
            setup.testlistStrings{idx}{count} = tmp(i,1:idx_loc-1);
            count = count+1;
        end
    end
     
    % %%%%%%%%%%% sanity check
    if ( def.interleaved && def.OLSA_SplitTestlistAcrossTracks )
        error('Continue validating here...')
        % if (size(setup.testlistStrings{idx},1) ~= def.maxiter*def.interleavenum )
        %     error('Items in testlist do not match def.maxiter. def.OLSA_SplitTestlistAcrossTracks is set to 1, def.maxiter must equal half the number of items in the testlist.');
        % end
    else
        if (size(setup.testlistStrings{idx},2) ~= def.maxiter )
            error('Items in testlist do not match def.maxiter');
        end
    end
	%%%%%%%%%%%%%%%%%%%%%%
    for i=1:size(setup.testlistStrings{idx},2)
        str = setup.testlistStrings{idx}(i);
        str =str{1}(1:end-2);
        
        for j = 1:length(words_ref)
            if strcmp(str,words_ref{j})
                idx_loc = j;
            end
        end
        work.Intellitest_list{idx}(i) = idx_loc;
        idx_loc = [];
    end
    % work.OLSA_ResponseVector{idx} = zeros(1,5); % intellitest_wordbuttonfcn
end

%%% Cosine ramp to be applied to the background noise
dur_ramp     = 50e-3; 
dur_ramp_samples = def.samplerate*dur_ramp;
win_both     = hannfl(2*dur_ramp_samples,dur_ramp_samples,dur_ramp_samples);
setup.win_up = win_both(1:dur_ramp_samples);
setup.win_dn = win_both(dur_ramp_samples:end);
%%%

if isfield(global_vars,'dBFS')
    dBFS = global_vars.dBFS;
else
    dBFS = 100;
    warning('AO: Calibration assumed to be equal to AMT convention...')
end

work.currentCalLevel = dBFS; 
%soundmexpro('show');
% eof