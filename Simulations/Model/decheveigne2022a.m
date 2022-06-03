function [outsig, fc, params] = decheveigne2022a(insig, fs, varargin);
%DECHEVEIGNE2022A   Simple auditory nerve model
%
%   Usage: [outsig, fc] = decheveigne2022a(insig,fs);
%          [outsig, fc] = decheveigne2022a(insig,fs,...);
%          [outsig, fc, mfc] = decheveigne2022a(insig,fs,...);
%          [outsig, fc, mfc, params] = decheveigne2022a(insig,fs,...);
%
%   Input parameters:
%     insig  : input acoustic signal.
%     fs     : sampling rate.
%  
%   `decheveigne2022a(insig,fs)` computes the internal representation of 
%   the signal *insig* sampled with a frequency of *fs* Hz.
%  
%   `[outsig,fc,params]=decheveigne2022a(...)` additionally returns the 
%   some extra parameters obtained from the auditory-nerve synapse model.
%  
%   The model consists of the following stages:
%   
%   1) an outer- and middle-ear filtering as used by Jepsen et al. 2008
%
%   2) a gammatone filter bank with 1-erb spaced filters.
%
%   3) an envelope extraction stage done by half-wave rectification
%      followed by low-pass filtering to 770 Hz as used by Breebaart et al. 2001
%
%   4) an adaptation stage modelling nerve adaptation by a cascade of 5
%      loops using a limiter factor of 5 (Osses and Kohlrausch, 2021).
%
%   5) a modulation filterbank
%
%   Any of the optional parameters for |auditoryfilterbank|,
%   |ihcenvelope| and |adaptloop| may be optionally specified for this
%   function. They will be passed to the corresponding functions.
%
%   See also: auditoryfilterbank, ihcenvelope, exp_osses2021, exp_osses2022,
%             breebaart2001, lopezpoveda2001, dau1997
%
%   References: dau1997mapI dau1997mapII

%   #StatusDoc: Perfect
%   #StatusCode: Perfect
%   #Verification: Unknown
%   #Requirements: M-Signal
%   #Author : Alejandro Osses.

% ------ Checking of input parameters ------------

if nargin<2
  error('%s: Too few input arguments.',upper(mfilename));
end;

if ~isnumeric(insig) 
  error('%s: insig must be numeric.',upper(mfilename));
end;

if ~isnumeric(fs) || ~isscalar(fs) || fs<=0
  error('%s: fs must be a positive scalar.',upper(mfilename));
end;
% load defaults from arg_auditoryfilterbank, arg_ihcenvelope, arg_adaptloop, arg_modfilterbank and arg_osses2021
definput.import={'auditoryfilterbank_local','ihcenvelope_local','adaptloop','modfilterbank_local','osses2021','decheveigne1990'}; 
definput.importdefaults={'afb_osses2021','ihc_breebaart2001', 'adt_osses2021','mfb_osses2022a'}; 
definput.keyvals.subfs=[];

[flags,keyvals]  = ltfatarghelper({'flow','fhigh'},definput,varargin);

fc  = [];
params = [];
% ------ do the computation -------------------------

insig = gaindb(insig,keyvals.dboffset-94); % from here on, the input signal is
                               % assumed to be at a dboffset of 100 dB (default AMT)
M = size(insig,2); 
if M ~= 1
    error('%s: The input signal should be monaural (each column is one signal)',upper(mfilename));
end
if flags.do_outerear
    hp_fir = headphonefilter(fs); % Getting the filter coefficients at fs
    N = floor(length(hp_fir)/2); % group delay for a FIR filter of order length(hp_fir)
    insig = [insig; zeros(N,M)]; % group delay compensation: step 1 of 2. 
    insig = filter(hp_fir,1,insig); % filtering
    insig = insig(N+1:end,1:M); % group delay compensation: step 2 of 2
end

if flags.do_no_outerear
    if keyvals.silent_mode == 0
        %fprintf('\t%s: outer-ear filtering by-passed\n',upper(mfilename));
        amt_disp([upper(mfilename),': outer-ear filtering by-passed']);
        amt_disp();
    end
end 

if flags.do_middleear || flags.do_no_middleear
    filtertype = 'lopezpoveda2001';
elseif flags.do_jepsen2008
    filtertype = 'jepsen2008';
end

if flags.do_middleear || flags.do_jepsen2008
    me_fir     = middleearfilter(fs,filtertype,'zero');
    me_gain_TF = max( 20*log10(abs(freqz(me_fir,1,round(fs/2)))) ); % max of the filter response - 'K=round(fs/2)': arbitrary

    N = floor(length(me_fir)/2); % group delay for a FIR filter of order length(me_fir)
    insig = [insig; zeros(N,M)]; % group delay compensation: step 1 of 2.
    insig = filter(me_fir,1,insig); % filtering
    insig = insig(N+1:end,1:M); % group delay compensation: step 2 of 2. 
    insig = gaindb(insig,-me_gain_TF); % if me_fir is a non-unit gain filter, 
                                       % the gain of the FIR filter is compensated.
    if keyvals.silent_mode == 0
        amt_disp();
        amt_disp([upper(mfilename),': middle-ear filter was adjusted to have a 0 dB bandpass gain (gain applied=',-me_gain_TF,' dB']);
        amt_disp();
    end
end                                    

if flags.do_no_middleear
    if keyvals.silent_mode == 0
        amt_disp();
        amt_disp([upper(mfilename),': middle-ear filtering by-passed']);
        amt_disp();
    end
end

if flags.do_afb
    % Apply the auditory filterbank
    [outsig, fc] = il_auditoryfilterbank(insig,fs,'argimport',flags,keyvals);
end

if flags.do_no_afb
    outsig = insig;
    if keyvals.silent_mode == 0
        %fprintf('\t%s: Gammatone filter bank by-passed\n',upper(mfilename));
        amt_disp();
        amt_disp([upper(mfilename),': Gammatone filter bank by-passed']);
        amt_disp();
    end
end

if flags.do_ihc || flags.do_adt
    if flags.do_ihc_hwr
        outsig = max( outsig, keyvals.ihc_minlvl );
    else
        % 'haircell' envelope extraction
        outsig = ihcenvelope(outsig,fs,'argimport',flags,keyvals);
    end
else
    if keyvals.silent_mode == 0
        %fprintf('\t%s: ihcenvelope processing by-passed\n',upper(mfilename));
        amt_disp();
        amt_disp([upper(mfilename),': ihcenvelope processing by-passed']);
        amt_disp();
    end
end

if flags.do_adt
    % outsig = il_decheveigne1990(outsig,fs,flags,keyvals);
    outsig = il_decheveigne1990_raster(outsig,fs,flags,keyvals);
    disp('')
    
    params.n_fibres = keyvals.n_fibres;
else
    if keyvals.silent_mode == 0
        %fprintf('\t%s: adaptation loops by-passed\n',upper(mfilename));
        amt_disp();
        % amt_disp([upper(mfilename),': adaptation loops by-passed']);
        amt_disp();
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [outsig, fc] = il_auditoryfilterbank(insig, fs, varargin)
%IL_AUDITORYFILTERBANK   Linear auditory filterbank with extra configuration
%   parameters.
%
%   Input parameters:
%     insig  : input acoustic signal.
%     fs     : sampling rate.
%  
% This script is identical to auditoryfilterbank.m, with the exception that
%   'betamul' can be used as an input parameter to allow the use of wider
%   auditory filters.
%
%   AUTHOR : Peter L. Søndergaard
%   Modified by Alejandro Osses
  
% ------ Checking of input parameters ------------

definput.import={'auditoryfilterbank'};
[flags,keyvals,flow,fhigh]  = ltfatarghelper({'flow','fhigh'},definput,varargin);

% ------ do the computation -------------------------

% find the center frequencies used in the filterbank, 1 ERB spacing (if bwmul=1)
fc = erbspacebw(flow, fhigh, keyvals.bwmul, keyvals.basef);

if ~isempty(keyvals.betamul)
    if keyvals.bwmul ~= keyvals.betamul
        fprintf('%s: \tBetamul and bwmul are different, this mught mean that you are\n\t\tunconsciously filters that do not cross at their -3 dB points\n',upper(mfilename))
        fprintf('\t\t(betamul=%.4f; bwmul=%.4f)\n',keyvals.betamul, keyvals.bwmul);
    end
end

% Calculate filter coefficients for the gammatone filter bank.
[gt_b, gt_a]=gammatone(fc, fs, 'complex','betamul',keyvals.betamul);

% Apply the Gammatone filterbank
outsig = 2*real(ufilterbankz(gt_b,gt_a,insig));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function spike_pulses = il_decheveigne1990(insig,fs,flags,keyvals)

gain      = keyvals.gain;
dead_time = keyvals.dead_time;
n_fibres  = keyvals.n_fibres;

nchans = size(insig,2);
dur    = size(insig,1)/fs;
spike_pulses = [];
for iChan = 1:nchans
    s = spike_train(gain*insig(:,iChan),fs,dead_time,n_fibres);
    p = spike_to_pulse(s,fs,dur);
    spike_pulses=[spike_pulses,p];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function spike_pulses = il_decheveigne1990_raster(insig,fs,flags,keyvals)

gain      = keyvals.gain;
dead_time = keyvals.dead_time;
n_fibres  = keyvals.n_fibres;

nchans = size(insig,2);
dur    = size(insig,1)/fs;
spike_pulses = [];
for iChan = 1:nchans
    for i = 1:n_fibres 
        s = spike_train(gain*insig(:,iChan),fs,dead_time,1);
        p = spike_to_pulse(s,fs,dur);
        spike_pulses=[spike_pulses,p];
    end
end