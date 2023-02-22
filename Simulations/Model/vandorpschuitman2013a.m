function [output, info, Psi_dir, Psi_rev] = vandorpschuitman2013a(insig, fs, varargin)
%VANDORPSCHUITMAN2013A   Binaural auditory model
%
%   Usage: outsig = vandorpschuitman2013a(insig,fs);
%          [outsig, info] = vandorpschuitman2013a(insig,fs);
%          [outsig, info, Psi_dir, Psi_rev] = vandorpschuitman2013a(insig,fs);
%          [outsig, info, Psi_dir, Psi_rev] = vandorpschuitman2013a(insig,fs,...);
%
%   Input parameters:
%        insig  : input (acoustic) signal
%        fs     : sampling frequency [Hz].
%  
%   Output parameters:
%       output  : struct containing the model outputs. Some of the returned
%                 fields include outputs from (1) the binaural central 
%                 processor, (2) raw outputs of the binaural front-end, and
%                 (3) Characterisation of the input signals. More specifically:
% 1- output.prev: Estimate of the perceived reverberation in Model Units (MU)
%  - output.prev_frame: Estimate of the perceived reverberation using a 
%      frame-based analysis, as suggested by Osses et al. (2017, 2020).
%  - output.prev_minmax: Minimum and maximum estimates of perceived reverberation.
%  - output.pcla: Estimate of the perceived clarity in MU       
%  - output.pcla_frame: Estimate of the perceived (frame-based) clarity
%  - output.pcla_minmax: Minimum and maximum estimates of perceived clarity
%  - output.frames_num: Number of analysis frames. frames_num is 1 for the 
%      analysis as suggested by van Dorp et al. (2013).
%  - output.frames_idx: Index of the analysed (non-silent) frames, between 
%      1 and frames_num.
% 2- output.outsig: Binaural internal representation scaled in Model Units.
%  - output.subfs: Sampling frequency (Hz) of the internal representation
%  - output.fc: Centre frequency of the Gammatone filter bands
%  - output.mfc: Centre modulation frequency of the low-pass modulation 
%      filter, which is equal to 8 Hz
% 3- output.insig_level: Calibrated level of the input signal, as assumed by
%      the model. It is related to the keyvalue 'dboffset'
%  - output.insig_level_description: Indicates the details of the values 
%      given by insig_level, dB(A) and dB(Z), using equivalent levels and 
%      max values.
%
%       info    : Extra information obtained from the binaural model.
%       Psi_dir : Internal representation containing only the foreground 
%                 segregated amplitudes 
%       Psi_rev : Internal representation containing only the background
%                 segregated amplitudes 
%
%   The van Dorp et al. (2013) model consists of the following stages:
%   
%   1) Outer- and middle- ear filtering: band-pass filter with pass-band 
%      between 1000 and 4000 Hz. The slopes are approximately 6 dB/Oct.
%
%   2) a gammatone filter bank with 1-ERB spaced filters.
%
%   3-4) an envelope extraction stage done by half-wave rectification
%      followed by low-pass filtering to 770 Hz.
%
%   5) Absolute threshold of hearing: in this implementation this was omitted
%
%   6) an adaptation stage modelling nerve adaptation by a cascade of 5 loops.
%
%   7) Signal smoothing by applying a single-pole low-pass filter with a time
%      constant of 20-ms (fcutoff = 8 Hz; as in dau1996, breebaart2001). 
%
%   8) Central processor
%
%   Any of the optional parameters for |auditoryfilterbank|,
%   |ihcenvelope| and |adaptloop| may be optionally specified for this
%   function. They will be passed to the corresponding functions.
% 
%   See also: vandorpschuitman2013_outmiddlefilter, auditoryfilterbank, 
%             ihcenvelope, adaptloop
%
%   References: breebaart2001binaural, osses2017a

%   #Author: Alejandro Osses 

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

definput.import = {'auditoryfilterbank_local','ihcenvelope','adaptloop','vandorpschuitman2013a'};
definput.importdefaults={'afb_vandorpschuitman2013','ihc_breebaart2001','adt_breebaart2001'};
definput.keyvals.subfs=[];

[flags,kv,flow,fhigh,basef]  = ltfatarghelper({'flow','fhigh','basef'},definput,varargin);

% ------ do the computation -------------------------
insig = gaindb(insig,kv.dboffset-100); % from here on, the input signal is
                               % assumed to be at a dboffset of 100 dB (default AMT)
                               
[insig, kv] = local_signal_preproc(insig,fs,flags,kv);
bDiotic = kv.bDiotic;

if flags.do_outerear || flags.do_middleear
    % This is a combined outer- and middle-ear module. If one of these modules
    %   is requested, then the filter is applied.
    
    [b,a] = vandorpschuitman2013a_outmiddlefilter(fs,kv.order);
    N_cascade = max(size(a,1), size(b,1)); % two cascaded filters
    
    outsig = filter(b(1,:),a(1,:),insig);
    for n = 2:N_cascade
        outsig = filter(b(n,:),a(n,:),outsig);
    end    
else 
    % if both outer and middle ear are by-passed:
    outsig = insig;
end

% Apply the auditory filterbank
[outsig, fc] = auditoryfilterbank(outsig,fs,'argimport',flags,kv);

% 'haircell' envelope extraction
outsig = ihcenvelope(outsig,fs,'argimport',flags,kv);

% non-linear adaptation loops
outsig = adaptloop(outsig,fs,'argimport',flags,kv); % default is with no limitation

% 8-Hz Low-pass filter
mfc = 8; % 'modulation filter frequency Hz'
[mlp_b, mlp_a] = local_IRIfolp(mfc,fs); 
outsig = filter(mlp_b,mlp_a,outsig);

% Downsampling, if resquested:
%   Saves more computer resources for the central processor:
if ~isempty(kv.subfs)
    outsig = fftresample(outsig,round(length(outsig)/fs*kv.subfs));
    fs = kv.subfs;

    N = size(outsig,1);

    if flags.do_exclude
        kv.N            = N;
        kv.N_no_silence = round(kv.ratioFrames*N); 
    else
        kv.N            = N;
        kv.N_no_silence = N;
    end
end;

% Central processor: We set psi following RAA nomenclature
par = [];

par.bDiotic    = bDiotic;
par.fc         = transpose(fc);
par.subfs      = fs;
par.subfs_description = 'Sampling frequency [Hz] of the internal representation in ''outsig''';
par.t_psi      = transpose(( 1:size(outsig,1) )/fs);

if bDiotic == 0
    par.PsiL       = outsig(:,:,1);
    par.PsiR       = outsig(:,:,2);
else
    par.PsiL       = outsig(:,:);
    par.PsiR       = outsig(:,:);
end
par.numBands   = size(par.PsiL,2);
par.numSamples = kv.N_no_silence; % it excludes silent segments if found

switch nargout
    case 3
        [par, Psi_dir] = vandorpschuitman2013a_centralprocessor(par,fs,flags,kv);
    case 4 
        [par, Psi_dir, Psi_rev] = vandorpschuitman2013a_centralprocessor(par,fs,flags,kv);
    otherwise
        par = vandorpschuitman2013a_centralprocessor(par,fs,flags,kv);
end

if isfield(kv,'Level')
    par.Level = kv.Level;
    par.Level_description = kv.Level_description;
end

%%% Outputs from the central processor
output.prev   = par.pRev; par = rmfield(par,'pRev');
output.prev_frame = par.pRev_frame; par = rmfield(par,'pRev_frame');
if isfield(par,'pRev_minmax')
    output.prev_minmax = par.pRev_minmax; par = rmfield(par,'pRev_minmax');
end
output.prev_units = 'MU';
output.prev_description = par.pRev_description; par = rmfield(par,'pRev_description');

output.pcla   = par.pClar; par = rmfield(par,'pClar');
output.pcla_frame = par.pClar_frame; par = rmfield(par,'pClar_frame');
if isfield(par,'pClar_minmax')
    output.pcla_minmax = par.pClar_minmax; par = rmfield(par,'pClar_minmax');
end
output.pcla_units = 'MU';
output.pcla_description = par.pClar_description; par = rmfield(par,'pClar_description');

output.frames_num = par.numFrames;   par = rmfield(par,'numFrames');
output.frames_idx = par.IdxFrames;   par = rmfield(par,'IdxFrames');
if output.frames_num == 1
    output.frames_description = 'full-length insig used to derive prev and pcla';
else
    output.frames_description = 'frame-based analysis of insig to derive prev and pcla as suggested by Osses et al. (2017, 2020)';
end

%%% Output of the front-end (after the 8-Hz low-pass modulation filter)
output.outsig = outsig;
% sampling frequency of the internal representations:
output.subfs = par.subfs;   par = rmfield(par,'subfs');
output.fc    = par.fc;      par = rmfield(par,'fc');
output.mfc   = mfc;

if isfield(par,'numBands')
    par = rmfield(par,'numBands');
end
if isfield(par,'numSamples')
    par = rmfield(par,'numSamples');
end

par = rmfield(par,'t_psi'); % time stamp of the internal representation

%%% Characterisation of the input signal
if isfield(par,'Level')
    output.insig_level = par.Level;
    output.insig_level_type = par.Level_description;
    output.insig_level_description = 'Calibrated levels using an embedded sound level meter';
    
    par = rmfield(par,'Level');
    par = rmfield(par,'Level_description');
end

%%% Needed
% subfs: 11025
% numFrames: 5
% FL: 41.9609
% BL: 9.1292
% TL: 48.0539
% TL_energy: [1×16 double]
% BL_energy: [1×16 double]
% FL_energy: [1×16 double]
% FL_frame: [55.8145 53.2458 50.7480 41.4852 31.1182]
% BL_frame: [11.3044 11.4419 11.1895 9.6359 7.6795]
% TL_frame: [63.2984 60.7995 58.2109 47.9988 36.4062]
% t_frame: [0 121550625 243101250 364651875 486202500]
% Psimin: 0.3400
% Psimin_dip: 0.0604
% LpsiL: [1×16 double]
% LpsiR: [1×16 double]
% IdxFrames: [1 2 3 4 5]
% pRev: 9.1292
% pRev_frame: [11.3044 11.4419 11.1895 9.6359 7.6795]
% pClar: 4.5964
% pClar_frame: [4.9374 4.6536 4.5353 4.3053 4.0521]
% pRev_osses2017: [11.1895 7.6795 11.4419 5]
% pRev_osses2017_description: 'Prev total derived from the frame based analysis (median, min, max, Nr of frames)'
% pClar_osses2017: [4.5353 4.0521 4.9374 5]
% pClar_osses2017_description: 'Pcla total derived from the frame based analysis  (median, min, max, Nr of frames)'
% Level: [66.1294 73.8390 65.8973 74.2191]
% Level_description: {'LAeq'  'LAmax'  'LZeq'  'LZmax'}

par = rmfield(par,'Psimin');     % centralproc configuration: already in keyvals
par = rmfield(par,'Psimin_dip'); % centralproc configuration: already in keyvals

info.extra_outputs = par;
info.flags = flags;
info.keyvals = kv;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [insig, kv] = local_signal_preproc(insig,fs,flags,kv);
% Truncates the signal to multiples of one second:
%    This makes the global estimates equivalent to the mean of the framed-based
%    values (if framelen is 5 s and hopsize 1 s, then the mean of the estimates
%    1, 6, 11... until the end (in steps of framelen/hopsize) provide this
% equivalence.

if size(insig,2) == 2
    % if the input signal is stereo (normal config) then is NOT diotic
    bDiotic = 0;
elseif size(insig,2) == 1 
    % if the input signal is mono then the 'Left' channel will be duplicated
    % and used as 'Right' channel. This is known as diotic
    bDiotic = 1;
end
kv.bDiotic = bDiotic;

Nactual = size(insig,1);
% if flags.do_multiple_hopsize
%     error('Not re-validated yet')
%     %%% Still to validate (draft implementation in Alejandro's local computer)
% else
N = Nactual; 
% end

% ------ do the computation -------------------------
if isempty(kv.hopsize)
    % i.e., if no hopsize (s) is specified:
    kv.hopsize = size(insig,1)/fs;
end
if isempty(kv.framelen)
    % i.e., if no hopsize (s) is specified:
    kv.framelen = size(insig,1)/fs;
end
if flags.do_exclude
    %%% Still to validate (draft implementation in Alejandro's local computer)
	dBFS = 100; % it is always this convention at this point
    dB = local_SLM(insig,fs,'A','f',dBFS); 
    
    if bDiotic == 0
        dB = 10*log10(0.5*abs(10.^(dB(:,1)/10)+10.^(dB(:,2)/10)));
    end
    Level(1) = local_get_Leq(dB); % LAeq
    Level(2) = max(dB); % LAmax
    Level_description = {'LAeq','LAmax'};
    
    LAeq     = local_get_Leq( dB,fs);
    LAeq1sec = local_get_Leq( dB,fs, kv.hopsize,kv.framelen); % 1-sec LAeq

    %%% Now flat:
    dB = local_SLM(insig,fs,'Z','f',100);
    if bDiotic == 0
        % If signal is stereo, the SPL are averaged into one channel
        dB = 10*log10(0.5*abs(10.^(dB(:,1)/10)+10.^(dB(:,2)/10)));
    end
    Level(3) = local_get_Leq(dB); % LZeq
    Level(4) = max(dB); % LZmax
    Level_description(end+1:end+2) = {'LZeq','LZmax'};
    %%%
    
    idx_no_silence = find( (LAeq1sec) >  LAeq-30 ); % looks for analysis frames
                                                    % with average levels greater than
                                                    % LAeq-30 [dB]
    NFrames = length(LAeq1sec); % NFrames = (siglen - framelen)/hopsize + 1
    NFrames_no_sil = length(idx_no_silence);
    ratioFrames = NFrames_no_sil/NFrames;
    
    kv.N              = N; % multiples of 1-s
    kv.ratioFrames    = ratioFrames;
    kv.N_no_silence   = round(ratioFrames*N); 
    kv.idx_no_silence = idx_no_silence;
    
    kv.Level = Level;
    kv.Level_description = Level_description;
end
if flags.do_no_exclude
    kv.N              = N;
    kv.N_no_silence   = N;
    kv.idx_no_silence = 1:size(insig,1)/(kv.hopsize*fs);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [b,a,timeconstant] = local_IRIfolp(f0,fs)
% function [b,a,timeconstant] = local_IRIfolp(f0,fs)
%
% Borrowed from the Two!Ears repository.

a = exp( -(2*pi*f0)/fs );
b = 1 - a;
a = [1, -a];

timeconstant = 1/(2*pi*f0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [outsig_dB, dBFS] = local_SLM(insig,fs,weight_freq,weight_time,dBFS)
% function [outsig_dB, dBFS] = local_SLM(insig,fs,weight_freq,weight_time,dBFS)
%
% 1. Description:
%       outsig is the weighted pressure in [Pa]
% 
% 2. Stand-alone example:
%   fs = 44100;
%   dur = 1; % s
%   insig = Create_sin(1000,dur,fs);
%   Do_SLM(insig,fs,'A','f');
% 
%   [insig, fs] = Wavread('D:\Databases\dir03-Speech\dutch\LISTman\jwz551.wav');
%   Do_SLM(insig,fs,'Z','f');
% 
% 3. Additional info:
%       Tested cross-platform: No
%       See also: Get_Leq
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Created on    : 12/07/2016
% Last update on: 12/07/2016 
% Last use on   : 19/07/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[b,a] = local_weighting_filters(fs,weight_freq);

% same processing, but without AMT, as done in PsySound:
dBoffset = 0.93; % determined empirically on 13/07/2016 to obtain the same values
                 % with this implementation and the one in PsySound.
cal = 10.^((dBFS+dBoffset-94)/20); % new internal dboffset = 93.07
insig = cal*insig; % same as: insig = setdbspl(insig,rmsdb(insig)+dBFS,'dboffset',94);

outsig = filter(b,a,insig);
outsig = local_time_integrator(outsig,fs,weight_time);

outsig_dB = 20*log10(abs(outsig)/2e-5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [b,a] = local_weighting_filters(fs,weightingType)
% function [b,a] = Gen_weighting_filters(fs,weightingType)
% function [b,a] = il_gen_weighting_Filters(fs,weightingType)
%
% 1. Description:
%		Generates the weighting filters with the follwing poles:
%	f1 =    20.6 Hz -> w1 =   129.43 rad/s % low frequency pole, see e.g. IEC 61672-1:2002, 5.4.6, 5.4.11
%	f2 =   107.7 Hz -> w2 =   676.70 rad/s
% 	f3 =   737.9 Hz -> w3 =  4636.36 rad/s
%	f4 = 12194   Hz -> w4 = 76617.16 rad/s % high frequency pole
%  	f5 =   158.5 Hz -> w5 =   995.88 rad/s (source Osses2010)
% 	A-curve uses: f1(x2),2(x1),3(x1),4(x2)
% 	B-curve uses: f1(x2),5(x1),4(x2)
%	C-curve uses: f1(x2),4(x2)
%
% 2. Stand-alone example:
%       fs = 44100; % Hz
%       Gen_weighting_filters(fs,'A');
% 
%       fs = 44100; % Hz
%       Gen_weighting_filters(fs,'C');
% 
% 3. Additional info:
%       References: Osses2010 (thesis), section 2.2.6; IEC 61672-1:2002
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Created on    : 14/07/2016
% Last update on: 14/07/2016 
% Last use on   : 14/07/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 16384; % arbitrary value
% the frequency of each Fourier bin
f = (0:fs/2/(N-1):fs/2);

% take a perceptual sampling of the frequency domain - for design purposes 
eMin   = freqtoaud(min(f)); % ERB min
eMax   = freqtoaud(max(f)); % ERB max
eScale = eMin:(eMax-eMin)/(length(f)-1):eMax; % ERB scale
fScale = audtofreq(eScale); % frequencies sample according to a linear ERB scale

% fLinear = f; % save the linear frequency scale
f = fScale;  % switch the reference frequencies to be f

s = i*2*pi*f; % set up the s-plane variable

% determine the weighting filter frequency responses convienient to 
% accurately set the desired filter orders (n,m)
  
w1 =   129.43; % rad/s
w2 =   676.70; % rad/s
w3 =  4636.36; % rad/s
w4 = 76617.16; % rad/s
% w5 =   995.88; % rad/s
switch weightingType
    case 'A' % A-weighting filter
	    K = 7.3901e9; 
        % freqResp = K*s.^4./((s+w1).^2 .* (s+w2).*(s+w3).*(s+w4).^2);
    
        n = 4; % at most we need a 4th order filter
        m = n; 
  
        zrs =  [0; 0; 0; 0];
        pls = -[w1; w1; w2; w3; w4; w4];
   
    case 'Z' % unweighted
        a = 1;
        b = 1;
        return
    
    otherwise % unknown request
        error('Unknown weighting curve. weightingType can take the values ''A'' or ''Z''')
end
  
m = m+1;
n = n*2;
m = m*2;
  
% generate the filter
if (2*f ~= fs) % correct small frequency error on the last fourier sample.
    f(end) = fs/2;
end
% Use the bilinear transformation to discretise the above transfer function.

[Zd, Pd, Kd] = bilinear(zrs, pls, K, fs);
[b, a] = zp2tf(Zd, Pd, Kd);
warning_state = warning('off', 'MATLAB:nearlySingularMatrix');
warning(warning_state);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function outsig = local_time_integrator(insig,fs, weightingType)
% function outsig = local_time_integrator(insig,fs, weightingType)
%  
% Generates and applies the temporal integration filter to 'insig'
%     
% tau is a time constant for the implemented leaky integrator. In other words,
% this implementation acts as a low-pass filter. It has the following transfer
% function :
%                      1
%        H(s) =  ---------------
%                tau s  +   1

% Filter coeffecients
switch weightingType
    case 'f' % fast leak - time constant = 125 ms
        tau = 125e-3;
    case 's' % slow leak - time constant = 1 s
        tau = 1;
    case 'i'
        tau = 35e-3; % impulse
    otherwise
        error(['integrator: unknown leak case ' char(fastOrSlow)]);
end

% Exponential term
E = exp(-1/(tau*fs));

b = 1 - E;  % Filter numerator - with gain adjustment
a = [1 -E]; % Filter denominator

% State vector
Z = [];

% Create run function handle
% Use filter to perform the integration
outsig = filter(b, a, abs(insig), Z, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Leq = local_get_Leq(levels,fs,dt,framelen_s)
% function Leq = local_get_Leq(levels,fs,dt,framelen_s)
% 
% 1. Description:
%       Computes the equivalent continuous sound level Leq of a sound
% 
% 2. Stand-alone example:
%       % Obtain a one-second Leq:
%       lvls = Do_SLM(insig,fs,'A','f',100); % A-weighted, 'fast'
%       dt = 1; %s
%       Leq = Get_Leq(lvls,fs,dt); % Make sure you enter only mono signals
% 
% 3. Additional info:
%       Tested cross-platform: Yes
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2017
% Created on    : 13/07/2016
% Last update on: 03/02/2017 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin > 2
    if nargin < 4
        framelen_s = dt;
    end
    
    H = round(dt*fs);
    N = round(framelen_s*fs);
    
    if size(levels,2) ~= 1
        error('Calculation validated only for mono inputs');
    end
    Ov = N-H; % Overlap
    levels = buffer(levels,N,Ov,'nodelay');
end

for i = 1:size(levels,2)
    idx   = find(~isnan(levels(:,i)));
    lvls  = levels(idx,i);
    lvls = (10.^(lvls/10));
    Leq(i)= 10*log10( mean(lvls) );
end
