function [outsig,fc,mfc,outsig_afb] = relanoiborra2019_preproc_debug(insig, fs, varargin)
%RELANOIBORRA2019_PREPROC_DEBUG Creates internal representation based on Relano-Iborra et al. (2019)
%   This is a script based on RELANOIBORRA2019_FEATUREEXTRACTION but including
%       improvements in the coding, reusing existing functions from AMT and a
%       reformating of the modulation filter bank outputs so that the outputs
%       are arranged in a cell array (as happens in dau1997.m and osses2021.m).
%       The current implementation from AMT (relanoiborra2019_featureextraction)
%       removes columns of the modulation filter bank in its back-end script 
%       (relanoiborra2019_decision.m: see 'rule4th') something that is conceptually
%       wrong, as the rule4th (see Verhey et al. 1999) is a limitation of the
%       number of modulation filters per cochlear filter that is related to the
%       'preprocessing'.
%
%   Other issues to solve by the AMT team:
%       - It is conceptually wrong to adopt the name 'outsig_afb' (output 
%         auditory filter bank, because here it is used for 'everything that
%         comes before the modulation filter bank'
%       - It is wrong to have the 50 dB gain after the DRNL, because this gain
%         should amplify the output of the DRNL with the hearing-threshold 
%         internal noise. Note that in Helia's implementation the 50-dB gain 
%         comes before the expansion stage and NOT after the filter bank.
%
%   Usage: [out, varargout] = relanoiborra2019_preproc_debug(insig, fs, varargin)
%          [out, varargout] = relanoiborra2019_preproc_debug(insig, fs, flow, fhigh, varargin)
%
%   Input parameters:
%     insig      :  signal to be processed
%     fs         :  Sampling frequency
%     flow       :  lowest center frequency of auditory filterbank
%     fhigh      :  highest center frequency of auditory filterbank
%     sbj        :  subject profile for drnl definition
%
%   Output parameters:
%     out           : correlation metric structure inlcuding
%
%   The out struct contains the following fields:
%
%     '.dint'      correlation values for each modulation band
%
%     '.dsegments' correlation values from each time window and mod. band.
%
%     '.dfinal'    final (averaged) correlation
%
%   Description:
%     This script builds the internal representations of the template and target signals
%     according to the CASP model (see references).
%     The code is based on previous versions of authors: Torsten Dau, Morten
%     Loeve Jepsen, Boris Kowalesky and Peter L. Soendergaard
%
%     The model has been optimized to work with speech signals, and the
%     preprocesing and variable names follow this principle. The model is
%     also designed to work with broadband signals. In order to avoid undesired
%     onset enhancements in the adaptation loops, the model expects to recive a
%     prepaned signal to initialize them.
%
%   References: jepsen2008 relanoiborra2019

%   REFERENCES:
%
%   Jepsen, M. L., Ewert, S. D., & Dau, T. (2008). A computational model
%   of human auditory signal processing and perception. Journal of the
%   Acoustical Society of America, 124(1), 422-438.
%
%   Relano-Iborra, H., Zaar, J., & Dau, T. (2019). A speech-based computational
%   auditory signal processing and perception model (sCASP). The Journal of the
%   Acoustical Society of America, 146(5), 3306-3317.

%   #Author: Helia Relano Iborra (March 2019): v4.0 provided to the AMT team
%   #Author: Clara Hollomey (2021): adapted to the AMT
%   #Author: Piotr Majdak (2021): adapted to the AMT 1.0
%   #StatusDoc: Good
%   #StatusCode: Good
%   #Verification: Unknown
%   #Requirements: M-Stats M-Signal
%
%   #Author 'improvements': Alejandro Osses (2022)
mfc = [];
fc = [];

definput.import={'auditoryfilterbank','ihcenvelope','adaptloop','modfilterbank', ...
    'relanoiborra2019_preproc_debug'}; % my local extra settings
definput.importdefaults={... %'drnl_relanoiborra2019',
    'ihc_relanoiborra2019','adt_relanoiborra2019','mfb_jepsen2008', ...
    'afb_relanoiborra2019_preproc'}; % settings to be loaded
definput.keyvals.subfs=[];

[flags,kv] = ltfatarghelper({'flow','fhigh'},definput,varargin);

%% Auditory filtering:
if isoctave
   warning(['Currently this model is only fully functional under MATLAB.']); 
end

%% Calculate fc's
if flags.do_erbspace,
  fc = erbspace(kv.flow, kv.fhigh, kv.fcnt); % method from Relano-Iborra et al. (2019)
end
if flags.do_erbspacebw,
  % find the center frequencies used in the filterbank, 1 ERB spacing
  fc = erbspacebw(kv.flow, kv.fhigh, kv.bwmul, kv.basef); % method from Osses et al. (2022)
end
nchannels = size(fc,2);

%% Outer- and middle-ear filtering
if flags.do_outerear
    b_hp  = headphonefilter(fs); % calc headphone filtercoeffs
    insig = filter(b_hp,1,insig); % Outer-ear filterring
end
if flags.do_middleear
    b_me  = middleearfilter(fs);
    insig = filter(b_me,1,insig); % middle-ear-ear filterring
end

%% Cochlear filter bank, the DRNL:
% Pre-allocate memory for auditory bands
outsig_afb = zeros(length(insig),nchannels);

for n = 1:nchannels % Filter
    %%% AMT 1.1:    
    
    %%% AO:
    hearing_profile = kv.subject; % 'NH' or 'HIx'
    % outsig_afb(:,n) = relanoiborra2019_drnl_NH(insig,fc(n),fs,hearing_profile); % using the local drnl function
    outsig_afb(:,n) = relanoiborra2019_drnl(insig',fc(n),fs, hearing_profile); % using the AMT function
    % Small details about the AMT implementation (relanoiborra2019_drnl.m):
    %     - the input is a row array: this is completely against the AMT conventions
end

if flags.do_internalnoise % flags.do_internalnoise_hearingthresholds
    % Internal noise related to absolute threshold of hearing
    N_samples = size(outsig_afb,1);
    int_noise = relanoiborra2019_internalnoise(N_samples, fc); % Alejandro's local function   
    outsig_afb    = outsig_afb  +int_noise;
end

%% 'Haircell' envelope extraction
if flags.do_ihc
    outsig_afb = ihcenvelope(outsig_afb,fs,'argimport',flags,kv);
end

% For use with DRNL, this cannot be turned off:
outsig_afb   = outsig_afb * 10^(50/20); % Gain to compensate for the Outer/middle ear attenuation

%% Expansion (Auditory-nerve inspired) and Non-linear adaptation loops:
if flags.do_adt || flags.do_mfb
    % Expansion - Motivation, from Jepsen et al. 2008: '' However, a squaring 
    %   expansion was introduced in the model after hair-cell transduction, reflecting 
    %   the square-law behavior of rate-versus-level functions of the neural 
    %   response in the AN (Yates et al., 1990; Muller et al., 1991)''.
    outsig_afb = outsig_afb.^2;
    
    % non-linear adaptation loops
    outsig_afb = adaptloop(outsig_afb,fs,'argimport',flags,kv);
end

if flags.do_mfb
    %% Modulation processing:
    % modfilterbank.m is equivalent to 'relanoiborra2019_mfbtd.m' and 'mfbtdpp.m' togheter:
    [outsig,mfc] = modfilterbank(outsig_afb,fs,fc,'argimport',flags,kv);
end