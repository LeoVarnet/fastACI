function output = relanoiborra2019_backend2Dst_ENS_version(int_rep_speech, int_rep_mix, fs, cf_mod)
% function output = relanoiborra2019_backend2Dst_ENS_version(int_rep_speech, int_rep_mix, fs, cf_mod)
%
% Back-end function for the model of Relano-Iborra et al. (2019).
% This script was modified to accept internal representations stored as
% cell arrays (as it is the AMT standard), where each cell is one auditory 
% filter that contains the time signals (Nsamp samples long) of up to 
% 12 modulation filters, already limited to have modulation centre frequencies
% such that: mfc < 1/4 fc.
%
% Authors: Helia Relano-Iborra and DTU colleagues. Modifications by Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nsamp = size(int_rep_mix{1},1);
Naud_ch = length(int_rep_mix);

WinDurs = 1./cf_mod; % The window duration is the inverse of the centerfrequency of the modulation channel
WinDurs(1) = 1/2.5;

WinLengths = floor(WinDurs * fs);
Nsegments = floor(Nsamp./WinLengths)+ ones(1,length(cf_mod)); % The total number of segments is Nframes plus any additional "leftover"
    
if find(WinLengths == Nsamp)% If the duration of the stimulus is exactly equal to the window duration
        segIdx = find(WinLengths == Nsamp);
        Nsegments(segIdx) =  Nsegments(segIdx)-1;
end

for m=1:length(cf_mod)
    
    %%% Rule 4th applied now in modfilterbank.m: 
    % rule4th=find(cf_aud > 4*cf_mod(m)); % Apply rule of cf_mod < 1/4 cf_aud
    speech = [];
    mix = [];
    if m == 5
        disp('')
    end
    for i = 1:Naud_ch
        N_mfc_here = size(int_rep_speech{i},2);
        if m <= N_mfc_here
            speech(:,end+1) = int_rep_speech{i}(:,m);
            mix(:,end+1)    = int_rep_mix{i}(:,m);
        end
    end
        
    tmp_ssnn = zeros(WinLengths(m), size(mix, 2), Nsegments(m)); % Allocate memory for multi-resolution
    tmp_ss = tmp_ssnn;
        
    segLengths = zeros(1,Nsegments(m)) ;
    
    % Find starting and ending points of the segments:
    
     for i = 1:Nsegments(m) % For each temporal segment of the signal
                               % find the start and end index of the frame
            if i > (Nsegments(m)-1)
                startIdx = 1 + (i-1)*WinLengths(m);
                endIdx = Nsamp;
            else
                startIdx = 1 + (i-1)*WinLengths(m);
                endIdx = startIdx + WinLengths(m)-1;
            end
            
            segment = startIdx:endIdx;
            segLengths(i) = length(segment);
            
         % internal representation of the temporal segments (samplesPerSegment x all bands x  number of segments)
         
         tmp_ss(1:segLengths(i), :, i) = speech(segment,:);
         tmp_ssnn(1:segLengths(i), :, i) = mix(segment,:);
     
         dint_temp(i, m) = corr2(tmp_ss(1:segLengths(i), :, i), tmp_ssnn(1:segLengths(i),:, i));
         dint_temp(dint_temp<0)=0; % Remove negative correlations
  
    end
    dmod(m)=nansum(dint_temp(1:Nsegments(m), m))/Nsegments(m);
end

output.dint=dmod;
output.dsegments=dint_temp;
output.dfinal=mean(dmod);
