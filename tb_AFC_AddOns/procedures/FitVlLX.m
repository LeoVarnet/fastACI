function out = FitVlLX(inLevels, inTrialIntelligibilities, inTargetIntelligibility);
% function out = FitVlLX(inLevels, inTrialIntelligibilities, inTargetIntelligibility);
%
% This is the same function as FitOlLX.m used in the OLSA test.
%
% SE 10.01.2013 14:54
% Matlab conversion of FitOlLX.cpp
% optimized version using vextor operations
%
% FitOlLXm(Levels, TrialIntelligibilities, TargetIntelligibility)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dL50 = [];
dSlope = [];
m_padLikeliHoods = [];
%dResult = 1.0; % SE needed?
dTargetIntelligibility = inTargetIntelligibility; 

for i = 1:length(inLevels)         
	[dL50, dSlope, m_padLikeliHoods] = il_FitOlL50(inLevels(i), inTrialIntelligibilities(i), dL50, dSlope, m_padLikeliHoods);
end
            
out = il_Ol_LevFitLevelAtTargetIntelligbility(dL50, dSlope, dTargetIntelligibility);

%%%%
% End-of-file
%%%

% Inline functions:
% calculates level at particular threshold from psychometric function
function thres = il_Ol_LevFitLevelAtTargetIntelligbility( dL50, dSlope, dTargetIntelligibility )

   if (dTargetIntelligibility == 0.5)
      thres = dL50;
      return;
   end


% This formulas are taken from:
% 
% Brand, T. and Kollmeier, B. (2002),
% "Efficient adaptive procedures for threshold and concurrent slope estimations
% for psychophysics and speech intelligibility tests",
% J. Acoust. Soc. Am. 111(6), 2801-2810
%
% Formula (A). Discrimination function. This corresponds to formula (1)
% in the paper:
% 
% p(L,L_50,s_50) =  1 / ( 1 + exp(4*s_50*(L_50-L)) )
%
% with the following equivalences:
%    dLevel:  L
%    dL50:    L_50
%    dSlope:  s_50
%   double dPSentence = 1.0/(1.0 + exp(4.0*dSlope*(dL50-dLevel)));

% Thus the reverse function is
dLevel = -1000.0;
dLog = 1.0 /  dTargetIntelligibility - 1.0;

if (dLog > 0.0)
      dLevel = dL50 - log ( 1.0 /  dTargetIntelligibility - 1.0 ) / (4.0 * dSlope);
end
   
thres = dLevel;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%------------------------------------------------------------------------------
% \brief calculates the logarithmic likelihood for a result of a sentence in a
% sentence test with the parameters of the discrimination function.
% \param[in] dLevel presentation level
% \param[in] dResult intelligibility
% \param[in] dL50 SRT
% \param[in] dSlope slope at SRT (1/dB)
% \param[in] dPMax maximum intelligibility
% \retval logarithmic likelihood
%------------------------------------------------------------------------------
function logLikeOut = il_LogLikeSentence( dLevel,  dResult,  dL50,  dSlope,  dPmax)
	% SE 10.01.2013 14:48 handles dL50 as vector in Matlab version
   
    FIT_EPS = 10e-9;
    
  if (dPmax ~= 1.0)
      error('PMax == 1.0 expected');
	end
		
		
   % This formulas are taken from:
   % 
   % Brand, T. and Kollmeier, B. (2002),
   % "Efficient adaptive procedures for threshold and concurrent slope estimations
   % for psychophysics and speech intelligibility tests",
   % J. Acoust. Soc. Am. 111(6), 2801-2810
   %
   % Formula (A). Discrimination function. This corresponds to formula (1)
   % in the paper:
   % 
   % p(L,L_50,s_50) =  1 / ( 1 + exp(4*s_50*(L_50-L)) )
   %
   % with the following equivalences:
   %    dLevel:  L
   %    dL50:    L_50
   %    dSlope:  s_50
   
   dPSentence = 1.0./(1.0 + exp(4.0*dSlope*(dL50-dLevel)));

   % Formula (B). Logarithmic likelihood function, corresponding to
   % formula (8) from paper. Here no product is calculated, but only
   % the logarithmic likelihood for a given intelligibility (dResult).
   % The accumulation (building of product) is done below in FitOlL50,
   % where these logarithmic (!) likelihoods are added!
   
   logLikeOut =  dResult * log(dPSentence + FIT_EPS) + (1.0 - dResult) .* log(1.0 - dPSentence + FIT_EPS);
%------------------------------------------------------------------------------


%------------------------------------------------------------------------------
% \brief performs 'Oldenburg' maximum likelihood fit.
% \param[in] reference to header where to read/write measurement data, here
% - [in|out] 'DataInstance' pointer to allocated likelihood data
% - [out] L50 SRT
% - [out] Slope slope
% - [out] PMax maximum intelligibility
% \param[in] bWithSlope flag, if its an 'SRT' or 'SRT and slope' measurement
% NOTE: function allocates data for the likelihoods and stores the pointer to the data
% in the internal string m_strAllocatedDataInstances for cleanup, and in the passed
% header for later reusage.
% NOTE: no exceptions are caught, only an additional error message is printed before
% re-throw!
%------------------------------------------------------------------------------
function [dL50, dSlope, m_padLikeliHoods] = il_FitOlL50( dLevel, dResult, dL50, dSlope, m_padLikeliHoods)
	
	FIT_HUGE_VAL = 10e9;
	
   % retrieve all necessary parameters
   % - Par1 values
   dPar1      = -80.0;
   dPar1Max   = 120.0;
   dPar1Step  = 0.1;

   % - number of steps for all three parameters
   nPar1 = floor((dPar1Max-dPar1 )/dPar1Step)+1;
   nPar2 = 20;     % No slope: only one variable parameter !!
   nPar3 = 1;    % No slope: only one variable parameter !!
   % set initial Likelihood
   dLike = -FIT_HUGE_VAL;
   %dLogLikeSentence, dPMax;


   if isempty(m_padLikeliHoods)
      m_padLikeliHoods  = zeros(1,nPar1*nPar2*nPar3);
		end

      % fill Par2 and Par arrays
     
      % No slope: only one variable parameter (i.e. ony one value, no arrays)
      %for (j = 0; j < nPar2; j++)
      %   padPar2[j] = 0.15; //(double)rhdConfig[STLEVFIT_FITCTRL_OL50_LATTICE][STLEVFIT_FITCTRL_OL50_PAR2_LIST][j].GetValue();
      
      % SE for faster vector calculation
      padPar1 = [dPar1:dPar1Step:dPar1Max];
      
      
      padPar2(1) = 0.02;
      padPar2(2) = 0.03;
      padPar2(3) = 0.04;
      padPar2(4) = 0.05;
      padPar2(5) = 0.06;
      padPar2(6) = 0.07;
      padPar2(7) = 0.08;
      padPar2(8) = 0.09;
      padPar2(9) = 0.1;
      padPar2(10) = 0.12;
      padPar2(11) = 0.14;
      padPar2(12) = 0.16;
      padPar2(13) = 0.18;
      padPar2(14) = 0.2;
      padPar2(15) = 0.22;
      padPar2(16) = 0.24;
      padPar2(17) = 0.28;
      padPar2(18) = 0.32;
      padPar2(19) = 0.4;
      padPar2(20) = 0.5;

      for (k = 1:nPar3)
         padPar3(k) = 1.0; %(double)rhdConfig[STLEVFIT_FITCTRL_OL50_LATTICE][STLEVFIT_FITCTRL_OL50_PAR3_LIST][k].GetValue();
        end

      % now start searching the complete 3d-lattice
      %for (i = 1:nPar1) SE replaced by vector operation 10.01.2013 14:45
         
         for (j = 1:nPar2) % SE if further speedup is required, replace this by vector and pass a matrix (nPar1 x nPar2) to LogLikeSentence
            
            for (k = 1:nPar3)
               
               % SE padPar1 vector in one chunk
               dLogLikeSentence = il_LogLikeSentence(dLevel, dResult, padPar1, padPar2(j), padPar3(k));
               
               iLikeliHoods = (j-1)*nPar3*nPar1 + (k-1)*nPar1 + 1;
               m_padLikeliHoods(iLikeliHoods:iLikeliHoods+nPar1-1) = m_padLikeliHoods(iLikeliHoods:iLikeliHoods+nPar1-1) + dLogLikeSentence;
               
               [maxElement, maxElementIndex] = max(m_padLikeliHoods(iLikeliHoods:iLikeliHoods+nPar1-1));
               
               if (maxElement > dLike)
                  
                  dLike    = maxElement;
                  dL50     = padPar1(maxElementIndex);
                  dSlope   = padPar2(j);
                  dPMax    = padPar3(k);
                end
             end
         end
         
  

   if (dLike == -FIT_HUGE_VAL)
      dL50 = -1000.0;
      dSlope = -1000.0;
   end
%------------------------------------------------------------------------------
