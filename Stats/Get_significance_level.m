function [sig_lev, sig_lev2] = Get_significance_level(num_trials, num_choices)
% function [sig_lev, sig_lev2] = Get_significance_level(num_trials, num_choices)
%
%   Gets the value sig_lev. If the scores are above this value then they
%   are significantly above Choice Level.
%   To get this value, it is assumed a t-Student distribution with a
%   p-value less than 0.05. The constant 1.65 is an approximate value of
%   T(t0, alpha = 0.95, N = inf)
%
%   num_choices is used to determine the mean related to the chance level.
%   If the experiment was 2I-2AFC, then num_choices should be 2
%
%%% Example:
%   num_trials  = 50;
%   num_choices = 2;
%   [sig_lev, sig_lev2] = Get_significance_level(num_trials, num_choices);
%
%   num_trials  = 5000;
%   num_choices = 2;
%   [sig_lev, sig_lev2] = Get_significance_level(num_trials, num_choices);
%
%   num_trials  = 4000;
%   num_choices = 2;
%   [sig_lev, sig_lev2] = Get_significance_level(num_trials, num_choices);
%
% Programmed by Matthias Milczynski, commented by Alejandro, ExpORL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
	return;
end

Mean1   = 1/num_choices; % Chance level
n       = num_trials;
S1      =      Mean1 *sqrt(1-Mean1); % Standard deviation considered for sig_lev
S2      = sqrt(Mean1)*sqrt(1-Mean1);

sig_lev = 1.65*S1/sqrt(n) + Mean1;
sig_lev2= 1.65*S2/sqrt(n) + Mean1;

sig_lev = 100*sig_lev;
sig_lev2= 100*sig_lev2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end