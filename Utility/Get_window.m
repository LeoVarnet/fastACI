function [win,wtype] = Get_window(nwtype,N,opts)
% function [win,wtype] = Get_window(nwtype,N,opts)
%
% 1. Description:
%       Generates a window in a column vector with length N. Default value
%       of M is 1. If you have a signal stored in a column vector use M to 
%       get the window immediately repeated in columns.
%
%       nwtype - window type (number arbitrarily assigned)
%       wtype  - window type (window name)
%       win    - array containing window
% 
%       nwtype      wtype                   Tested  Tested on
%       0           'rectangular'           No
%       1           'hanning'               Yes     24/06/2014
%       2           'triangular'            No
%       3           'bartlett'              No
%       4           'hamming'               Yes     03/01/2015
%       5           'blackman'              No
%       6           'blackman-harris'       No
%       7           'gaussian(alpha=2.5)'   No
%       8           'gaussian(alpha=3.5)'   No
%       9           'gaussian(alpha=4.5)'   No
%       10          'flattop'               No
% 
% 2. Stand-alone example:
%       N = 4096;
%       [win,wtype] = Get_window('hanning',N);
%       
%       % Other way of obtaining the same Hanning window:
%       [win,wtype] = Get_window(1,N);
% 
% 3. Additional info:
%   Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on   : 16/06/2014
% Last updateon: 04/12/2015 
% Last used on : 04/12/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    opts = 'periodic'; % periodic for analysis purposes, symmetric for filtering purposes
end

M = 1;

if length(N) ~= 1 % then N is assumed to be a 'y'-vector
    N = length(N); 
end

switch nwtype
    case {0, 'rectangular'}
        win = ones(N,1); % 'rectangular window'
        wtype = 'rectangular';
    case {1, 'hanning'}
        win = hanning(N,opts); % periodic for analysis purposes, symmetric for filtering purposes
        wtype = 'hanning';
    case {2, 'triangular'} 
        win = triang(N);
        wtype = 'triangular';
    case {3, 'bartlett'}
        win = bartlett(N);
        wtype = 'bartlett';
    case {4, 'hamming'}
        win = hamming(N,opts);
        wtype = 'hamming';
    case {5, 'blackman'}
        win = blackman(N);
        wtype = 'blackman';
    case {6, 'blackman-harris'}
        H     = sigwin.blackmanharris(N);
        win   = generate(H);
        wtype = 'blackman-harris';
    case {7, 'gaussian(alpha=2.5)'}
        win   = gausswin(N,2.5);
        wtype = 'gaussian(alpha=2.5)';
    case {8, 'gaussian(alpha=3.5)'}
        win   = gausswin(N,3.5);
        wtype = 'gaussian(alpha=3.5)';
    case {9, 'gaussian(alpha=4.5)'}
        win   = gausswin(N,4.5);
        wtype = 'gaussian(alpha=4.5)';
    case {10, 'flattop'}
        win = flattopwin(N);
        wtype = 'flattop';
    otherwise
        error('Not recognised window')
end

win = repmat(win,1,M);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end