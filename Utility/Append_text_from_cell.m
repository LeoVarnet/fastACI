function text2show = Append_text_from_cell(text_in_cell, pre_message, during_message, post_message)
% function text2show = Append_text_from_cell(text_in_cell, pre_message, during_message, post_message)
%
% The text in 'pre_message' will be contained in the first element of 
%   text2show, if text_in_cell has X elements, they will be collated
%   as separate cell elements to the text 'during_message', and as last
%   element, 'post_message' will be appended.
%
% Author: Alejandro Osses
% Date: 30/09/2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    pre_message = 'to play the stim again';
end
if nargin < 3
    during_message = 'to play a '; 
end
if nargin < 4
    post_message = 'to leave the warm-up phase';
end

if ~isoctave
    % If MATLAB is used, then we can use the embedded function append.m
    text2show = [pre_message append(during_message, text_in_cell) post_message];
else
    % If GNU Octave is used:
    text2show = [];
    text2show{end+1} = pre_message;
    for i = 1:length(text_in_cell)
        text2show{end+1} = [during_message text_in_cell{i}];
    end
    text2show{end+1} = post_message;
end