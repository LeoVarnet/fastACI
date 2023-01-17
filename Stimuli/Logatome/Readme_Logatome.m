if length(y) <= 2 && length(y) ~= 0
    if ~exist('bShow_header','var');
        bShow_header = 1;
    end

    if bShow_header
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        fprintf('%s.m:\n',upper(mfilename))
        fprintf('\tSpeech sounds from the Logatome corpus.\n')
        fprintf('\tThe files here have been edited by Leo Varnet\n')
        fprintf('\t(e.g., level adjusted, level balanced across phonemes).\n')
    end

    if bShow_header
        fprintf('\tThe following files are used:\n')
    end
    
    try
        for i = 1:length(y)
            fname = strsplit(y{i},'_');
            Speaker_ID = fname{1};
            Sample = fname{2};
            Sample2 = fname{2}([2 1]);

            extra_text = '';
            switch Speaker_ID
                case 'S41F'
                    fprintf('\t\tab_ba, ad_da from Female French speaker S41F (sounds L007_V1_M1_N1, L001_V1_M1_N1)\n') % As-received-20210526, email
                case 'S43M'
                    switch Sample
                        case {'ab','ad'}
                            extra_text = ' (as used in Osses2022b)';
                    end
                    fprintf('\t\t%s_%s from  Male  French speaker S43M%s\n',Sample,Sample2,extra_text);
                case 'S46M'
                    fprintf('\t\tap_pa, at_ta from  Male  French speaker S46M (sounds L008_V6_M1_N2, L002_V6_M1_N2)\n') % As-received-20210409, Slack
                otherwise
                    % Nothing to do
            end
        end
    end
    
end
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
