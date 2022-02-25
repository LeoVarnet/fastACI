function intellitest_win(action)
% function intellitest_win(action)
%
% Based on 'olsa_win.m' - custom window function for OLSA procedure, version
% 1.30.0 (last modified 01.03.2013 09:36) taken from the AFC toolbox, 
% programmed by Stephan Ewert
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global def
global work
global msg

% if the window is not enabled just leave here
if ( def.afcwinEnable == 0 )
	return;
end

% check whether window is still existing
% if we attempt to call a dead window with any action other than 'open'
% shut the whole thing down
if ( strcmp(action, 'open') == 0 && isempty(findobj('Tag','afc_win')) )
	work.terminate = 1;
	work.abortall = 1;
	%warning('AFC:afcwin', 'Attempt to access non existing response window. Run terminated');
	%warning('Attempt to access non existing response window. Run terminated');
	return;	
end

[numWordRows, numWordCols] = size(msg.buttonString);

switch action
   case 'open'
       
% define number of words       
[num, numWordAlternatives] = size(msg.buttonString); %  msg.buttonString

h = figure('BusyAction','cancel','KeyPressFcn',...			% create modal dialog window
   ['afc_pressfcn(' num2str(num) ',0)'],'Tag','afc_win',...
   'menubar','none',...
   'Windowstyle', 'modal',...
   'Color',[0.75 0.75 0.75],...
   'Interruptible','off',...
   ... % added secure close
   'CloseRequestFcn','afc_close', ...
   ... % added Truetype
   'defaultUIControlFontname','Arial', ...
   ... % change size once to avoid Matlab 6 modal window refresh bug
   'position',def.winFigurePosition - [0 0 1 1], ...
   'Name',['AFC-measurement (' def.version ')'] );

% change to default size again
set(h,'position',def.winFigurePosition);


% 08.11.2006 09:43
% some position adjustments
yShift = 0;
yShiftButtons = 0; 

switch def.winButtonConfiguration
case 1
	% open start and end button, too
	
	b = uicontrol('Parent',h, ...										% create start button
		'Units','normalized', ...
		'FontUnits','normalized', ...
		'BusyAction','cancel', ...
		'Callback',['afc_pressfcn(' num2str(num) ',''s'')'], ...
		'FontSize',0.4, ...
		'Position',[0.91 0.925 0.07 0.07], ...
		'BackgroundColor',[0.75 0.75 0.75],...
		'String',msg.startButtonString, ...
	   'Tag','afc_startbutton');
	      
	 % 15-04-2005 14:15 check version for >=7 and add keypressfcn to buttons
   if ( work.matlabVersion > 6 ) 
   	set(b,'KeyPressFcn',['afc_pressfcn(' num2str(num) ',0)']);
   end
   
   if ( (def.mouse == 0) )
   		set(b,'Enable','off');
   end
	
	b = uicontrol('Parent',h, ...										% create end button
		'Units','normalized', ...
		'FontUnits','normalized', ...
		'BusyAction','cancel', ...
		'Callback',['afc_pressfcn(' num2str(num) ',''e'')'], ...
		'FontSize',0.4, ...
		'Position',[0.91 0.825 0.07 0.07], ...
		'BackgroundColor',[0.75 0.75 0.75],...
		'String',msg.endButtonString, ...
		'Tag','afc_endbutton');
		
		% 15-04-2005 14:15 check version for >=7 and add keypressfcn to buttons
   if ( work.matlabVersion > 6 ) 
   	set(b,'KeyPressFcn',['afc_pressfcn(' num2str(num) ',0)']);
   end
   
   if ( (def.mouse == 0) )
   		set(b,'Enable','off');
   end
   
end

switch def.windetail
case 0
b = uicontrol('Parent',h, ...										% create message display
	'Units','normalized', ...
	'FontUnits','normalized', ...
   'BackgroundColor',[0.9 0.9 0.9], ...
  	'FontSize',0.2, ...
   'Position',[0.1 0.825 0.8 0.15-yShift], ...
	'Style','text', ...
	'Tag','afc_message');
	
	% some position adjustments
   yShift = 0;
   yShiftButtons = 0; 
case 1
   t1 = uicontrol('Parent',h, ...										% create message display
	'Units','normalized', ...
	'FontUnits','normalized', ...
   'BackgroundColor',[0.9 0.9 0.9], ...
  	'FontSize',0.6, ...
   'Position',[0.1 0.925-yShift 0.35 0.05], ...
   'Style','text', ...
   'String',' ', ...
   'Tag','afc_t1');

	t2 = uicontrol('Parent',h, ...										% create message display
	'Units','normalized', ...
	'FontUnits','normalized', ...
   'BackgroundColor',[0.9 0.9 0.9], ...
  	'FontSize',0.6, ...
   'Position',[0.55 0.925-yShift 0.35 0.05], ...
   'Style','text', ...
   'String',' ', ...
	'Tag','afc_t2');

   b = uicontrol('Parent',h, ...										% create message display
	'Units','normalized', ...
	'FontUnits','normalized', ...
   'BackgroundColor',[0.9 0.9 0.9], ...
  	'FontSize',0.4, ...
   'Position',[0.1 0.825-(yShift-yShiftButtons) 0.8 0.075-yShiftButtons], ...
   'Style','text', ...
  	'Tag','afc_message');
end 

% create word matrix

FontSize = 0.3; 

bwidth = 0.8/(num+(0.5*(num-1)));
bsepar=bwidth/4;

bheight = 0.8/(numWordAlternatives+(0.5*(numWordAlternatives-1)));
          
str = strsplit(mfilename,'_');
str = str{1};

for i = 1:num
    for k = 1:numWordAlternatives
        % Setting the buttons
        kWordAlternative = def.WordAlternativeMapping(k,i);
       
        b = uicontrol(  'Parent',h, ...
                        'Units','normalized', ...
                        'FontUnits','normalized', ...
                        'BusyAction','cancel', ...
                        'Callback',[str '_wordbuttonfcn(' num2str(i) ',' num2str(kWordAlternative) ')'], ...
                        'FontSize',FontSize, ...
                        'Position',[0.1+(i-1)*(1.5*bwidth) 0.1-yShiftButtons+(k-1)*(1.3*bheight) bwidth bheight], ...
                        'BackgroundColor',[0.75 0.75 0.75],...
                        'String',msg.buttonString{kWordAlternative,i}, ...
                        'enable','off', ...
                        'Tag',['afc_button' num2str(i)]);
       
        % 15-04-2005 14:15 check version for >=7 and add keypressfcn to buttons
        if ( work.matlabVersion > 6 ) 
            set(b,'KeyPressFcn',[str '_wordbuttonfcn(' num2str(i) ',' num2str(kWordAlternative) ')']);
        end
    end
    
    if i == num
        % word selection field
        ts = uicontrol( 'Parent',h, ... % create message display
                        'Units','normalized', ...
                        'FontUnits','normalized', ...
                        'BackgroundColor',[0.9 0.9 0.9], ...
                        'FontSize',FontSize, ...
                        'Position',[0.1+(i-1)*(1.5*bwidth) 0.1-yShiftButtons-0.08 bwidth bheight], ...
                        'Style','text', ...
                        'String',' ', ...
                        'Tag',['afc_wordSelection' num2str(i)]);
    end
end


% acceptbutton

  	b = uicontrol ...
   ('Parent',h, ...
   'Units','normalized', ...
	'FontUnits','normalized', ...
	'BusyAction','cancel', ...
   'Callback',[str '_pressfcn(5,1)'], ...   % vlmatrix_pressfcn calls function that communicates with AFC
	'FontSize',0.5, ...
	'Position',[0.91 0.02 0.07 bheight], ...
   'BackgroundColor',[0.75 0.75 0.75],...
	'String',msg.okButtonString, ...
   'Tag','afc_okButton');
   
    % 15-04-2005 14:15 check version for >=7 and add keypressfcn to buttons
   if ( work.matlabVersion > 6 ) 
   	set(b,'KeyPressFcn',[str '_pressfcn(5,0)']);
   end



% refresh(h)

%------------------------ end of action 'open'----------------------------

case 'close'
	h=findobj('Tag','afc_win');
	delete(h);
	%close(h);
   
case 'start_ready'
  hm=findobj('Tag','afc_message');					% handle to message box
	h=findobj('Tag','afc_win');						%
	ht2=findobj('Tag','afc_t2');
	ht1=findobj('Tag','afc_t1');
   
   % Was Andy request
   % 01-08-2005 10:30 SE introduce config var for this screen skip because the start_msg
   % might be used for different msg during the experiment run
   % are we really starting the exp or do we start the next measurement?
   % 08.09.2005 14:34 added def.skipStartMessage
   if ( (work.terminate > 0) & (def.skipStartMessage == 1) )
   	set(h,'UserData',0);
		pause( 0.25);	
   else
   	set(hm,'string',msg.start_msg);
		set(h,'UserData',4);
   end
   
   set(ht1,'string',sprintf(msg.experiment_windetail, work.filename));
   set(ht2,'string',sprintf(msg.measurementsleft_windetail, size(work.control,1) + 1 - work.numrun, size(work.control,1)));
   
case 'start'
   hm=findobj('Tag','afc_message');
   ht2=findobj('Tag','afc_t2');

   set(hm,'string','');
   set(ht2,'string',sprintf(msg.measurement_windetail, work.numrun, size(work.control,1)));
   
case 'finished'
   hm=findobj('Tag','afc_message');
   h=findobj('Tag','afc_win');
	set(hm,'string',msg.finished_msg); % called if experiment is already finished
	set(h,'UserData',-2);              % ready to get only end command

case {'correct','false'}
    hm=findobj('Tag','afc_message');
    
    %SE 04.02.2013 14:42 clear buttons
    for clearb = 1:numWordCols
        hb=findobj('Tag',['afc_button' num2str(clearb)]);
        set(hb,'enable','off');
				
        hb=findobj('Tag',['afc_wordSelection' num2str(clearb)]);
        set(hb,'string','');				
    end
    drawnow; % flush event que 15-04-2005 13:38

    if def.feedback
       hm=findobj('Tag','afc_message');
       switch action
           case 'correct'
               set(hm,'string',msg.correct_msg);
           case 'false'
               set(hm,'string',msg.false_msg);
       end
       drawnow; 
    end
    
case 'clear'
    hm=findobj('Tag','afc_message');
    set(hm,'string','');

    %SE 04.10.2012 21:58:31
    for clearb = 1:numWordCols
        hb=findobj('Tag',['afc_button' num2str(clearb)]);
        set(hb,'enable','off');
				
        hb=findobj('Tag',['afc_wordSelection' num2str(clearb)]);
        set(hb,'string','');
    end
   	drawnow; % flush event que 15-04-2005 13:38

case 'measure'
   hm=findobj('Tag','afc_message');
   set(hm,'string',msg.measure_msg);
   
case 'maxvar'
	hm=findobj('Tag','afc_message');
   set(hm,'string',msg.maxvar_msg);
   
case 'minvar'
   hm=findobj('Tag','afc_message');
	set(hm,'string',msg.minvar_msg);
 
case 'markint'						% must be a blocking action !!!
   tic;

	inter=max(0,def.intervallen/def.samplerate-0.02);
	paus=max(0,def.pauselen/def.samplerate-0.02);
	
	% 1.00.1
	if ( ~(strcmp(def.externSoundCommand, 'sndmex') & def.sndmexmark) & ~((strcmp(def.externSoundCommand, 'soundmex') | strcmp(def.externSoundCommand, 'soundmexfree') | strcmp(def.externSoundCommand, 'soundmex2') | strcmp(def.externSoundCommand, 'soundmex2free') | strcmp(def.externSoundCommand, 'soundmexpro') ) & def.soundmexMark) )
	%if ((def.sndmex == 0) | (def.sndmexmark == 0))
		switch def.markinterval
	   	case 1
			pause(def.presiglen/def.samplerate)
	
			for i=1:def.intervalnum-1
	   			h=findobj('Tag',['afc_button' num2str(i)]);
	   			if ( ~isempty(h) )
	   				set(h,'backgroundcolor',[1 0 0])
	   			end
	   			pause(inter)
	   			if ( ~isempty(h) )
					set(h,'backgroundcolor',[0.75 0.75 0.75])
	   			end	
	   			pause(paus)
			end
			h=findobj('Tag',['afc_button' num2str(def.intervalnum)]);
			if ( ~isempty(h) )
				set(h,'backgroundcolor',[1 0 0])
	   		end
	   		pause(inter)
	   		if ( ~isempty(h) )
	   			set(h,'backgroundcolor',[0.75 0.75 0.75])
	  		end
	  		pause(def.postsiglen/def.samplerate)
	   
	   	case 0
	      	%pause((def.presiglen+def.postsiglen+def.intervalnum*def.intervallen+ ...
	      	%(def.intervalnum-1)*def.pauselen)/def.samplerate)   
		end
	end

	elapsed = toc;												% blocking until end of signal presentation is reached
	%while elapsed < def.bgsiglen/def.samplerate+0.1	% plus 0.1 sec safety margin 
   	while elapsed < work.blockButtonTime
   		pause(0.02);
   		elapsed = toc;
		end
   
% ------------------------- end of action 'markint' -------------------
   
case 'markpressed'
   hb=findobj('Tag',['afc_button' num2str(work.answer{work.pvind}(work.stepnum{work.pvind}(end)))]);
   
   if ( ~isempty(hb) )
   	set(hb,'backgroundcolor',[0 0 0.75]);
   	pause(0.2);
	set(hb,'backgroundcolor',[0.75 0.75 0.75]);
	pause(0.0001);
   end

case 'markcorrect'
   hb=findobj('Tag',['afc_button' num2str(work.position{work.pvind}(work.stepnum{work.pvind}(end)))]);
   
   if def.markpressed
   	pause(0.2);
   end
   
   if ( ~isempty(hb) )
	   set(hb,'backgroundcolor',[0.5 0.5 0]);
	   pause(0.2);
	   set(hb,'backgroundcolor',[0.75 0.75 0.75]);
	   pause(0.0001);
   end

case 'response_ready'
   hm=findobj('Tag','afc_message');
	set(hm,'string',msg.ready_msg);
    
    %SE 04.10.2012 21:58:31
    for ( clearb = 1:5 )
				hb=findobj('Tag',['afc_button' num2str(clearb)]);
				set(hb,'enable','on');
				drawnow;
		end

end	% end switch

% eof
