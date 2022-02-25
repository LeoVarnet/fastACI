% Intellitest_srt_savefcn - custom save function for Intellitest procedure

function Intellitest_srt_savefcn

global def
global work

str=[def.result_path work.filename];

%%%%%%%%% build headerline %%%%%%%%
if def.interleaved > 0
   headerl = '% track#';
else
   headerl = '%';
end
for i = 1:def.expparnum
   eval(['parunit = def.exppar' num2str(i) 'unit;']);
   headerl = [headerl '   exppar' num2str(i) '[' parunit ']'];
end
headerl = [headerl '   expvar[' def.expvarunit ']'];
%%%%%%%% end headerline %%%%%%%%%

% write experimental variable at all reversals during the measurement phase to disk
%res=[work.exppar work.expvarrev(end - def.reversalnum+1:end)];
%
%r=[];
%for i=1:def.reversalnum
%   r=[r '   %8.8f'];
%end

% write mean and standard deviation to disk  

%res = work.expvarrev(end - def.reversalnum+1:end);
%res=[work.exppar mean(res) std(res,1)];

tmpinterleavenum = 1;
if def.interleaved > 0
   tmpinterleavenum = def.interleavenum;
end

for i = 1:tmpinterleavenum
	 	
    %SE 22.01.2013 18:18 save SNR and target intelligibility
    thresVl = FitVlLX(work.expvar{i}(1:end), work.answer{i}(1:end), def.TargetIntelligibility{i});
	res{i}=[thresVl def.TargetIntelligibility{i}];
    
	r='   %8.8f   %8.8f';
	
end

% write median, lower and upper quartiles to disk
% not yet implemented

ex = exist([str '.dat']);

if ex == 0
   %dat=['% exppar[' def.exppar1unit ']   expvar[' def.expvarunit ']'];
   fid=fopen([str '.dat'],'w');
   fprintf(fid,['%s\n'],headerl);
   for k=1:tmpinterleavenum
      if def.interleaved > 0
         fprintf(fid,'%i',k);
      end
   	for i=1:def.expparnum
      	eval(['tmp = work.int_exppar' num2str(i) '{' num2str(k) '};']);
      	if def.exppartype(i) == 0
         	fprintf(fid,['   %8.8f'],tmp);
      	else
         	fprintf(fid,['   %s'],tmp);
      	end
      end
      fprintf(fid,[r '\n'],res{k});
   end
   
   %fprintf(fid,['%8.8f' r '\n'],res);
   fclose(fid);
else
   fid=fopen([str '.dat'],'a');
   for k=1:tmpinterleavenum
      if def.interleaved > 0
         fprintf(fid,'%i',k);
      end
   	for i=1:def.expparnum
      	eval(['tmp = work.int_exppar' num2str(i) '{' num2str(k) '};']);
      	if def.exppartype(i) == 0
         	fprintf(fid,['   %8.8f'],tmp);
      	else
         	fprintf(fid,['   %s'],tmp);
      	end
      end
      fprintf(fid,[r '\n'],res{k});
   end

   % fprintf(fid,['%8.8f' r '\n'],res);
   fclose(fid);
end

if def.debug == 1
    save([str '.mat'], 'work'); 
end

% eof
