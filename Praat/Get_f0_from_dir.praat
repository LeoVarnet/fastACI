# Praat version : 5.3.76 for Windows
# Created on    : 11/08/2014
# Last update on: 11/08/2014
# Description   : This file computes F1 and F2 for each segment defined previously in TextGrid (manually created)

form exportar 
	text dir D:\MATLAB\Output\tmp-VoD_MIRtoolbox\
	positive timestep 0.01
	positive pitchfloor 100
	positive pitchceiling 500
endform

# Make a list consisting of the names of files in the directory
Create Strings as file list... list 'dir$'*.wav

# count the number of names and store the number
number_files = Get number of strings

for i from 1 to number_files
	
	select Strings list
	name$ = Get string... 'i'
	#name$ = selected$("wav")
	#Load files from directory
	
	Read from file... 'dir$''name$'
	
	nlength = length(name$)
	name_out$ = left$ (name$,nlength-4)

	outfile$ = name_out$ + "_f0.txt"
	#If the output file already exists, delete it
	filedelete 'dir$'\'outfile$'

	#In output file, add a line with name, duration, f0 values
	fileappend 'dir$'\'outfile$' 'name' 'tab$' 'time' 'tab$' 'f0' 'newline$'

	sound = selected ("Sound")
	tmin = Get start time
	tmax = Get end time
	To Pitch... timestep pitchfloor pitchceiling
	clearinfo
	#echo Here are the results:
	for j to (tmax-tmin)/timestep
		timeinitial = tmin + j * timestep
		durms = timestep*1000
		# name without extension:
		select Pitch 'name_out$'
		pitch = Get value at time... 'timeinitial' Hertz linear
		fileappend 'dir$'\'outfile$' 'timeinitial:3' 'tab$' 'pitch:1' 'newline$'
		endif

	endfor 

endfor
#select TextGrid 'name$'
#plus Sound 'name$'
#plus Formant 'name$'
#Remove 
