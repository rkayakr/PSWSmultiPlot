# PSWSmultiPlot
	program to plot up to 10 PSWS "rawdata" files and average value
	based on "WWV_plt2.py" by dkazdan & jgibbons
	windows version for hard coded directory PSWS
	files in a PSWS subdirs, leaves plot in Mplot dir
	can plots files from multiple subdir to compare different node results for same day
	plot title taken from first file in plotfiles.txt
	uses WWV_utility2.py
		Bob Benedict, KD8CGH, 7/29/2021
 
	create text file "plotfiles" in PSWS directory
		keyword ('Doppler' or 'Power')
		Average (if "Average" will plot average of data, else not
		subdir/filename1 
		subdir/filename2
		...
 
	if found 'Doppler' will plot Doppler shifts, if not, will plot Power 
		loads file names in list
		plots first file and create axis and title info
		plots rest in loop as curves on first plot
		writes date for each curve in legend
		calculates average and plots
