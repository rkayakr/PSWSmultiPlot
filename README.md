# PSWSmultiPlot
# mod to plot up to 10 PSWS "rawdata" files and average value
# files in a PSWS subdir, leaves plot in Splot
# can plots files from multiple subdir to compare different node results for same day
# plot title from first file
# windows version for directories that mirror Pi
# uses WWV_utility2.py
# Bob Benedict, KD8CGH, 7/29/2021
# 
# create text file "plotfiles" in PSWS directory
#    keyword ('Doppler' or 'Power')
#    subdir/filename1 
#    subdir/filename2
#    ...
# 
# if found 'Doppler' will plot Doppler shifts, if not, will plot Power
# 
# loads file names in list
# plots first file and create axis and title info
# plots rest in loop as curves on first plot
# writes date for each curve in legend
# calculates average and plots
