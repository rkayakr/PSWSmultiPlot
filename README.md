# PSWSmultiPlot
multiPlot version 1.2 cleans up to make Pi and Windows use easier
multiPlot version 1.1 adds sunrise sunset times for location of lat long in first file
requires suntime library https://github.com/SatAgro/suntime
multiPlot version v1.0 plots up to 10 PSWS "rawdata" files and average value
modified from WWV_plt2.py @authors dkazdan jgibbons
expects a homepath directory with raw files in subdirs, leaves plot in Mplot directory
plots files from multiple subdir to compare node results
plot title from first file

windows version hardcoded homepath directory location
for Pi comment out windows homepath and uncomment Pi lines

uses WWV_utility2.py
Bob Benedict, KD8CGH, 7/29/2021

create text file "plotfiles.txt" in homepath directory
  keyword ('Doppler' or 'Power')
  subdir/filename1 
  subdir/filename2
  ...
if found 'Doppler' will plot Doppler shifts, else will plot Power
loads file names in list
plots first file and create axis and title info
plots rest in loop as curves on first plot
calculates average and plots
leaves plotfile in Mplot directory

uses WWV_utility2.py
20 February 2020
WWV utility file
Routines and classes used in WWV file management and graphing
David Kazdan, AD8Y
John Gibbons, N8OBJ - mods to plot header 2/3/2020
