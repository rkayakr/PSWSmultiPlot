# PSWSmultiPlot
mod to plot up to 10 PSWS "rawdata" files and average value
files in PSWS subdir, leaves plot in Splot
plots files from multiple subdir to compare node results
plot title from first file
windows version for directories that mirror Pi
Bob Benedict, KD8CGH, 7/29/2021

create text file "files" in PSWS directory
  keyword ('Doppler' or 'Power')
  subdir/filename1 
  subdir/filename2
  ...

if found 'Doppler' will plot Doppler shifts, if not, will plot Power

loads file names in list
plots first file and create axis and title info
plots rest in loop as curves on first plot
calculates average and plots
