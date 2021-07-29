#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon July 1 2020
Full Beacon (WWV / CHU) plotting of 2 input fomats
added hours ticks, removed time from UTC dat in header, added zero ref line for doppler freq shift, Elv now Elev  7-31=20 JCG
added autoplot capability 11/21/20 jgibbons
added autoxfer, autoplot on/off capability, fixed .png filename format problem  (RadioID <-> GridSqr) 1-29-2021
@authors dkazdan jgibbons

**********************************************
mod to plot up to 10 PSWS "rawdata" files and average value
files in PSWS subdir, leaves plot in Splot
plots files from multiple subdir to compare node results
plot title from first file
windows version for directories that mirror Pi
Bob Benedict, KD8CGH, 7/29/2021

create text file in PSWS directory
  keyword ('Doppler' or 'Power')
  subdir/filename1 
  subdir/filename2
  ...

if found 'Doppler' will plot Doppler shifts, if not, will plot Power

loads file names in list
plots first file and create axis and title info
plots rest in loop as curves on first plot
calculates average and plots

"""

import os
from os import path
import sys
import csv
import shutil
from datetime import date, timedelta
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import legend, show, grid, figure, savefig
import matplotlib.colors as mcolors
from scipy.signal import filtfilt, butter
import subprocess
from WWV_utility2 import time_string_to_decimals, graph_Doppler_and_power_data
import maidenhead as mh

names = open("E:\Documents\PSWS\\DirTestfiles.txt", "r")
PlotTarget = names.readline()
PlotTarget = PlotTarget.strip('\n')

Filenames=['a' for a in range (10)]
Filedates=['a' for a in range (10)]
PrFilenames=['a' for a in range (10)]

nfiles = 0
print(nfiles)

colors=['b','g','r','c','m','y','tab:orange','tab:gray','tab:purple','tab:brown']

while True:
    temp = names.readline()
    if  len(temp) == 0:
        break
    Filenames[nfiles]=temp.strip("\n")
    Filedates[nfiles]=temp[17:23]
    nfiles=nfiles + 1
        
print(Filenames[0:9])
print(Filedates[0:9])
print(nfiles)
if nfiles > 10 :
    print('10 file limit')
    sys.exit(0)
#
# ~ points to users home directory - usually /home/pi/
#homepath = os.path.expanduser('~')

# imbed the trailing / in the home path
#homepath = homepath + "/PSWS/"
homepath = "E:\\Documents\\PSWS\\"
#print('Home Path = ' + homepath)

# Define the base processing directory
#PROCESSDIR = homepath + 'Srawdata/'
PROCESSDIR = homepath
#directories for temp storing processed data files
DATADIR = homepath + 'Stemp/'

#saved data directrory
SAVEDIR = homepath + 'Sdata/'

#saved plot directrory
PlotDir = homepath + 'Splot/'

#saved info about system directory
InfoDir = homepath + 'Sinfo/'
#print('InfoDir = ' + InfoDir)

#transfer directory  (files to be sent to server node)
XferDir = homepath + 'Sxfer/'
#print('XferDir = ' + XferDir)

# flag file for auto plot on / off
autoplotfile = homepath + 'Scmd/autoplot'
#print('AutoPlot file path = ' + autoplotfile)

# flag file for auto transfer on / off
autoxferfile = homepath + 'Scmd/autoxfer'
#print('AutoPlot file path = ' + autoplotfile)

#-------------------------------------------------
'''
read first file
'''
PrFilenames=(PROCESSDIR + Filenames[0])

#-------------------------------------------------------------------
# Check to see if the file from yesterday exists; if not, exit the program
if (path.exists(PrFilenames)):
    print('File ' + PrFilenames + ' found!\nProcessing...')
else:
    print('File ' + PrFilenames + ' not available.\nExiting disappointed...')
    sys.exit(0)

with open(PrFilenames, 'r') as dataFile:
    dataReader=csv.reader(dataFile)
    data = list(dataReader)
    Header = data.pop(0)

    #Figure out which header format reading
    NewHdr = 'Unknown'
    print('Header to check=',Header)
    # Check if First header line is of new format example
    #,2020-05-16T00:00:00Z,N00001,EN91fh,41.3219273, -81.5047731, 284.5,Macedonia Ohio,G1,WWV5
    if (Header[0] == "#"):
        print('New Header String Detected')
        # Have new header format - pull the data fields out
        NewHdr = 'New'
        UTCDTZ = Header[1]
        print('\nUTCDTZ Original Header from file read = ' + UTCDTZ)
        UTC_DT = UTCDTZ[:10] # Strip off time and ONLY keep UTC Date
        print('\nExtracted UTC_DT only = ' + UTC_DT)
        UTCDTZ=UTCDTZ.replace(':','') # remove the semicolons
        print('\ncorrected UTCDTZ =', UTCDTZ)
        node= Header[2]
#        print('Node =', node)
        GridSqr = Header[3]
#        print('GridSqr =', GridSqr)
        Lat = Header[4]
#        print('Lat =', Lat)
        Long = Header[5]
#        print('Long =', Long)
        Elev = Header[6]
#        print('Elev =', Elev)
        citystate = Header[7]
#        print('City State =', citystate)
        RadioID = Header[8]
#        print('Radio ID =', RadioID)
        beacon = Header[9]
#        print('Beacon =', beacon)

    # Try using original FLdigi format w/o info line fake all data
#    if (Header[0] == "UTC"):
#        print('Detected Original FLDigi Header Format')
#        UTCDTZ = "2020-00-00T00:00:00Z"
#        node= "N00000"
#        Lat = '00.00000'
#        Long = '-00.00000'
#        Elev = '000'
#        print('GridSqr =', GridSqr)
#        citystate = 'NOcity NOstate'
#        RadioID = SRID
#        beacon = "Unknown"
#        NewHdr = 'FLDigi'

    if (NewHdr == 'Unknown'):
        ChkDate = Header[0]  # load in first row entry
        Cent = ChkDate[:2] # check first 2 digits  = 20?
        print( ChkDate, 'Header Yields Century of', Cent)  # diag printout

        if  Cent == "20":
            print('Old Header String Detected')
            # Have old header format - pull the data fields out
            #2020-05-15,N8OBJ Macedonia Ohio EN91fh,LB GPSDO,41.3219273, -81.5047731, 284.5
            UTCDTZ = Header[0]
            UTC_DT = UTCDTZ[:10] # Strip off time and ONLY keep UTC Date
            UTCDTZ=UTCDTZ.replace(':','') # remove the semicolons
            print('UTCDTZ =', UTCDTZ)
            #get this stations Node #
            Lat = Header[3]
            #print('Lat =', Lat)
            Long = Header[4]
            #print('Long =', Long)
            Elev = Header[5]
            #print('Elev =', Elev)
            GridSqr = mh.to_maiden(float(Lat), float(Long))
            print('GridSqr =', GridSqr)
            citystate = Header[1]
            #print('City State =', citystate)
            RadioID = 'G1'
            #print('Radio ID =', RadioID)
#            beacon = "Unknown"
#            print('Beacon =', beacon)
            NewHdr = 'Old'

    print('Header Decode =',NewHdr)
    #print('Scanning for UTC header line')

if (NewHdr == 'Unknown'):
    print('Unknown File header Structure - Aborting!')
    sys.exit(0)

print('Ready to start processing records')

# Prepare data arrays
hours=[[],[],[],[],[],[],[],[],[],[]]
Doppler=[[],[],[],[],[],[],[],[],[],[]]
Vpk=[[],[],[],[],[],[],[],[],[],[]]
Power_dB=[[],[],[],[],[],[],[],[],[],[]] # will be second data set, received power 9I20
filtDoppler=[[],[],[],[],[],[],[],[],[],[]]
filtPower=[[],[],[],[],[],[],[],[],[],[]]

LateHour=False # flag for loop going past 23:00 hours

# eliminate all metadata saved at start of file - Look for UTC (CSV headers)
#find first row of data0
FindUTC = 0
recordcnt  = 0
freqcalc = 0
calccnt = 0
plot=0

for row in data:
    if (FindUTC == 0):
        #print('looking for UTC - row[0] =',row[0])
        if (row[0] == 'UTC'):
            FindUTC = 1
#            print('UTC found =', row[0])
    else:
        #print('Processing record')
        decHours=time_string_to_decimals(row[0])
        if (NewHdr != 'New'):
            if (calccnt  < 101):
                calcnt = calcnt+1
                freqcalc = freqcalc + (float(row[1])/100)
#        if decHours > 23:
#            LateHour=True # went past 23:00 hours
        if (not LateHour) or (LateHour and (decHours>23)): # Otherwise past 23:59:59.  Omit time past midnight.
            hours[0].append(decHours) # already in float because of conversion to decimal hours.
            Doppler[0].append(float(row[2])) # frequency offset from col 2
            Vpk[0].append (float(row[3])) # Get Volts peak from col 3
            Power_dB[0].append (float(row[4])) # log power from col 4


print('nf ',0,'hours',len(hours[0]))

###############################################################################################
# Find max and min of Power_dB for graph preparation:
min_power=np.amin(Power_dB[0]) # will use for graph axis min
max_power=np.amax(Power_dB[0]) # will use for graph axis max
min_Vpk=np.amin(Vpk[0]) # min Vpk
max_Vpk=np.amax(Vpk[0]) # max Vpk
min_Doppler=np.amin(Doppler[0]) # min Doppler
max_Doppler=np.amax(Doppler[0]) # max Doppler

print('\nDoppler min: ', min_Doppler, '; Doppler max: ', max_Doppler)
print('Vpk min: ', min_Vpk, '; Vpk max: ', max_Vpk)
print('dB min: ', min_power, '; dB max: ', max_power)
#sys.exit(0) 
 
#%% Create an order 3 lowpass butterworth filter.
# This is a digital filter (analog=False)
# Filtering at .01 to .004 times the Nyquist rate seems "about right."
# The filtering argument (Wn, the second argument to butter()) of.01
# represents filtering at .05 Hz, or 20 second weighted averaging.
# That corresponds with the 20 second symmetric averaging window used in the 1 October 2019
# Excel spreadsheet for the Festival of Frequency Measurement data.
#FILTERBREAK=.005 #filter breakpoint in Nyquist rates. N. rate here is 1/sec, so this is in Hz.
FILTERBREAK=0.005 #filter breakpoint in Nyquist rates. N. rate here is 1/sec, so this is in Hz.
FILTERORDER=6
b, a = butter(FILTERORDER, FILTERBREAK, analog=False, btype='low')
#print (b, a)
#%%
# Use the just-created filter coefficients for a noncausal filtering (filtfilt is forward-backward noncausal)

#print ('Filter Doppler shift data')
filtDoppler[0] = filtfilt(b, a, Doppler[0])

#print ('Filter power data')
filtPower[0] = filtfilt(b, a, Power_dB[0])
#sys.exit(0)
##%% modified from "Double-y axis plot,
## http://kitchingroup.cheme.cmu.edu/blog/2013/09/13/Plotting-two-datasets-with-very-different-scales/

##################################################################################################


# set up x-axis with time
fig = plt.figure(figsize=(19,10)) # inches x, y with 72 dots per inch
ax = fig.add_subplot(111)
ax.set_xlabel('UTC Hour')
ax.set_xlim(0,24) # UTC day
ax.set_xticks([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24], minor=False)

# plot first curve
if (PlotTarget == 'Doppler'):    
    ax.plot(hours[0], filtDoppler[0], colors[0],label=Filedates[0]) # color k for black
    ax.set_ylabel('Doppler shift, Hz '+ Filedates[0])
    ax.set_ylim([-1.0, 1.0]) # -1 to 1 Hz for Doppler shift
    plt.axhline(y=0, color="gray", lw=1) # plot a zero freq reference line for 0.000 Hz Doppler shift
else:
    ax.plot(hours[0], filtPower[0], colors[0],label=Filedates[0]) # color k for black
    ax.set_ylabel('Power, dB '+ Filedates[0])
    ax.set_ylim(-90, 0)    

# add grid lines - RLB
plt.grid(axis='both')
'''
######################################################################
read and plot files loop
'''

for nf in range(1, nfiles):
# splot second curve
# read second file, skip header
    print('process file ',nf, Filenames[nf])
    PrFilenames=(PROCESSDIR + Filenames[nf])
    with open(PrFilenames, 'r') as dataFile: # read second set
        dataReader=csv.reader(dataFile)
        data = list(dataReader)
    FindUTC = 0
    recordcnt  = 0
    freqcalc = 0
    calccnt = 0

    for row in data:
        if (FindUTC == 0):
            #print('looking for UTC - row[0] =',row[0])
            if (row[0] == 'UTC'):
                FindUTC = 1
#            print('UTC found =', row[0])
        else:
            decHours=time_string_to_decimals(row[0])
            hours[nf].append(decHours) # already in float because of conversion to decimal hours.
            Doppler[nf].append(float(row[2])) # frequency offset from col 2
            Vpk[nf].append (float(row[3])) # Get Volts peak from col 3
            Power_dB[nf].append (float(row[4])) # log power from col 4    

# filter second file data
    filtDoppler[nf] = filtfilt(b, a, Doppler[nf])
    filtPower[nf] = filtfilt(b, a, Power_dB[nf])

#    print('nf ',nf,'hours',len(hours[nf]))
#    print('filtDoppler',len(filtDoppler[nf]))
#    print(filtDoppler[nf][0:9])
#    print(hours[nf][0:9])

#ax2 = ax1.twinx()
    if (PlotTarget == 'Doppler'):    
        ax.plot(hours[nf], filtDoppler[nf], colors[nf], label=Filedates[nf]) # color k for black
    else:
        ax.plot(hours[nf], filtPower[nf], colors[nf], label=Filedates[nf]) # color k for black

# following lines set ylim for power readings in file
#ax2.set_ylim(min_power, max_power) #as determined above for this data set
#for tl in ax2.get_yticklabels():
#    tl.set_color('r')

'''
#############################################################################
end for read and plot loop, start average
'''
# find shortest data set, limit average to that
al=1000000
ak=0

for k in range(nfiles):
    templ=len(hours[k])
    if templ < al:
        al=templ
        ak=k

avg=[]

if (PlotTarget == 'Doppler'):
    for i in range(al):
        temp=0.0
        for j in range(nfiles):
            temp=temp+filtDoppler[j][i]
        temp=temp/nfiles
        avg.append(temp)
else:
    for i in range(al):
        temp=0.0
        for j in range(nfiles):
            temp=temp+filtPower[j][i]
        temp=temp/nfiles
        avg.append(temp)



print('avg',len(avg))


ax.plot(hours[ak], avg, 'k', label='Average') # color k for black


'''
end average
'''

ax.legend(loc="lower right", title="Legend Title", frameon=False)

#2.5MHz WWv
if (beacon == 'WWV2p5'):
    print('Final Plot for Decoded 2.5MHz WWV Beacon\n')
    beaconlabel = 'WWV 2.5 MHz'

#5MHz WWV
elif (beacon == 'WWV5'):
    print('Final Plot for Decoded 5MHz WWV Beacon\n')
    beaconlabel = 'WWV 5 MHz'

#10MHz WWV
elif (beacon == 'WWV10'):
    print('Final Plot for Decoded 10MHz WWV Beacon\n')
    beaconlabel = 'WWV 10 MHz'

#15MHz WWV
elif (beacon == 'WWV15'):
    print('Final Plot for Decoded 15MHz WWV Beacon\n')
    beaconlabel = 'WWV 15 MHz'

#20MHz WWV
if (beacon == 'WWV20'):
    print('Final Plot for Decoded 20MHz WWV Beacon\n')
    beaconlabel = 'WWV 20 MHz'

#25MHz WWV
elif (beacon == 'WWV25'):
    print('Final Plot for Decoded 25MHz WWV Beacon\n')
    beaconlabel = 'WWV 25 MHz'

#3.33MHz CHU
if (beacon == 'CHU3'):
    print('Final Plot for Decoded 3.33MHz CHU Beacon\n')
    beaconlabel = 'CHU 3.330 MHz'

#7.85MHz CHU
elif (beacon == 'CHU7'):
    print('Final Plot for Decoded 7.85MHz CHU Beacon\n')
    beaconlabel = 'CHU 7.850'

#14.67MHz CHU
elif (beacon == 'CHU14'):
    print('Final Plot for Decoded 14.67MHz CHU Beacon\n')
    beaconlabel = 'CHU 14.670 MHz'

elif (beacon == 'Unknown'):
    print('Final Plot for Decoded Unknown Beacon\n')
    beaconlabel = 'Unknown Beacon'

#PlotDir = homepath + 'Splot/'
# Create Plot Title
plt.title(beaconlabel + ' Grape Data Plot\nNode:  ' + node + '     Gridsquare:  '+ GridSqr + '\nLat= ' + Lat + '    Long= ' + Long + '    Elev= ' + Elev + ' M\n' )

# Create Plot File Nam
#GraphFile = yesterdaystr + '_' + node + '_' + RadioID + '_' + GridSqr + '_' + beacon + '_graph.png'
GraphFile = PlotTarget + '_' + node + '_' + RadioID + '_' + GridSqr + '_' + beacon + '.png'
PlotGraphFile = PlotDir + GraphFile
XferGraphFile = XferDir + GraphFile

# create plot
#plt.savefig(PlotDir + yesterdaystr + '_' + node + '_' +  GridSqr + '_' +  RadioID + '_' +  beacon + '_graph.png', dpi=250, orientation='landscape')
plt.savefig(PlotDir + GraphFile, dpi=250, orientation='landscape')
# =============================================================================

print('Plot File: ' + GraphFile + '\n')  # indicate plot file name for crontab printout


#-------------------------------------------------------------------
# Check to see if the autoxfer file exists; if not, skip copying the plot
if (path.exists(autoxferfile)):
    print('autoxfer enable File found\n')
    # copy plot file to transfer directory
    print('Copying file to Sxfer directory for upload ',XferGraphFile)  # Copy final processed data file also to /Sxfer/ directory for transfer
    shutil.copy(PlotGraphFile, XferGraphFile)

else:
    print('No autoxfer enable File found - exiting')

#-------------------------------------------------------------------
# Check to see if the autoplot file exists; if not, skip the plot
if (path.exists(autoplotfile)):
    print('autoplot enable File found - Processing Plot...\n')
    subprocess.call('gpicview ' + PlotGraphFile +' &', shell=True)   #create shell and plot the data from graph file
else:
    print('No autoplot enable File found - exiting')
#-----------------------------------------------

#subprocess.call('gpicview ' + GraphFile +' &', shell=True)   #create shell and plot the data from graph file

print('Exiting python combined processing program gracefully')
