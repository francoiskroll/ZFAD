# 210413_psen1LONGTRACK

## Experiment

### Info

dof 09/04/2021

Started tracking 13/04/2021 at 4 dpf

BOX12/13

* Sensitivity 20
* Burst 200
* Freezing 3
* Misc 400

DAY = 95% (range Up to 1040 lux)  
NIGHT = 0% (range Up to 1040 lux)  
DIM FREERUNNING = 5% (range Up to 1040 lux)

Irrigation device from Amazon  
Set to 50 mL every 4 hours  


### 21/04/2021 12dpf  
09.41  
BOX13: one fish escaped; it is swimming below the plate  

### 22/04/2021 13dpf  

STOP at 10.48 = 208h30m15s  

## Analysis

For purpose of thesis, plotting BOX13 SCR only.  

### Zebralab rawoutput  

3,553 xls files with 1M rows each  
i.e. 3553000000 rows total  
/ 192 wells / 25 frames-per-second / 60 seconds / 60 minutes = ~ 205 hours  

### Fixing ordering errors

With Vp_Sorter.R (ran on Rihel lab common computer, Windows), ran on each box at a time.  

Writes (for each box)
* 210413_12_RAWs.csv + 210413_13_RAWs.csv
* 210413_12_VpSorterlog.md + 210413_13_VpSorterlog.md

### Frame-by-frame data to trace

See 210413_PSEN1_main.R


## 210907_PSEN2

! used backup RNPs from -80C freezer
(not used before)

Ran test for one night 06/09/2021 to 07/09/2021 in Box 12/13 with Hobo  

Hobo data looks perfect: transitions at correct times and temperature fluctuates between 25.8 and 26.2˚C

### tracking

Stopped 10/09/2021 at 10.32  
total 67 hours 48 minutes  
All fish looking OK on video  
Water level OK in both boxes  
Temperature in room 27˚C  
In baths 27–27.5˚C  

! BOX12 is written as box1 in Zebralab, but prefix of wells is 2-  
and vice-versa for BOX13  
Check looking at the data of the right box (left well 96 of BOX12 empty for this purpose)  