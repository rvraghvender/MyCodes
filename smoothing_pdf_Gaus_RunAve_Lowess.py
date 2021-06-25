import os
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
from statsmodels.nonparametric.kernel_regression import KernelReg

#!!!!!!!!!!_modify_here_if_needed_!!!!!!!!!!!!!!!
#Choose method                                 #!
Gaussian = True                                #!
SIGMALISS = 0.05                               #!
                                               #!
RunningAverage = True                          #!
movingwidth = 4                                #!
                                               #!
Lowess = False                                  #!
Lwidth = 0.02                                  #!
                                               #!
showplot = True                                #!
printsmoothed = True                           #!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#Check that you are doing something useful!
try:
    print(3/0)
except:
    if not showplot and not printsmoothed:
        print('ERROE: I am not doing anything!')
        exit()

#Open input file 
try:
    file = open('gr-xrays.dat', 'r') 
except:
    print('ERROE: gr-xrays.dat not found')
    exit()

#Open output file to save the smoothed data
if printsmoothed:
    fileout = open('gr-xrays_smoothed.dat', 'w') 
    
#define lists for r, gr
r = list()
gr = list()

#read r, gr from input
for line in file:
    if 'diffraction' not in line:
        line = line.rstrip()
        line = line.split()
        r.append(float(line[0]))
        gr.append(float(line[1]))
if showplot:
    plt.scatter(r, gr, label='Original data', marker='.',c='k')

#Method 1: Gaussian smoothing
if Gaussian:
    DIMTOLISS = len(r) - 1 
    FACTLISS = 1 / (SIGMALISS * np.sqrt(2 * np.pi ))
    smoothedR = list()

    for INDA in range(DIMTOLISS):
        xx = 0
        RRR = []
        for INDB in range(DIMTOLISS):
            if INDB == 0:
                DQ = r[1] - r[0]
            elif INDB == DIMTOLISS:
                DQ = r[DIMTOLISS] - r[DIMTOLISS-1]
            else:
                 DQ = 0.5 * (r[INDB + 1] - r[INDB - 1])
            xx = np.exp(-(r[INDB] - r[INDA])**2 / ( 2 * SIGMALISS **2)) * gr[INDB] * DQ 
            RRR.append(xx)
        smoothedR.append(sum(RRR) * FACTLISS)
    if showplot:
        plt.plot(r[:DIMTOLISS-1],smoothedR[:DIMTOLISS-1], label='Gaussian Smoothing')
    if printsmoothed:
        print('#Gaussian smoothing', file=fileout)
        print('# r(ang) gr', file=fileout)
        for c, cc in list(zip(np.array(r),np.array(smoothedR))):
            print(c,cc, file=fileout)

#Method 2: Moving average box (by convolution)
if RunningAverage:
    def smooth(y, box_pts):
        box = np.ones(box_pts)/box_pts
        y_smooth = np.convolve(y, box, mode='same')
        return y_smooth
    SmoothedGR = list()
    SmoothedGR = smooth(gr, movingwidth)
    if showplot:
        plt.plot(r, SmoothedGR, label='Moving average box')
    if printsmoothed:
        print(' ', file=fileout)
        print('#Moving average box ', file=fileout)
        print('# r(ang) gr', file=fileout)
        for cx, ccx in list(zip(np.array(r),np.array(SmoothedGR))):
            print(cx,ccx, file=fileout)

#Method 3: Lowess similar to moving average
if Lowess:
    lowess = sm.nonparametric.lowess(gr, r, frac=Lwidth)
    if showplot:
        plt.plot(lowess[:, 0], lowess[:, 1], label='Lowess Smoothing')
    if printsmoothed:
        print(' ', file=fileout)
        print('#Lowess Smoothing ', file=fileout)
        print('# r(ang) gr', file=fileout)
        for cxx, ccxx in list(zip(lowess[:, 0], lowess[:, 1])):
            print(cxx,ccxx, file=fileout)

#set legend and show plot
if showplot:
    plt.xlabel('r (Ang)' )
    plt.ylabel('gr-xrays-tot')
    plt.legend()
    plt.show()
