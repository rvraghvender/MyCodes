#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 25 09:59:01 2021

@author: Raghvender
"""

import re
import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats
from decimal import Decimal

InitialValue = 5 # ps
FinalValue = 8  # ps


def statistics(y):
    slope, intercept, r_value, p_value, std_err = stats.linregress(MSD_Fit['Time'], MSD_Fit[y])    
    Diffusion = f"{Decimal( float(slope/60000)     ):.2E}"
    line = slope * MSD_Fit['Time'] + intercept
    return Diffusion, line



with open('msd.dat','r') as ReadFile:
    msd = ReadFile.readlines()
    
    # Reading all the elemets
    elements = re.findall(r"\[([A-Za-z ]+)\]", msd[0])[1:]
    
    # Removing blank space from each element of list
    elements = [name.strip() for name in elements]    
    

# Creating empty dictionary
MSD_data = {'Time':[]}
for element in elements:
    MSD_data[element] = []
    
    
for data in msd[1:]:
    x = [float(value) for value in data.split()]

    # Appending values to dictionary
    i = 0
    for key in MSD_data.keys():
        MSD_data[key].append(x[i])
        i=i+1
        
        
# Converting dictionary to DataFrame and Plotting      
MSD_data = pd.DataFrame(MSD_data)

MSD_Fit = MSD_data.loc[(MSD_data['Time'] > InitialValue) & (MSD_data['Time'] < FinalValue) ]



ax = plt.gca()
for element in elements:
    Diffusion, line = statistics(element)
    MSD_data.plot(kind='line',x='Time',y=element, ax=ax, label=f"{Diffusion} cm2/s : {element}")
    plt.plot(MSD_Fit['Time'], line, 'k--')
plt.xlabel('Time (ps)')
plt.ylabel('MSD')
plt.legend()
plt.show()

