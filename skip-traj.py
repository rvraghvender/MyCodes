#!/usr/bin/python3.7
#program to calculate the partial distribution funtion and plot it                    ##
import platform                                                                       ##
import sys                                                                            ##
print()
print()                                                                               ##
print("System Information : ",sys.version, sys.platform)                              ##
print("Build Info : ",platform.machine(),platform.system())                           ##
print()                                                                               ##
print()                                                                               ##
                                                                                      ##
import os                                                                             ##
cwd = os.getcwd()  # Get the current working directory (cwd)                          ##
files = os.listdir(cwd)  # Get all the files in that directory                        ##
print("Available Files in directory %r : " % (cwd))                                   ##
print()                                                                               ##
for counter, file_items in enumerate(files,1):                                        ##
    print(counter,"--" ,file_items)                                                   ##
print()                                                                               ##
print("--------------------------------------------------------------------------")   ##
                                                                                      ##
print()                                                                               ##
print()                                                                               ##
                                                                                      ##
########################################################################################





import numpy as np


def parse_input(filename='',skip = 10):  # for getting coordinates from the trajectory file

    with open(filename, 'r') as fo: # opening file
        data = fo.readlines()   # reading each line

    n_atoms = int(data[0].split()[0]) 
    print(n_atoms)
    n_steps_total = int(len(data) / (n_atoms +2)) 
  
    n_steps = n_steps_total // skip

    
    countt = 0

    atom_list = [line.split()[0] for line in data[2: n_atoms + 2]]  
    Atom_type = list(dict.fromkeys(atom_list)) ## types of atom


    for step in range(n_steps): # here i am creating a loop with nsteps

        coords = np.zeros((n_atoms, 4), dtype = object)
        i = step * skip * (n_atoms + 2) 
        for j, line in enumerate(data[i + 2 : i + n_atoms + 2]): 
            coords[j, :] = [value for value in line.split()[:]]

        
        fw.write("{}\n".format(n_atoms))
        fw.write("{}".format(data[i + 1]))

        for k in range(n_atoms):
            fw.write(" {}       {}      {}      {} \n ".format(coords[k,0],coords[k,1],coords[k,2],coords[k,3]))

fo1 = input("Enter the name of trajectory file (.xyz) format : ")
skip_ts = int(input("Enter the steps to jump : "))

fw = open("{}".format("traj.xyz"),'w')
parse_input(fo1, skip_ts)

print()

print("Your file has been created succesfully")

print()
print()
