#!/usr/bin/python3.7
#program to  sort xyz file into stochiometry pattern                                  ##


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


def parse_input(filename=''):  # for getting coordinates from the trajectory file

    with open(filename, 'r') as fo: # opening file
        data = fo.readlines()   # reading each line

    n_atoms = int(data[0].split()[0]) 
    
    n_factor = n_atoms / ( icoeff1 + icoeff3 + icoeff4)

    atom_list = [line.split()[0] for line in data[2: n_atoms + 2]]  
    Atom_type = list(dict.fromkeys(atom_list)) ## types of atom



    coords = np.zeros((n_atoms, 4), dtype = object)

    for j, line in enumerate(data[2 : n_atoms + 2]): 
        coords[j, :] = [value for value in line.split()[:]]

    fw.write("{}\n".format(n_atoms))
    fw.write("{}".format(data[1]))


    count = 0
    while count <= n_factor:
        itera1 = 1
        itera2 = 1
        itera3 = 1
        itera4 = 1
        for k in range(n_atoms):
            if coords[k,0] == "Ti" and itera1 <= icoeff1 :
                fw.write(" {}       {}      {}      {}\n".format(coords[k,0],coords[k,1],coords[k,2],coords[k,3]))
                coords[k,0] = [0,0,0,0]
                itera1 = itera1 + 1

        for k in range(n_atoms):
            if coords[k,0] == "Tl" and itera2 <= icoeff2 :
                fw.write(" {}       {}      {}      {}\n".format(coords[k,0],coords[k,1],coords[k,2],coords[k,3]))
                coords[k,0] = [0,0,0,0]
                itera2 = itera2 + 1


        for k in range(n_atoms):
            if coords[k,0] == "Te" and itera3 <= icoeff3 :
                fw.write(" {}       {}      {}      {}\n".format(coords[k,0],coords[k,1],coords[k,2],coords[k,3]))
                coords[k,0] = [0,0,0,0]
                itera3 = itera3 + 1

        for k in range(n_atoms):
            if coords[k,0] == "O" and itera4 <= icoeff4 :
                fw.write(" {}       {}      {}      {}\n".format(coords[k,0],coords[k,1],coords[k,2],coords[k,3]))
                coords[k,0] = [0,0,0,0]
                itera4 = itera4 + 1

        count = count + 1

fo1 = input("Enter the name of trajectory file (.xyz) format : ")
fw = open("{}-sort".format(fo1),'w')

icoeff1=int(input("Enter the coefficient of Ti : "))
icoeff2=int(input("Enter the coefficient of Tl : "))
icoeff3=int(input("Enter the coefficient of Te : "))
icoeff4=int(input("Enter the coefficient of O  : "))





parse_input(fo1)
print()
print("Your files has been created succesfully")

print()
print()
