#!/usr/bin/python3.7
#program to  find the bond length of oxides                                           ##


import platform                                                                       ##
import sys                                                                            ##
print()
print()                                                                               ##
print("System Information : ",sys.version, sys.platform)                              ##
print("Build Info : ",platform.machine(),platform.system())                           ##
print()                                                                               ##
print()                                                                               ##
                                                                                      ##
#import os
import sys##
# cwd = os.getcwd()  # Get the current working directory (cwd)                          ##
# files = os.listdir(cwd)  # Get all the files in that directory                        ##
# print("Available Files in directory %r : " % (cwd))                                   ##
# print()                                                                               ##
# for counter, file_items in enumerate(files,1):                                        ##
#     print(counter,"--" ,file_items)                                                   ##
# print()                                                                               ##
# print("--------------------------------------------------------------------------")   ##
#                                                                                       ##
print()                                                                               ##
print()                                                                               ##
                                                                                      ##
########################################################################################


print()


import numpy as np


def parse_input(filename=''):  # for getting coordinates from the trajectory file

    global b
    global coords
    global clust


    with open(filename, 'r') as fo: # opening file
        data = fo.readlines()   # reading each line

    n_atoms = int(data[0].split()[0]) 


    atom_list = [line.split()[0] for line in data[2: n_atoms + 2]]  
    Atom_type = list(dict.fromkeys(atom_list)) ## types of atom

    print("Atoms present in the trajectory file are : " ,Atom_type)
    print ()

    coords = np.zeros((n_atoms, 4), dtype = object)

    for j, line in enumerate(data[2 : n_atoms + 2]): 
        coords[j, 0:4] = [value for value in line.split()[:]]

    atom1_oxide = []
    atom2_oxide = []
    for i, value in enumerate(coords):


        if coords[i,0] == "O" :

            xo = float(coords[i,1])
            yo = float(coords[i,2])
            zo = float(coords[i,3])

            # print(coords[i, :])

            
            for j, value in enumerate(coords):
                if coords[j,0] == "Tl":

                    x1 = float(coords[j,1])
                    y1 = float(coords[j,2])
                    z1 = float(coords[j,3])
                
                    dist = np.sqrt( (x1-xo) ** 2 + (y1-yo)**2 + (z1-zo)**2 )

                    # print(dist)

                    if dist <= 3.27:
                        atm_oxide = round(dist,3)
                        atom1_oxide.append(atm_oxide)



                if coords[j,0] == "Te":

                    x1 = float(coords[j,1])
                    y1 = float(coords[j,2])
                    z1 = float(coords[j,3])

                    dist = np.sqrt( (x1-xo) ** 2 + (y1-yo)**2 + (z1-zo)**2 )

                    # print(dist)

                    if dist <= 3.2:
                        atm_oxide2 = round(dist,3)
                        atom2_oxide.append(atm_oxide2)



    for i in range(len(atom1_oxide)):
        print('Tl-O',atom1_oxide[i])
    print()

    for i in range(len(atom2_oxide)):
        print('Te-O',atom2_oxide[i])
try:
    fo1=sys.argv[1]
except:
    print("Please write script name followed by trajectory name (xyz)  ")
    sys.exit()



parse_input(fo1)
print()
print("Done")

print()
print()
