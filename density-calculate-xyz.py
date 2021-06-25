#!/usr/bin/python3.7
# program to find the density from xyz file
import platform
import sys



print("System Information : ",sys.version, sys.platform)
print(platform.machine(),platform.system())
print()

import os
cwd = os.getcwd()  # Get the current working directory (cwd)
files = os.listdir(cwd)  # Get all the files in that directory
print("Available Files in directory %r : " % (cwd))

for counter, file_items in enumerate(files,1):
    print(counter,"--" ,file_items)

print("--------------------------------------------------------------------------")

print()
print()

##################################################################################

print("This program finds the number density (N/V).")
print()

f = open(input("Enter the name of the xyz file : "),'r')
print()
number_of_atom = 0

Atom_name = []
Atom_x = []
Atom_y = []
Atom_z = []


count = 0  #to read lattice vectors of xyz file
while True:
    count += 1
    lines = f.readline().split()
    if count == 2 :
        lx = float(lines[0])
        ly = float(lines[1])
        lz = float(lines[2])
        print()
        break


while True: # to read atomic coordinates
    line = f.readline().split()
    if len(line) == 0: break
    if len(line) < 4: continue
    number_of_atom += 1
    na, x, y, z = [a for a in line[:4]]
    Atom_name.append(na)
    Atom_x.append(float(x))
    Atom_y.append(float(y))
    Atom_z.append(float(z))

Atom_type = list(dict.fromkeys(Atom_name))

#max_x = max(Atom_x)
#min_x = min(Atom_x)

#max_y = max(Atom_y)
#min_y = min(Atom_y)

#max_z = max(Atom_z)
#min_z = min(Atom_z)

print("There are", number_of_atom , "atoms of type : \n", Atom_type)
print()
print()
print("Cell parameters are :")
print(lx,ly,lz)
#print( )
#print("Max length in x :", max_x," \n Min length in x :", min_x )
#print("Length (x) : ", max_x - min_x)
#print()
#print("Max length in y :", max_y ," \n Min length in y :", min_y )
#print("Length (y) : ", max_y - min_y)
#print()
#print("Max length in z :", max_z ," \n Min length in z :", min_z )
#print("Length (z) : ", max_y - min_y)


################# Density calculations####
print()
vol = lx * ly * lz
print("Volume is : ",  vol)

print()
print("Number density is : ", number_of_atom / vol)




f.close()
