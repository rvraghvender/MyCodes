#!/usr/bin/python3.7


import platform  ##
import sys  ##

print()
print()  ##
print("System Information : ", sys.version, sys.platform)  ##
print("Build Info : ", platform.machine(), platform.system())  ##
print()  ##
print()  ##
##
import os

print()  ##
print()  ##

A =  float(input("Enter lattice parameter A : "))  # Entering cell parameters A,B,C
B =  float(input("Enter lattice parameter B : "))
C =  float(input("Enter lattice parameter C : "))

print()
cutoff =  float(input("Enter the cutoff distance for the desired atom : "))


import numpy as np


def distance(a, b):
    dx = a[0] - b[0]
    abs_x = abs(dx)
    if abs_x > A:
        dx = abs_x % A
    x = min(abs_x, (A - abs_x))

    dy = a[1] - b[1]
    abs_y = abs(dy)
    if abs_y > B:
        dy = abs_y % B
    y = min(abs_y, (B - abs_y))

    dz = a[2] - b[2]
    abs_z = abs(dz)
    if abs_z > C:
        dz = abs_z % C
    z = min(abs_z, (C - abs_z))

    return np.sqrt(x ** 2 + y ** 2 + z ** 2)


def atoms_surround():
    global b
    global count_ti
    global total_tiO

    for k, value in enumerate(coords):

        if coords[k, 0] == "O" and coords[k, 4] == 0:

            x1 = float(coords[k, 1])
            y1 = float(coords[k, 2])
            z1 = float(coords[k, 3])
            a = [x1, y1, z1]

            if 0.1 < distance(a, b) < cutoff:
                total_tiO += 1
                 #       print(coords[k, 0], coords[k, 1], coords[k, 2], coords[k, 3])
                fw1.write("\n {}  {} {} {}".format(coords[k, 0], coords[k, 1], coords[k, 2], coords[k, 3]))

                coords[k, 4] = 1

                for j, value in enumerate(coords):

                    if coords[j, 0] == ChooseAtom and coords[j, 4] == 0:

                        x2 = float(coords[j, 1])
                        y2 = float(coords[j, 2])
                        z2 = float(coords[j, 3])
                        c = [x2, y2, z2]

                        if 0.1 < distance(a, c) < cutoff:
                            #                    print(coords[j, 0], coords[j, 1], coords[j, 2], coords[j, 3])
                            fw1.write("\n {}  {} {} {}".format(coords[j, 0], coords[j, 1], coords[j, 2], coords[j, 3]))
                            clust.append(c)
                            coords[j, 4] = 1
                            count_ti += 1
                            total_tiO += 1

                            atoms_surround()

    for s, value in enumerate(clust):
        d1 = float(clust[s][0])
        d2 = float(clust[s][1])
        d3 = float(clust[s][2])
        d = [d1, d2, d3]
        # print(d)
        b = d
        clust.remove(clust[s])
        atoms_surround()


class Trajectory:
    def __init__(self, filename, skip):
        self.filename = filename
        self.skip = skip
        self.parse_input()

    def parse_input(self):  # for getting coordinates from the trajectory file
        # global Atom_type

        global b
        global coords
        global clust
        global count_ti
        global total_tiO
        global ChooseAtom

        print()
        try:
            with open(self.filename, 'r') as fo:  # opening file
                data = fo.readlines()  # reading each line
        except FileNotFoundError as fnf_error:
            print(fnf_error)
            print()
            sys.exit()

        self.n_atoms = int(data[0].split()[0])  # reading first line as number of atoms which is expected to be integer

        self.n_steps_total = int(
            len(data) / (self.n_atoms + 2))  # we then find the total timestep by dividing the total
        # no. of lines in trajectory files by n_atoms (no. of atoms) +2 header (2 is for additional two lines (1. no.
        # of atoms 2. comment or blank line))

        self.n_steps = self.n_steps_total // self.skip  # this allows to check on every interval timestep as atoms does
        # not move appreciably in each step. Of course setting skip = 1, will give 100% accuracy.

        self.atom_list = [line.split()[0] for line in data[2:self.n_atoms + 2]]  # here i am trying to capture each
        # atoms, i will start reading from 3 rows and will go upto total no. of atoms + 2 that's why
        # data[2:self.n_atoms + 2]] also line.np.it()[0] will read the first element of the list which is indeed
        # atom name

        Atom_type = list(dict.fromkeys(self.atom_list))  ## types of atom
        print(Atom_type)
        print()
        ChooseAtom = input("Choose atom from above list : ")

        ti_average = np.zeros((200))
        # print(ti_average)

        begining_file = 0

        for step in range(self.n_steps):
            coords = np.zeros((self.n_atoms, 5), dtype=object)

            istep = step * self.skip * (self.n_atoms + 2)


            #

            if begining_file == 0:
                fw1.write("number_of_atoms \n")
                begining_file = 1
            elif begining_file == 1:
                fw1.write("\n number_of_atoms \n")

            fw1.write("")


            no_of_ti = 0

            for j, line in enumerate(data[istep + 2: istep + self.n_atoms + 2]):
                coords[j, 0:4] = [value for value in line.split()[:]]
                coords[j, 4] = 0
                if coords[j,0] == ChooseAtom :
                    no_of_ti += 1




            ncluster = 0
            total_tiO = 0
            #            print()

            for i, value in enumerate(coords):
                clust = []

                count_ti = 0

                if coords[i, 0] == ChooseAtom and coords[i, 4] == 0:
                    ncluster += 1
                    count_ti += 1
                    total_tiO += 1
                    # print("Cluster number : ", ncluster)

                    xo = float(coords[i, 1])
                    yo = float(coords[i, 2])
                    zo = float(coords[i, 3])
                    b = [xo, yo, zo]

                    # print(coords[i, 0], coords[i, 1], coords[i, 2], coords[i, 3])
                    fw1.write(" \n{}  {} {} {} ".format(coords[i, 0], coords[i, 1], coords[i, 2], coords[i, 3]))
                    coords[i, 4] = 1

                    atoms_surround()

                    # print()
                    ti_average[count_ti] += 1


        print()
        fw2.write("\n A = {} ".format(A))
        fw2.write("\n B = {} ".format(B))
        fw2.write("\n C = {} ".format(C))
        fw2.write("\n \n Cut-off = {} \n ".format(cutoff))
        fw2.write("\n \n")
        for ik in range(1, 100):
            if ti_average[ik] != 0:
                    fw2.write("{}    {}     {} \n".format(ik, round(ti_average[ik] / self.n_steps, 3),
                                                  round( (ik * ti_average[ik] * 100 ) / (self.n_steps * no_of_ti ),3) ))

                    print(ik, round(ti_average[ik] / self.n_steps, 3),
                                                  round((ik * ti_average[ik] * 100 ) / (self.n_steps * no_of_ti), 3) )

        fw2.write("\n \n \n")
        fw2.write("Total number of {} in the snapshot is {}  ".format(ChooseAtom,no_of_ti))
        fw2.write("\n \n \n")


try:
    fo1 = sys.argv[1]
    fw1 = open("cluster.xyz", 'w')
    fw2 = open("cluster.dat", 'w')
except:
    print("Please write script name followed by trajectory name (xyz)  ")
    sys.exit()


try:
    fo2 = int(sys.argv[2])
except:
    print("No steps to jump found. Taking default skip = 1")
    fo2 = 2

Trajectory(fo1,fo2)
fw1.close()
fw2.close()

import fileinput

for line in fileinput.input("cluster.xyz", inplace=True):
    # inside this loop the STDOUT will be redirected to the file
    # the comma after each print statement is needed to avoid double line breaks
    line = line.rstrip().replace("number_of_atoms", str(total_tiO))
    print(line)

print()
print("Your file [cluster.xyz]  has been created succesfully")

print()
print()
