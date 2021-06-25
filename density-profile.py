#!/usr/bin/python3.7
#program to calculate the partial distribution funtion and plot it                    ##
import platform                                                                       ##
import sys                                                                            ##
print()
print()
print("System Information : ",sys.version, sys.platform)                              ##
print("Build Info : ",platform.machine(),platform.system())                           ##
print()                                                                               ##
print()                                                                               ##
                                                                                      #                                                   ##
print()							                              ##
print("--------------------------------------------------------------------------")   ##
                                                                                      ##
print()                                                                               ##
print()                                                                               ##
                                                                                      ##
########################################################################################




import matplotlib.pyplot as plt
import numpy as np


A =   float(input("Enter lattice parameter A : ")) #Entering cell parameters A,B,C
B =   float(input("Enter lattice parameter B : "))
C =   float(input("Enter lattice parameter C : "))

print()


class Trajectory:
    def __init__(self,filename, skip = 1):
        self.filename = filename
        self.skip = skip
        self.parse_input()

    def parse_input(self):  # for getting coordinates from the trajectory file
        global Atom_type
        print()
        try:
            with open(self.filename, 'r') as fo: #opening file
                data = fo.readlines()   # reading each line
        except FileNotFoundError as fnf_error:
            print(fnf_error)
            print()
            sys.exit()

        self.n_atoms = int(data[0].split()[0]) # reading first line as number of atoms which is expected to be integer
        self.n_steps_total = int(len(data) / (self.n_atoms +2))  # we then find the total timestep by dividing the total
        # no. of lines in trajectory files by n_atoms (no. of atoms) +2 header (2 is for additional two lines (1. no.
        # of atoms 2. comment or blank line))


        self.n_steps = self.n_steps_total // self.skip # this allows to check on every interval timestep as atoms does
        #not move appreciably in each step. Of course setting skip = 1, will give 100% accuracy.


        self.atom_list = [line.split()[0] for line in data[2:self.n_atoms + 2]]  # here i am trying to capture each atoms
        # i will start reading from 3 rows and will go upto total no. of atoms + 2 that's why data[2:self.n_atoms + 2]]
        # also line.np.it()[0] will read the first element of the list which is indeed atom name

        Atom_type = list(dict.fromkeys(self.atom_list))  ## types of atom
        print()

        print("Available atom/s in the given Trajectory file is/are :", Atom_type)
        print()


        self.Choosen_atom = input("Enter one of the atom from the above list to find density profile : ")
        print()


        if self.Choosen_atom not in Atom_type:
            raise Exception(self.Choosen_atom, "is not in the above list.")
            sys.exit()


        self.n_bin = int(input("Enter the number of parts to divide the axis : "))
        print()

        self.xaxis = np.arange(self.n_bin + 1)
        self.delta = float(A / self.n_bin)

        self.population_x = np.zeros((self.n_bin+1))
        self.population_y = np.zeros((self.n_bin+1))
        self.population_z = np.zeros((self.n_bin+1))
        # print (self.population[0])
        # print()

        for step in range(self.n_steps): # here i am creating a loop with nsteps
            coords = np.zeros((self.n_atoms, 4), dtype = object) # for each step i am creating zero matrix with n_atoms rows and 3 colu.
            i = step * self.skip * (self.n_atoms + 2)  # now this i will iterate each atom index in increasing timestep.
            # and thus this i represents move forward in timestep starting from 0.

            self.atom_count = 0

            for j, line in enumerate(data[i + 2 : i + self.n_atoms + 2]): # here i am enumerating in each current
                # timestep "i". Thus this loop scans whole atoms in a single timestep.

                coords[j, :] = [value for value in line.split()[:]]# here data is a list where number of
                # element is equals to no. of lines in trajectory file. Also, line is temporary variable name
                # assigned to each element of data. And this line is alias for atom name, x ,y and z of each atom.
                #Therefore the list has started from 1 i.e. [1:] so that it ignores the atom name.


                if  coords[j,0] == self.Choosen_atom:
                    self.atom_count += 1
                    A0 = float(0)
                    AN = A0
                    n = 0
                    is_between = False
                    if float(coords[j,1]) < 0 : coords[j,1] = float(coords[j,1]) + A
                    if float(coords[j, 1]) > A: coords[j, 1] = float(coords[j, 1]) - A
                    while is_between == False:
                        n += 1
                        AI = A0 +  n  * self.delta
                        # print(float(coords[j,1]),n*self.delta,n)
                        is_between = AN < float(coords[j,1]) <= AI
                        # print(AI, AN, is_between, n)
                        AN = AI

                    # print()
                    self.population_x[n] += 1




                    A0 = float(0)
                    AN = A0
                    n = 0
                    is_between = False
                    if float(coords[j, 2]) < 0: coords[j, 2] = float(coords[j, 2]) + B
                    if float(coords[j, 2]) > B: coords[j, 2] = float(coords[j, 2]) - B
                    while is_between == False:
                        n += 1
                        AI = A0 + n * self.delta
                        # print(float(coords[j,1]))
                        is_between = AN < float(coords[j, 2]) <= AI
                        # print(AI, AN, is_between, n)
                        AN = AI

                    # print()
                    self.population_y[n] += 1

                    A0 = float(0)
                    AN = A0
                    n = 0
                    is_between = False
                    if float(coords[j, 3]) < 0: coords[j, 3] = float(coords[j, 3]) + C
                    if float(coords[j, 3]) > C: coords[j, 3] = float(coords[j, 3]) - C
                    while is_between == False:
                        n += 1
                        AI = A0 + n * self.delta
                        # print(float(coords[j,1]))
                        is_between = AN < float(coords[j, 3]) <= AI
                        # print(AI, AN, is_between, n)
                        AN = AI

                    # print()
                    self.population_z[n] += 1


                # print(n*delta,self.population[n])
        print("Total no. of {} atom are : {}".format(self.Choosen_atom, self.atom_count))
        # for n in range(11):
             # print(n,n*self.delta,self.population_z[n]/self.n_steps)

    def plot(self, filename=""):
        fig, (ax1, ax2, ax3) = plt.subplots(nrows=3, ncols=1, sharey=True, sharex=True)
        fig.suptitle('Density Profiles \n Number of Bins : {} '.format(self.n_bin))
        fig.subplots_adjust(hspace=0.5)

        plt.MaxNLocator(10)

        ax1.bar(self.xaxis * self.delta, self.population_x / self.n_steps, edgecolor='black', align='edge', width=-self.delta * 0.9)
        ax1.hlines(self.atom_count/ self.n_bin , 0, A, linestyles='dashed', colors='red')


        ax2.bar(self.xaxis * self.delta, self.population_y / self.n_steps, edgecolor='black', align='edge', width=-self.delta * 0.9)
        ax2.hlines(self.atom_count / self.n_bin, 0, B, linestyles='dashed', colors='red')

        ax3.bar(self.xaxis * self.delta, self.population_z / self.n_steps, edgecolor='black', align='edge', width=-self.delta * 0.9)
        ax3.hlines(self.atom_count / self.n_bin, 0, C, linestyles='dashed', colors='red')


        ax1.set_xticks(np.arange(0, int((self.n_bin + 1) * self.delta), 1))
        ax2.set_xticks(np.arange(0, int((self.n_bin + 1) * self.delta), 1))
        ax3.set_xticks(np.arange(0, int((self.n_bin + 1) * self.delta), 1))

        ax1.set_xlabel('x (A$^0$)')
        ax2.set_xlabel('y (A$^0$)')
        ax3.set_xlabel('z (A$^0$)')


        # fig.legend('Number of Bins : {}'.format(self.n_bin))

        print(self.xaxis*self.delta,self.population_x/self.n_steps)

        fig.text(0.04, 0.5, 'No. of {} atoms'.format(self.Choosen_atom), va='center', rotation='vertical')

        ax1.grid(True)
        ax2.grid(True)
        ax3.grid(True)

        plt.savefig('density-profile.png', dpi=300, bbox='tight', format='png')
        plt.clf()


try:
    fo1 = sys.argv[1]
except:
    print("Please enter script name, File name, Steps to Jump (optional: Default = 1) in the command line")
    print()
    sys.exit()

try:
    skip_ts = int(sys.argv[2])
except:
    skip_ts = 1



file = Trajectory(fo1,skip_ts)
file.plot(fo1)

print()



print()
print("Your file [ density-profile.png ] has been created succesfully")
