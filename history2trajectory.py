#!/usr/bin/python3.7
#program to generate the trajectory file out of HISTORY file in the present directory

fo = open('HISTORY','r')

fw = open('Trajectory.xyz','w')

cs = input("Press ' y ' to include both core-shell or  ' n ' for only cores : ")

if cs == 'y':
    denominator = 1   #if it is already a coreshell way written history files then its okay
elif cs == 'n':
    denominator = 2   # since we want no coreshell 2 will divide the total no. of atoms

count = 2    # any value other than 0 will work
while True:
    line = fo.readline().split()
    if len(line) == 0 : break  # so that we get a break from reading file in end
    if line[0] == 'timestep' :  # this will look for timestep line and will print after getting
        fw.write(" %i \n \n " % (int(line[2])/int(denominator))) # no. of atoms picked from history file
    if line[0] == 'O':  # the first atom
        count = 0
        fw.write('  O        ')
        continue
    if cs == 'y':  # this wont allow to write the coordinates of shell if no core shell is used
        if line[0] == 'O_sh':
            count = 0
            fw.write('  O_sh        ')
            continue
    if line[0] == 'Te':
        count = 0
        fw.write('  Te        ')
        continue
    if cs == 'y':  # same as above
        if line[0] == 'Te_sh':
            count = 0
            fw.write('  Te_sh        ')
            continue
    if count == 0:
        fw.write( " %8.9f  %8.9f %8.9f \n " % (float(line[0]),float(line[1]),float(line[2])) )
        count = 1
print()
print("Conversion into .xyz format with the file name  ( Trajectory.xyz ) has been done.")

fw.close()
fo.close()
