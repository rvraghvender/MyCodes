#!/usr/bin/python3.7
#program to extract STATIS file and get the step count

f = open('STATIS','r')

while True:
    line = f.readline().split()
    if len(line) == 0 : break
    if len(line) == 3:
        x, y, z = [float(a) for a in line[0:3]]
        print(x)
