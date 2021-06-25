#!/bin/bash

dl_analyser
#System: MD_time volume  pressure
#Energy: MD_time engcns  engsrc
#Temperature: MD_time temp

grep 'System:' statis.out > VolPress.dat0
sed 's/System://g' VolPress.dat0 > VolPress.dat1
sed 's/MD_time/#MD_time/g' VolPress.dat1 > VolPress.dat
rm -f VolPress.dat0 VolPress.dat1

grep 'Energy:' statis.out > Energies.dat0
sed 's/Energy://g' Energies.dat0 > Energies.dat1
sed 's/MD_time/#MD_time/g' Energies.dat1 > Energies.dat 
rm -f Energies.dat0 Energies.dat1 

grep 'Temperature:' statis.out > Temp.dat0
sed 's/Temperature://g' Temp.dat0 > Temp.dat1
sed 's/MD_time/#MD_time/g' Temp.dat1 > Temp.dat
rm -f Temp.dat0 Temp.dat1 

#grep 'A=' dl_analyser.output > Cell.dat0
#sed 's/A=//g' Cell.dat0 > Cell.dat1 
#sed 's/B=//g' Cell.dat1 > Cell.dat2 
#sed 's/C=//g' Cell.dat2 > Cell.dat3
#sed 's/Alpha=//g' Cell.dat3 > Cell.dat4
#sed 's/Beta=//g' Cell.dat4 > Cell.dat
#rm -f Cell.dat1 Cell.dat2 Cell.dat3 Cell.dat4 Cell.dat0 


