#!/bin/bash
LANG=en_US

######  PARAMETERS OF CONTROL FILE    ##########

RelaxTolerance=2
Count=1


RlxtolTest=$(grep rlxtol CONTROL)
sed -e  "s/$RlxtolTest/rlxtol           $RelaxTolerance/g" CONTROL > CONTROL.tmp1
cp CONTROL.tmp1 CONTROL

IterativeStep=$(grep steps CONTROL| awk '{print $2}')

while [[ $IterativeStep -le 100000 ]]
do
  mpirun -np 8 DLPOLY.Z
  if [ $? -eq 0 ]
  then
    echo "The script ran ok"
    mv REVIVE REVOLD
    mv REVCON CONFIG
    StepTest=$(grep steps CONTROL)
    StepNumber=$(grep steps CONTROL | awk '{print $2}')
    UpdatedStepNumber=$(echo "$StepNumber + 1000" | bc -l)
    sed -e "s/$StepTest/steps        $UpdatedStepNumber"/ CONTROL > CONTROL.tmp3
    cp CONTROL.tmp3 CONTROL
    IterativeStep=$UpdatedStepNumber
    RelaxtolTest=$(grep rlxtol CONTROL)
    sed -e "s/$RelaxtolTest/rlxtol           ${RelaxTolerance}/" CONTROL > CONTROL.tmp2
    cp CONTROL.tmp2 CONTROL
    continue
    exit 0
  else
    echo "The script failed" >&2
    grep "DL_POLY_4 terminated due to error" OUTPUT
    ErrorCode=$(grep "DL_POLY_4 terminated due to error" OUTPUT | awk '{print $6}')
    if [[ $ErrorCode == 474 ]]
    then
      CurrentRlxtol=$(grep "specified convergence tolerance:" OUTPUT | awk '{print $5}')  
      NeededRlxtol=$(grep "specified convergence tolerance:" OUTPUT | awk '{print $12}')
      UpdatedRlxtol=$(echo "$NeededRlxtol + 2" | bc -l)
      CurrentRlxtol=$UpdatedRlxtol
      RelaxtolTest=$(grep rlxtol CONTROL)
      sed -e  "s/$RelaxtolTest/rlxtol           ${UpdatedRlxtol}/" CONTROL > CONTROL.tmp2
      cp CONTROL.tmp2 CONTROL
      cat CONTROL
      if [ -f REVIVE ]; then
        mv CONFIG CONFIG-$Count
        mv REVOLD REVOLD-$Count
        mv REVIVE REVOLD
        mv REVCON CONFIG
        Count=$((Count+1))
      fi
      sleep 2
      continue
    else
      break
    fi
    exit 1
  fi
  break
done
