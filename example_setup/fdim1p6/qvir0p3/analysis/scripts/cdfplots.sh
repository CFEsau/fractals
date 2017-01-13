#!/bin/bash
#Great stack overflow thread about passing variables:
#http://stackoverflow.com/questions/9772036/pass-all-variables-from-one-shellscript-to-another

#TODO: change this to read info from file path
#echo -n "Enter binary fraction: "
#read fbin
#echo -n "Enter fractal dimension: "
#read fdim
#echo -n "Enter virial ratio: "
#read qvir

#Make 'plots' directory if it doesn't exist:
mkdir -p ../outputs/plots
snapnum=(50 76 105 125) #snapshot number
nmst=10 #number of stars in MST
simnum=07 #simulation number
nran=20 #number of random MSTs

echo ""
echo "*************************"
echo "*     Plotting CDFs     *"
echo "*************************"

for snap in ${snapnum[@]}; do

    echo ""
    echo "==================================="
    echo "  Snapshot "${snap}":"

    echo ""
    echo "-----------------"
    echo "-- Object CDFs --"
    echo "-----------------"
    # (Plotting CDFs of 'object' star edge lengths)
    plots/cdf.py ${snap} ${nmst} ${simnum}
    
    
    echo ""
    echo "-----------------"
    echo "-- Random CDFs --"
    echo "-----------------"
    # (Plotting CDFs of random star edge lengths)
    plots/cdf_ran.py ${snap} ${nmst} ${simnum} ${nran}
    
    
    echo ""
    echo "------------------"
    echo "-- Merging CDFs --"
    echo "------------------"
    
    . ./plots/mergecdfs.sh
    
    echo ""
    echo "       ... done."
    
    
    echo ""
    echo "-----------------------"
    echo "-- Plotting D1 vs D2 --"
    echo "-----------------------"
    # Difference between x values at given points in object & random CDFs
    
    plots/cdf_delta.py ${snap} ${nmst} ${simnum} ${nran}
    
    
    echo ""
    echo "       ... done."
    
done
