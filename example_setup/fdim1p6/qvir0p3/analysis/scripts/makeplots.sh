#!/bin/bash
#Great stack overflow thread about passing variables:
#http://stackoverflow.com/questions/9772036/pass-all-variables-from-one-shellscript-to-another

echo -n "Enter binary fraction: "
read fbin
echo -n "Enter fractal dimension: "
read fdim
echo -n "Enter virial ratio: "
read qvir

#Make 'plots' directory if it doesn't exist:
mkdir -p ../outputs/plots

echo ""
echo "**************************"
echo "*      Making plots      *"
echo "**************************"

#Plotting energies etc takes a while. Most of the time
#we're only interested in lambda. Option to skip other
#types of plots and just do lambda.

defaultskip="y"

echo -n "Make only lambda plots? [y/n]"
if ! read -t 6 skip; then
    echo "No response. Default is:"
    if [ "$defaultskip" = "y" ]; then
	skip="y"
	echo "        make only lambda plots"
    else
	skip="n"
	echo "        make all plots"
    fi
    sleep 0.8
fi

echo ""
echo "------------"
echo "-- Energy --"
echo "------------"

if [ "$skip" = "n" ]; then

    python plots/energyplots.py ${fbin} ${fdim} ${qvir}

    param=E
    cluster='all'
    . ./plots/merge.sh

else
    echo "   Skipping plots"
fi

echo ""
echo "------------------"
echo "-- Virial Ratio --"
echo "------------------"

if [ "$skip" = "n" ]; then

    python plots/virialplots.py ${fbin} ${fdim} ${qvir}

    param=Qvir
    cluster='all'
    . ./plots/merge.sh

else
    echo "   Skipping plots"
fi


echo ""
echo "----------------------"
echo "-- Half-mass radius --"
echo "----------------------"

if [ "$skip" = "n" ]; then
    python plots/lagrangeplots.py ${fbin} ${fdim} ${qvir}

    param=Rh
    cluster='all'
    . ./plots/merge.sh
    cluster='FoV5pc'
    . ./plots/merge.sh
    cluster='r2rhalf'
    . ./plots/merge.sh
    cluster='r3rhalf'
    . ./plots/merge.sh

else
    echo "   Skipping plots"
fi

echo ""
echo "------------"
echo "-- Lambda --"
echo "------------"
# (Plotting comparisons between different projections)

python plots/lambdaproj.py ${fbin} ${fdim} ${qvir}

param=lambar
cluster='all'
. ./plots/merge.sh
cluster='FoV5pc'
. ./plots/merge.sh
cluster='r2rhalf'
. ./plots/merge.sh
cluster='r3rhalf'
. ./plots/merge.sh

param=lamtil
cluster='all'
. ./plots/merge.sh
cluster='FoV5pc'
. ./plots/merge.sh
cluster='r2rhalf'
. ./plots/merge.sh
cluster='r3rhalf'
. ./plots/merge.sh

param=lamstar
cluster='all'
. ./plots/merge.sh
cluster='FoV5pc'
. ./plots/merge.sh
cluster='r2rhalf'
. ./plots/merge.sh
cluster='r3rhalf'
. ./plots/merge.sh

param=gamma
cluster='all'
. ./plots/merge.sh
cluster='FoV5pc'
. ./plots/merge.sh
cluster='r2rhalf'
. ./plots/merge.sh
cluster='r3rhalf'
. ./plots/merge.sh


echo ""
echo "       ... done."


echo ""
echo "****************************"
echo "* Comparing lambda methods *"
echo "****************************"
# (Plotting lambar, lamtilde, etc against each other)

plots/differentlambda.py ${fbin} ${fdim} ${qvir}

param=lam
cluster='all'
. ./plots/merge.sh
cluster='FoV5pc'
. ./plots/merge.sh
cluster='r2rhalf'
. ./plots/merge.sh
cluster='r3rhalf'
. ./plots/merge.sh

echo ""
echo "       ... done."


echo ""
echo "************************************"
echo "* Making comparison plots across k *"
echo "************************************"
# (Plotting all values of k together for each E, Q, lambda, etc)

if [ "$skip" = "n" ]; then
    plots/energy_comparek.py ${fbin} ${fdim} ${qvir}
    echo ""
    sleep 0.6

    plots/virial_comparek.py ${fbin} ${fdim} ${qvir}
    echo ""
    sleep 0.6

    plots/lagrange_comparek.py ${fbin} ${fdim} ${qvir}
    echo ""
    sleep 0.6
fi

plots/lambda_comparek.py ${fbin} ${fdim} ${qvir}
echo ""
sleep 0.6

echo "       ... done."
