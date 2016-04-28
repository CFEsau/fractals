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

echo ""
echo "------------"
echo "-- Energy --"
echo "------------"
python plots/energyplots.py ${fbin} ${fdim} ${qvir}

param=E
cluster='all'
. ./plots/merge.sh


echo ""
echo "------------------"
echo "-- Virial Ratio --"
echo "------------------"
python plots/virialplots.py ${fbin} ${fdim} ${qvir}

param=Qvir
cluster='all'
. ./plots/merge.sh


echo ""
echo "----------------------"
echo "-- Half-mass radius --"
echo "----------------------"
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


echo ""
echo "------------"
echo "-- Lambda --"
echo "------------"
python plots/lambdaplots.py ${fbin} ${fdim} ${qvir}

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
plots/energy_comparek.py ${fbin} ${fdim} ${qvir}
echo ""
sleep 0.6

plots/virial_comparek.py ${fbin} ${fdim} ${qvir}
echo ""
sleep 0.6

plots/lagrange_comparek.py ${fbin} ${fdim} ${qvir}
echo ""
sleep 0.6

plots/lambda_comparek.py ${fbin} ${fdim} ${qvir}
echo ""
sleep 0.6

echo "       ... done."
