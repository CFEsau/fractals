#!/bin/bash

echo -n "Enter binary fraction: "
read fbin
echo -n "Enter fractal dimension: "
read fdim
echo -n "Enter virial ratio: "
read qvir

echo ""
echo "**************************"
echo "*      Making plots      *"
echo "**************************"
python plots/energyplots.py ${fbin} ${fdim} ${qvir}
python plots/virialplots.py ${fbin} ${fdim} ${qvir}
python plots/lagrangeplots.py ${fbin} ${fdim} ${qvir}
python plots/lambdaplots.py ${fbin} ${fdim} ${qvir}
echo ""
echo "       ... done."


echo ""
echo "***************************"
echo "* Making comparison plots *"
echo "***************************"
plots/compareenergy.py ${fbin} ${fdim} ${qvir}
echo ""
sleep 0.6
plots/comparevirial.py ${fbin} ${fdim} ${qvir}
echo ""
sleep 0.6
plots/comparelagrange.py ${fbin} ${fdim} ${qvir}
echo ""
sleep 0.6
plots/comparelambda.py ${fbin} ${fdim} ${qvir}
echo ""
sleep 0.6
echo "       ... done."


echo ""
echo "****************************"
echo "* Comparing lambda methods *"
echo "****************************"
plots/differentlambda.py ${fbin} ${fdim} ${qvir}
echo ""
echo "       ... done."


echo ""
echo "**************************"
echo "*     Merging plots      *"
echo "**************************"
plots/merge.sh
echo ""
echo "       ... done."
