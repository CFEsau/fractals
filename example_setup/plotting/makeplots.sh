#!/bin/bash
#Great stack overflow thread about passing variables:
#http://stackoverflow.com/questions/9772036/pass-all-variables-from-one-shellscript-to-another

echo -n "Enter binary fraction: "
while read fbin; do
    if [ "$fbin" = "0%" ]; then fbindir="fbin0p0"; break
    elif [ "$fbin" = "50%" ]; then fbindir="fbin0p5"; break
    elif [ "$fbin" = "100%" ]; then fbindir="fbin1p0"; break
    else
	echo "     Warning: Binary fraction '$fbin' not recognised."
	echo "     Options are 0%, 50%, 100%."
	echo -n "     Enter new binary fraction: "
    fi
done

echo -n "Enter fractal dimension: "
while read fdim; do
    if [ "$fdim" = "1.6" ]; then fdimdir="fdim1p6"; break
    elif [ "$fdim" = "2.0" ]; then fdimdir="fdim2p0"; break
    elif [ "$fdim" = "2.6" ]; then fdimdir="fdim2p6"; break
    elif [ "$fdim" = "3.0" ]; then fdimdir="fdim3p0"; break
    else
	echo "     Warning: Fractal dimension '$fdim' not recognised."
	echo "     Options are 1.6, 2.0, 2.6, 3.0."
	echo -n "     Enter new fractal dimension: "
    fi
done

echo -n "Enter virial ratio: "
while read qvir; do
    if [ "$qvir" = "0.3" ]; then qvirdir="qvir0p3"; break
    elif [ "$qvir" = "0.5" ]; then qvirdir="qvir0p5"; break
    else
	echo "     Warning: Virial ratio '$qvir' not recognised."
	echo "     Options are 0.3, 0.5."
	echo -n "     Enter new virial ratio: "
    fi
done

outpath="../${fbindir}/${fdimdir}/${qvirdir}"
echo "Model path is $outpath"

#Make 'plots' directory if it doesn't exist:
mkdir -p ${outpath}/outputs/plots

echo ""
echo "**************************"
echo "*      Making plots      *"
echo "**************************"

#Plotting energies etc takes a while. Most of the time
#we're only interested in lambda. Option to skip other
#types of plots and just do lambda.

defaultskip="y"

echo -n "Skip to lambda plots? [y/n] "
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
    python plots/energyplots.py ${fbin} ${fdim} ${qvir} ${outpath}
    
    #param=E
    . ./plots/merge.sh

else
    echo "   Skipping plots"
fi


echo ""
echo "------------------"
echo "-- Virial Ratio --"
echo "------------------"

if [ "$skip" = "n" ]; then
    python plots/virialplots.py ${fbin} ${fdim} ${qvir} ${outpath}
    
    #param=Qvir
    . ./plots/merge.sh

else
    echo "   Skipping plots"
fi


echo ""
echo "----------------------"
echo "-- Half-mass radius --"
echo "----------------------"

if [ "$skip" = "n" ]; then
    python plots/lagrangeplots.py ${fbin} ${fdim} ${qvir} ${outpath}
    
    #param=Rh
    . ./plots/merge.sh

else
    echo "   Skipping plots"
fi


echo ""
echo "------------"
echo "-- Lambda --"
echo "------------"
# Plot comparisons between different projections for each lambda

python plots/lambdaproj.py ${fbin} ${fdim} ${qvir} ${outpath}

. ./plots/merge.sh

echo ""
echo "       ... done."


echo ""
echo "****************************"
echo "* Comparing lambda methods *"
echo "****************************"
# Plot lambar, lamtilde, etc against each other for each k

plots/differentlambda.py ${fbin} ${fdim} ${qvir} ${outpath}

#param=alllam
. ./plots/merge.sh

echo ""
echo "       ... done."


echo ""
echo "************************************"
echo "* Making comparison plots across k *"
echo "************************************"
# Plot all values of k together for each E, Q, lambda, etc

if [ "$skip" = "n" ]; then
    plots/energy_comparek.py ${fbin} ${fdim} ${qvir} ${outpath}
    echo ""
    sleep 0.6
    
    plots/virial_comparek.py ${fbin} ${fdim} ${qvir} ${outpath}
    echo ""
    sleep 0.6
    
    plots/lagrange_comparek.py ${fbin} ${fdim} ${qvir} ${outpath}
    echo ""
    sleep 0.6
fi

plots/lambda_comparek.py ${fbin} ${fdim} ${qvir} ${outpath}
echo ""
sleep 0.6

echo "       ... done."
