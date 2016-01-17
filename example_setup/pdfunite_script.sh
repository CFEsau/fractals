#!/bin/bash

echo "Combining plots..."
#Read in fdim, file names for all qvir, for each plot type.
#e.g.: fdim (1.6 2.0 2.4 3.0) qvir (0.3 0.5)
#and write out as 'all

echo -n "Number of fdim values: "
read nfdim
for i in $(seq 1 $nfdim); do
    echo -n "fdim${i}: "
    read fdim_val${i}
    this_fdim=fdim_val${i}
    fdimpath_1st="`echo ${!this_fdim} | cut -d. -f1`"
    fdimpath_2nd="`echo ${!this_fdim} | cut -d. -f2`"
#    echo "fdim${fdimpath_1st}p${fdimpath_2nd}"
    fdim_path[${i}]="fdim${fdimpath_1st}p${fdimpath_2nd}"
    echo "${this_fdim} is  ${!this_fdim}"
    echo "fdim_path${i} is ${fdim_path[i]}"
done

echo ""
echo -n "Number of qvir values: "
read nqvir
for j in $(seq 1 $nqvir); do

    echo -n "nqvir${j}: "
    read qvir_val${j}
    this_qvir=qvir_val${j}
    qvirpath_1st="`echo ${!this_qvir} | cut -d. -f1`"
    qvirpath_2nd="`echo ${!this_qvir} | cut -d. -f2`"
    qvir_path[${j}]="qvir${qvirpath_1st}p${qvirpath_2nd}"
    echo "${this_qvir} is  ${!this_qvir}"
    echo "qvir_path${j} is ${qvir_path[j]}"
done

echo ""
echo "Enter type of plot"
echo -n "(1 for energies, 2 for lagrangian, 3 for virial): "
read plottype
if plottype==1; then
    plotstring='energies'
elif plottype==2; then
    plotstring='lagrangian'
elif plottype==3; then
    plotstring='virial'
fi

list=""
tot=0
echo ""
echo "Combining the following files:"
for i in $(seq 1 $nfdim); do
    for j in $(seq 1 $nqvir); do
	#echo "i: ${i} j: ${j}"
	filepath="${fdim_path[i]}/compare${plotstring}_${qvir_path[j]}.pdf"
	echo "${filepath}"
	#echo "filepath: ${fdim_path[i]}/compare${plotstring}_${qvir_path[j]}.pdf"
#	echo "${!filename}"
#	eval "echo \$$filename"
#see http://unix.stackexchange.com/questions/41406/use-a-variable-reference-inside-another-variable for eval

	list="${list} ${filepath}"
	tot+=1
    done
done

echo "Total: ${tot}"

echo "Combined in all${plotstring}.pdf"

echo ""
echo "pdfunite ${list} all${plotstring}.pdf"

#pdfunite fdim1p6/compareenergies_qvir0p3.pdf fdim1p6/compareenergies_qvir0p5.pdf fdim2p0/compareenergies_qvir0p3.pdf fdim2p0/compareenergies_qvir0p5.pdf fdim2p4/compareenergies_qvir0p3.pdf fdim2p4/compareenergies_qvir0p5.pdf fdim3p0/compareenergies_qvir0p3.pdf fdim3p0/compareenergies_qvir0p5.pdf allenergies_comparefdim.pdf

#pdfunite fdim1p6/comparelagrangian_qvir0p3.pdf fdim1p6/comparelagrangian_qvir0p5.pdf fdim2p0/comparelagrangian_qvir0p3.pdf fdim2p0/comparelagrangian_qvir0p5.pdf fdim2p4/comparelagrangian_qvir0p3.pdf fdim2p4/comparelagrangian_qvir0p5.pdf fdim3p0/comparelagrangian_qvir0p3.pdf fdim3p0/comparelagrangian_qvir0p5.pdf alllagrangian_comparefdim.pdf

#pdfunite fdim1p6/comparevirial_qvir0p3.pdf fdim1p6/comparevirial_qvir0p5.pdf fdim2p0/comparevirial_qvir0p3.pdf fdim2p0/comparevirial_qvir0p5.pdf fdim2p4/comparevirial_qvir0p3.pdf fdim2p4/comparevirial_qvir0p5.pdf fdim3p0/comparevirial_qvir0p3.pdf fdim3p0/comparevirial_qvir0p5.pdf allvirial_comparefdim.pdf
