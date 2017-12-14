#!/bin/bash

echo "Copying plots from all kdum to one directory..."
filestring=analysis/outputs/runinv_0100
for dir in kdum*; do
    echo "  ${dir} plots:"
    for file in ${dir}/${filestring}/*.pdf; do
#	echo "In ${filestring}"
	filename="`echo ${file} | cut -d/ -f5 | cut -d. -f1`"
#	echo "Filestring is ${filestring}, file is ${file}, filename is ${filename}"
	cp ${dir}/${filestring}/${filename}.pdf ./${filename}_${dir}.pdf
	echo "      done ${filename}.pdf"
    done
done

echo ""
echo "  Merging plots:"
pdfunite $(ls -v energies_*.pdf) allenergies.pdf
echo "        allenergies.pdf done"
pdfunite $(ls -v lagrangian_*.pdf) alllagranian.pdf
echo "        alllagrangian.pdf done"
pdfunite $(ls -v virial_*.pdf) allvirial.pdf
echo "        allvirial.pdf done"
