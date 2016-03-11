#!/bin/bash

echo ""
echo "Combining plots..."
echo ""
sleep 0.6
outpath='../outputs'


flist=''
tot=0
for file in ${outpath}'/'k*energies.pdf; do
    echo ${file}
    flist="${flist} ${file}"
    let 'tot+=1'
    #echo "flist is ${flist}"
done

echo "Saved in all_energies.pdf"
pdfunite ${flist} ${outpath}/all_energies.pdf
echo ""
sleep 0.6


flist=''
tot=0
for file in ${outpath}'/'k*lagrangian.pdf; do
    echo ${file}
    flist="${flist} ${file}"
    let 'tot+=1'
    #echo "flist is ${flist}"
done

echo "Saved in all_lagrangian.pdf"
pdfunite ${flist} ${outpath}/all_lagrangian.pdf
echo ""
sleep 0.6


flist=''
tot=0
echo ""
for file in ${outpath}'/'k*virial.pdf; do
    echo ${file}
    flist="${flist} ${file}"
    let 'tot+=1'
    #echo "flist is ${flist}"
done

echo "Saved in all_virial.pdf"
pdfunite ${flist} ${outpath}/all_virial.pdf
echo ""
sleep 0.6


flist=''
tot=0
for file in ${outpath}'/'k*lambda.pdf; do
    echo ${file}
    flist="${flist} ${file}"
    let 'tot+=1'
    #echo "flist is ${flist}"
done

echo "Saved in all_lambda.pdf"
pdfunite ${flist} ${outpath}/all_lambda.pdf
echo ""
sleep 0.6


flist=''
tot=0
for file in ${outpath}'/'k*lambdabar.pdf; do
    echo ${file}
    flist="${flist} ${file}"
    let 'tot+=1'
    #echo "flist is ${flist}"
done

echo "Saved in all_lambdabar.pdf"
pdfunite ${flist} ${outpath}/all_lambdabar.pdf
echo ""
sleep 0.6


flist=''
tot=0
for file in ${outpath}'/'k*lambdatilde.pdf; do
    echo ${file}
    flist="${flist} ${file}"
    let 'tot+=1'
    #echo "flist is ${flist}"
done

echo "Saved in all_lambdatilde.pdf"
pdfunite ${flist} ${outpath}/all_lambdatilde.pdf
echo ""
sleep 0.6


flist=''
tot=0
for file in ${outpath}'/'k*lambdastar.pdf; do
    echo ${file}
    flist="${flist} ${file}"
    let 'tot+=1'
    #echo "flist is ${flist}"
done

echo "Saved in all_lambdastar.pdf"
pdfunite ${flist} ${outpath}/all_lambdastar.pdf
echo ""
sleep 0.6


flist=''
tot=0
for file in ${outpath}'/'k*gamma.pdf; do
    echo ${file}
    flist="${flist} ${file}"
    let 'tot+=1'
    #echo "flist is ${flist}"
done

echo "Saved in all_gamma.pdf"
pdfunite ${flist} ${outpath}/all_gamma.pdf
echo ""
sleep 0.6
