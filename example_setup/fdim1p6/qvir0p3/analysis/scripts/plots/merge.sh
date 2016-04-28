#!/bin/bash

#echo ""
#echo "Combining plots..."
echo ""
sleep 0.6
outpath='../outputs/plots'

flist=''
tot=0
for file in ${outpath}'/'k*_${cluster}_${param}.pdf; do
    echo ${file}
    flist="${flist} ${file}"
    let 'tot+=1'
    #echo "flist is ${flist}"
done

pdfunite ${flist} ${outpath}/${param}_${cluster}.pdf
echo "Saved in "${param}"_"${cluster}".pdf"
rm ${outpath}'/'k*_${cluster}_${param}.pdf
echo ""
sleep 0.6
