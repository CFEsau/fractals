#!/bin/bash

#echo ""
#echo "Combining plots..."
echo ""
sleep 0.6
outpath='../outputs/plots'

flist=''
tot=0
for file in ${outpath}/cdf_snap${snap}_*.png; do
    echo ${file}
    flist="${flist} ${file}"
    let 'tot+=1'
    #echo "flist is ${flist}"
done

convert ${flist} ${outpath}/cdf_snap${snap}.pdf
echo "Saved in cdf_snap"${snap}".pdf"
rm ${outpath}'/'cdf_snap${snap}_*.png
echo ""
sleep 0.6
