#!/bin/bash

#echo ""
#echo "Combining plots..."
echo ""
sleep 0.6

flist=''
for file in ${outpath}/cdf_snap${snap}_*.png; do
    echo ${file}
    flist="${flist} ${file}"
    #echo "flist is ${flist}"
done

convert ${flist} ${outpath}/cdf_snap${snap}.pdf
echo ""
echo "Saved in cdf_snap"${snap}".pdf"
rm ${cdfoutpath}cdf_snap${snap}_*.png
echo ""
sleep 0.6
