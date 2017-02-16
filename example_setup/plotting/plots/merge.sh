#!/bin/bash
#Merge plots for different values of k into one document.
#Read in cluster type & parameter from filenames of plots.
#Do this by stripping prefixes & suffixes from file strings:

#File format: {filepath}/k#_{cluster}_{param}.pdf

#To find {cluster} & {param} take 1st value of k and list these
#to get all {cluster} & {param} types. (This assumes {cluster} &
#{param} types are uniform across all k)

echo ""
sleep 0.6
filepath=${outpath}'/outputs/plots/'

#strip string after final '/' & before 1st '_' to get k value
for file in ${filepath}k*.pdf; do
    kval=''
    firstfile=${file##*/} #strip longest matching pattern from front end
    kval=${firstfile%%_*} #strip longest matching pattern from back end
    break 1
done

for file in ${filepath}${kval}*.pdf; do
    noext=''
    clustparam=''
    cluster=''
    param=''
    noext=${file%[.]*} #remove file extension (string from '.')
    clustparam=${noext#*[_]} #strip string before 1st '_' to remove 'k#_'
    cluster=${clustparam%[_]*} #strip {cluster}_{param} after '_'
    param=${clustparam#*[_]} #strip {cluster}_{param} before '_'
    echo ""
    echo "    Merging each k for '${cluster}', parameter '${param}'"

#Old script which used 'cut' to get the cluster type:
#Take cluster type from first file with string 'param' in:
#for clustype in ${filepath}k_*_${param}.pdf; do # format: k_cluster_param
#    cluster="`echo ${clustype} | cut -d_ -f2`" #all, FoV5pc, etc
#    break
#done

#Make a list of files to merge:
    flist=''
    for thisk in ${filepath}*_${cluster}_${param}.pdf; do
	flist="${flist} ${thisk}"
    done
    
    pdfunite ${flist} ${filepath}${param}_${cluster}.pdf
    echo "         Saved in "${param}"_"${cluster}".pdf"
    rm ${filepath}k*_${cluster}_${param}.pdf
    echo " "
    echo "         Removed ${filepath}*_${cluster}_${param}.pdf"
    sleep 0.6
done
