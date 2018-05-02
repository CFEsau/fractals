#!/bin/bash

#fvals=(1 2 3 4)
#qvals=(1 2)

#printf "%s\n" "${fvals[@]}"

for f in {3..4}; do
    for q in {1..2}; do
	#for k in {1..10}; do
	    
	    echo "Rscript 2Dpositions_sciencey_mst_test2.r $f $q $k"
	    Rscript 2Dpositions_sciencey_mst_test2.r $f $q #$k
	    
	#done
    done
done
