#!/bin/bash

#Filepath: ../outputs/runinv_k##/cluster_[type]/escaped[proj]

outdir=("../outputs")
    printf "\n"

for kdir in ${outdir}"/"runinv_k*; do
    #printf " \nIn ${kdir}\n"

    #dont do cluster_all since no 'escaped' files here

    for cdir in ${kdir}/cluster_FoV*; do
	#printf "   \nIn ${cdir}\n"
	tar -czvf ${cdir}/escaped.tar.gz ${cdir}/escaped*.txt
	rm ${cdir}/escaped*.txt
	printf "\n"
    done

    for cdir in ${kdir}/cluster_r*; do
	#printf "   \nIn ${cdir}\n"
	tar -czvf ${cdir}/escaped.tar.gz ${cdir}/escaped*.txt
	rm ${cdir}/escaped*.txt
	printf "\n"
    done

    sleep 0.8

done
