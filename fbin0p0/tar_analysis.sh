#!/bin/bash
echo -n "Tarring analysis files... tar together on completion? [y/n]:"
if ! read -t 6 tarfiles; then
	echo "No response, I'll tar them just in case!"
	tarfiles="y"
	sleep 0.8
fi

sleep 0.8
for dir in kdum*; do
	echo "   In ${dir}"
	sleep 0.8	
	tar -czvf ${dir}_analysis.tar.gz ${dir}/analysis
	echo "         ...${dir} done."
	sleep 1.
done

if [ "$tarfiles" = "y" ]; then
	echo "tarring kdum#.tar.gz files as allkdum.tar:"
	tar -cvf allkdum.tar kdum*.tar.gz
fi

echo "Finished!"
echo "All .tar files:"
ls *.tar*
