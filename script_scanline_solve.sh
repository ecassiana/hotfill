#! /bin/bash

ALG=$1

for FILENAME in $(ls ./tests/in/); do
	DIRECTORY=./tests/out/${FILENAME%.*}_${ALG}
	OUTFILE=${FILENAME%.*}_${ALG}/out
	if [ ! -d $DIRECTORY ]; then	
		mkdir -p ${DIRECTORY}/
		touch ${DIRECTORY}/out.txt
		echo -n "Running: ${FILENAME} | Algorithm: ${ALG} | Start Time: "
		date +"%H:%M:%S"
		timeout 30m python3 main_sc.py ${FILENAME} ${ALG} ${OUTFILE} 
		echo -n "!!! Done | End Time: "
		date +"%H:%M:%S"
	
	else
		echo "Skipping: ${FILENAME} | Algorithm: ${ALG}"
	
	fi
done

echo "!!!!!!! END | Algorithm: ${ALG}"
