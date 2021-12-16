#! /bin/bash

ALG=$1
CASE=$2

for FILENAME in $(ls ./tests/in_gcodes/); do
	if [[ ${FILENAME} =  *${CASE}* ]]; then		
		DIRECTORY=./tests/out/${FILENAME%.*}_${ALG}
		OUTFILE=${FILENAME%.*}_${ALG}/out
		if [ ! -d $DIRECTORY ]; then
			mkdir -p ${DIRECTORY}/
			touch ${DIRECTORY}/out.txt
			echo -n "Running: ${FILENAME} | Algorithm: ${ALG} | Start Time: "
			date +"%H:%M:%S"
			timeout 30m python3 main_gcode.py ${FILENAME} ${OUTFILE} 
			echo -n "!!! Done | End Time: "
			date +"%H:%M:%S"
		
		else
			echo "Skipping: ${FILENAME} | Algorithm: ${ALG}"
		
		fi
	fi
done

echo "!!!!!!! END | Algorithm: ${ALG}"
