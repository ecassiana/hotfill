#! /bin/bash

DELTA=$1
LIMIT=$2

for FILENAME in $(ls ./tests/in/); do
	DIRECTORY=./tests/out/${FILENAME%.*}_delta${DELTA}_maxband${LIMIT}
	OUTFILE=${FILENAME%.*}_delta${DELTA}_maxband${LIMIT}/out
	if [ ! -d $DIRECTORY ]; then	
		mkdir -p ${DIRECTORY}/
		touch ${DIRECTORY}/out.txt
		echo -n "Running: ${FILENAME} | Delta: ${DELTA} | Limit Time: ${LIMIT} | Start Time: "
		date +"%H:%M:%S"
		timeout 4h python3 main.py ${FILENAME} ${DELTA} ${LIMIT} ${OUTFILE} 
		echo -n "!!! Done | End Time: "
		date +"%H:%M:%S"
	
	else
		echo "Skipping: ${FILENAME} | Delta: ${DELTA} | Limit Time: ${LIMIT}"
	
	fi
done

echo "!!!!!!! END | Delta: ${DELTA} | Limit Time: ${LIMIT}"
