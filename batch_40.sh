#! /bin/bash

clear

LIMIT=40
for DELTA in 2.0 4.0 8.0 16.0 32.0; do
    ./script_hotfill_solve.sh ${DELTA} ${LIMIT} 
done

