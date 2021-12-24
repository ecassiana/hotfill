#! /bin/bash

LIMIT=20
for DELTA in 2 4 8; do
    ./script_hotfill_solve.sh ${DELTA} ${LIMIT} 
done
