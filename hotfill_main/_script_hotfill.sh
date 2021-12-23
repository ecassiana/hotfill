#! /bin/bash

DELTA=1.4
for LIMIT in 5 15 20 25 30 80; do
    ./_script_hotfill_solve.sh ${DELTA} ${LIMIT} 
done

DELTA=2.0
for LIMIT in 5 15 20 25 30 80; do
    ./_script_hotfill_solve.sh ${DELTA} ${LIMIT} 
done

DELTA=2.3
for LIMIT in 20; do
    ./_script_hotfill_solve.sh ${DELTA} ${LIMIT} 
done

DELTA=2.8
for LIMIT in 5 15 20 25 30 80; do
    ./_script_hotfill_solve.sh ${DELTA} ${LIMIT} 
done

DELTA=4.0
for LIMIT in 5 15 20 25 30 80; do
    ./_script_hotfill_solve.sh ${DELTA} ${LIMIT} 
done

DELTA=5.7
for LIMIT in 5 15 20 25 30 80; do
    ./_script_hotfill_solve.sh ${DELTA} ${LIMIT} 
done

DELTA=8.0
for LIMIT in 5 10 15 20 25 30 40 50 60 70 80; do
    ./_script_hotfill_solve.sh ${DELTA} ${LIMIT} 
done

DELTA=11.3
for LIMIT in 5 15 20 25 30 80; do
    ./_script_hotfill_solve.sh ${DELTA} ${LIMIT} 
done

DELTA=16.0
for LIMIT in 5 10 15 20 25 30 40 50 60 70 80; do
    ./_script_hotfill_solve.sh ${DELTA} ${LIMIT} 
done

DELTA=22.6
for LIMIT in 5 15 20 25 30 80; do
    ./_script_hotfill_solve.sh ${DELTA} ${LIMIT} 
done

DELTA=32.0
for LIMIT in 5 15 20 25 30 80; do
    ./_script_hotfill_solve.sh ${DELTA} ${LIMIT} 
done

DELTA=45.3
for LIMIT in 5 15 20 25 30 80; do
    ./_script_hotfill_solve.sh ${DELTA} ${LIMIT} 
done

DELTA=64.0
for LIMIT in 5 10 15 20 25 30 40 50 60 70 80; do
    ./_script_hotfill_solve.sh ${DELTA} ${LIMIT} 
done

DELTA=90.5
for LIMIT in 5 15 20 25 30 80; do
    ./_script_hotfill_solve.sh ${DELTA} ${LIMIT} 
done

DELTA=128.0
for LIMIT in 5 15 20 25 30 80; do
    ./_script_hotfill_solve.sh ${DELTA} ${LIMIT} 
done

DELTA=181.0
for LIMIT in 5 15 20 25 30 80; do
    ./_script_hotfill_solve.sh ${DELTA} ${LIMIT} 
done

DELTA=256.0
for LIMIT in 5 15 20 25 30 80; do
    ./_script_hotfill_solve.sh ${DELTA} ${LIMIT} 
done