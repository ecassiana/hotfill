#! /bin/bash

clear

if [[ ! -e ./tests/out/ ]]; then
    mkdir ./tests/out/
fi

./script_gcode_solve.sh RP3 rp3

./script_gcode_solve.sh SLIC3R slic3r


./script_scanline_solve.sh SCANLINE

./script_scanline_solve.sh ALTSCANLINE


./script_hotfill_fast.sh


python3 ./summary_results_main.py False
