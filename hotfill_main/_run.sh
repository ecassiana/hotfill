#! /bin/bash

clear

if [[ ! -e ./tests/out/ ]]; then
    mkdir ./tests/out/
fi

./_script_gcode_solve.sh RP3 rp3

./_script_gcode_solve.sh SLIC3R slic3r


./_script_scanline_solve.sh SCANLINE

./_script_scanline_solve.sh ALTSCANLINE


./_script_hotfill.sh


python3 ./summary_results_main.py False
