#! /bin/bash

clear

./script_hotfill.sh


./script_scanline_solve.sh SCANLINE

./script_scanline_solve.sh ALTSCANLINE


./script_gcode_solve.sh RP3 rp3

./script_gcode_solve.sh SLIC3R slic3r


./script_select_models.sh


python3 ./summary_results_main.py
