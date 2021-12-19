#! /bin/bash
# Last edited on 2021-06-07 21:25:32 by jstolfi

# Reads a data file produced by the test program of module {move_parms}
# with tag "times", containing extrusion and jump times as function of distance.
# Outputs a PNG file to standard output.

cmd="$0"; cmd="${cmd##*/}"

dfile="$1"; shift  # Data file name.

# END COMMAND LINE PARSING
# ----------------------------------------------------------------------

show=1

# Prefix for temporary file names
tmp="/tmp/$$"

hPix=1200
vPix=1200

tfile="${tmp}.png"
export GDFONTPATH="tt-fonts"

gnuplot <<EOF
  set terminal png truecolor size ${hPix},${vPix} font "arial,18"
  set output "${tfile}"
  set nokey
  
  set xlabel "distance (mm)"
  set ylabel "time (s)"

  plot \
    "${dfile}" using 1:2 notitle with linespoints lt 1 lw 1 pt 7 ps 0.75 lc rgb '#0022cc'
EOF

if [[ -s ${tfile} ]]; then
  pfile="${tmp}-r.png"
  convert ${tfile} -resize '50%' ${pfile}
  if [[ ${show} -ne 0 ]]; then display ${pfile}; fi
  rm ${tfile}
  cat ${pfile}
else
  echo "** plot failed" 1>&2 ; exit 1
fi
