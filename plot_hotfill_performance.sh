#! /bin/bash
# Last edited on 2021-11-28 00:50:58 by stolfi

title="$1"; shift    # Plot title.
rawfile="$1"; shift  # Name of input data file.
pngfile="$1"; shift  # Name of output plot PNG file.

tmp="/tmp/$$"

datfile="${tmp}.txt"

pngtemp="${tmp}.png"

cat ${rawfile} \
  | tr -d '\015' \
  | sort -k1,1 -k2,2g \
  | gawk \
      ' BEGIN { part_prev = "?"; np = 0 }
        /[0-9]/ { 
          part = $1; Delta = $2; Trast = $3; Tfab_RP3 = $4; Tfab_HF = $5
          if ((part_prev != "?") && (part != part_prev)) {
            printf "\n"; np++
          }
          if (Tfab_HF == "-") {
            r = -99;
          } else {
            r = (Tfab_HF - Trast)/(Tfab_RP3 - Trast)
          }
          printf "%-20s %5.1f %5.3f %3d\n", part, Delta, r, np;
          part_prev = part
        }
      ' \
  > ${datfile}

export GDFONTPATH=".:${HOME}/ttf" 
gnuplot <<EOF
set term png enhanced font arialb 28 size 1200,1400

set output "${pngtemp}" 
set title "${title}"
unset title

set cbrange [0:24]
set palette defined ( 0 "#ff0000", 0.200 "#dd5500", 0.400 "#008800", 0.600 "#0077ff", 0.800 "#2200ff", 1.000 "#cc0088" )
unset colorbox
# set key top right  
unset key

# set ytics 100
# set mytics 5

set logscale x
set xrange [1.8:170]; set xlabel "cooling time limit (s)"
set xtics ( 2 0, 4 0, 8 0, 16 0, 32 0, 64 0, 128 0, 256 0 )
set grid xtics,mxtics  lt 1 lw 3 lc rgb '#ffddaa', lt 1 lw 2 lc rgb '#ffddaa'

set logscale y
set yrange [69:810]; set ylabel "relative connection time (%)"
set ytics ( 60 0, 70 0, 80 0, 90 0, 100 0, 120 0, 140 0, 160 0, 180 0, 200 0, 240 0, 260 0, 300 0, 400 0, 500 0, 600 0, 700 0, 800 0, 900 0, 1000 0 )
set grid ytics,mytics  lt 1 lw 3 lc rgb '#ffddaa', lt 1 lw 2 lc rgb '#ffddaa'

vmap(k) = (column(k) < 0 ? 0/0 : 100*column(k))

plot "${datfile}" using 2:(vmap(3)):4 title "Tcon_{HF}/Tcon_{RP3}" with linespoints pt 7 ps 2 lw 2 palette
quit
EOF

if [[ -s ${pngtemp} ]] ; then 
  convert ${pngtemp} -resize '50%' ${pngfile}
  eom ${pngfile}
fi

# rm -f ${datfile} ${pngtemp}
