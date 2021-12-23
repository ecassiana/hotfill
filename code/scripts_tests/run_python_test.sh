#! /bin/bash
# Last edited on 2021-10-26 22:22:19 by stolfi

# Usage: "run_python_test.sh [ -noshow ] {MODULE}..."
# Runs "tests/{MODULE}_TST.py" with the proper {PYTHONPATH},
# and sends the standard output (if not empty)
# to tests/out/{MODULE}_TST.txt

if [[ "/$1" == "/-noshow" ]]; then
  show=0; shift
else
  show=1
fi

modules=( "$@" ); shift  # 
opref="tests/out/${modulo}"
ofile="${opref}.html"

while [[ ( $PWD =~ .*/tests ) || ( $PWD =~ .*/src ) ]]; do cd .. ; done
mkdir -p tests/out
# ( cd tests && if [[ ! ( -r images ) ]]; then ln -s ../images; fi )
# ( cd tests && if [[ ! ( -r Makefile ) ]]; then ln -s ../Makefile; fi )
( cd tests && if [[ ! ( -r tests ) ]]; then ln -s ./ tests; fi )
outexts=( ext in txt eps dat gcode gcd png pgm ppm )

for module in "${modules[@]}" ; do 
  echo "=== testing module ${module} =============================" 1>&2
  echo "running tests/${module}_TST.py" 1>&2
  # Check for functions that are not being tested:
  egrep -e '^def ' ${module}.py \
    | sed \
        -e 's/^def *//g' \
        -e 's/[ ]*[\\(].*$//g' \
    | sort \
    > .funcs
  find_python_uses.sh tests/${module}_TST.py ${module} > .uses
  bool 1-2 .funcs .uses > .missing
  if [[ -s .missing ]]; then
    echo "!! warning: these functions are not called explicitly by the test program:" 1>&2
    cat .missing | sed -e 's:^:  :g' 1>&2 
  fi
  cat tests/${module}_TST.py | egrep -e '^[#]+ *test' > .commented
  if [[ -s .commented ]]; then
    echo "!! warning: these lines in the test program are commented out:" 1>&2
    cat .commented | sed -e 's:^:  :g' 1>&2 
  fi
  
  # Remove output files of previous runs:
  opref="tests/out/${module}_TST"
  rm -f ${opref}*
  
  # Run the test program, save output:
  tfile="${opref}.txt" # Text output.
  export PYTHONPATH=".:tests:.." ; \
    time python3 tests/${module}_TST.py > ${tfile}
    exstatus=$?
    echo "status = ${exstatus}" 1>&2  
  
  # Remove empty output files:
  for ext in ${outexts[@]} ; do
    for ff in ${opref}*.${ext}; do
      if [[ ( -e ${ff} ) && ( ! ( -s ${ff} ) ) ]]; then 
        if [[ "${ff}" != "tests/out/${module}_TST.txt" ]]; then
          echo "!! ${ff} is empty" 1>&2
        fi
        rm -fv ${ff}
      fi
    done
  done

  if [[ ( ${exstatus} -eq 0 ) && ( $show -ne 0 ) ]]; then
    # Show EPS outputs, if any:
    for efile in ${opref}*.eps; do
      if [[ -s ${efile} ]]; then 
        evince ${efile}
      fi
    done

    # Plot graphs, if any:
    touch ${opref}_BOGUS.dat # Damn shell...
    for dfile in ${opref}*.dat; do
      # Name of the plot is the first underscore-delimited word after 'TST":
      plotname="${dfile}"
      plotname="${plotname##tests/out/*_TST_}"
      plotname="${plotname%%.dat}"
      plotname="${plotname%%_*}"
      if [[ "/${plotname}" != "/BOGUS" ]]; then
        gscript="plot_${module}_${plotname}.sh"
        # Name of gnuplot output file:
        pfile="${dfile/.dat/.png}"
        if [[ ! ( -s ${dfile} ) ]]; then 
          echo "** data file ${dfile} is empty or missing" 1>&2
        elif [[ ! ( -s ${gscript} ) ]]; then 
          echo "** plot script ${gscript} not found" 1>&2
        else
          ${gscript} ${dfile} > ${pfile}
        fi
      fi
    done
  fi
  
  echo "=== done testing ${module} ==========================="
done
