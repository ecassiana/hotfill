#! /bin/bash
# Last edited on 2021-11-09 01:27:22 by stolfi

check_dir_tail.sh "2021-10-03-js"

todir="tests/out"
topref="${todir}/paper_figures_B_TST_p"

today=`yyyy-mm-dd-hhmmss`
svdir=figures/${today}

mkdir -p ${svdir}

for ff in ${topref}*.{png,txt,tex} ; do 
  oname=${ff/${topref}_/}
  nname=${oname//_/-}
  # echo "fname = ${fname}" 2>&2
  mv -vi "${topref}_${oname}" "${svdir}/${nname}"
done

