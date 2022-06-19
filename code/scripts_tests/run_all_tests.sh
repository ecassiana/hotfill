#! /bin/bash

# Run tests on all modules, in a logical sequence.
# Usage: {$0} [ -noshow ] { "all" | "bug" | "inc" | "gud" }

if [[ "/$1" == "/-noshow" ]]; then
  showop=( "-noshow" ); shift
else
  showop=()
fi

which="$1" # Either "bug", "gud", "inc", or "all" 

if [[  "/${which}" == "/all" ]]; then
  filt="[?+XI]"
elif [[  "/${which}" == "/bug" ]]; then
  filt="[?X]"
elif [[  "/${which}" == "/inc" ]]; then
  filt="[?I]"
elif [[  "/${which}" == "/gud" ]]; then
  filt="[+]"
else
  echo "** invalid option" 1>&2 ; exit 1
fi
echo "filt = '${filt}'"

modules=( `cat 00-MODULES.txt | egrep -e "^${filt}" | cut -d' ' -f2` )

for module in ${modules[@]}; do 
  if [[ ! ( -s ${module}.py ) ]] ; then
    echo "** module ${module}.py does not exist" 1>&2;
  elif [[ ! ( -s tests/${module}_TST.py ) ]]; then
    echo "** test program tests/${module}_TST.py does not exist" 1>&2;
  else
    run_python_test.sh ${showop[@]} ${module}
  fi
done

