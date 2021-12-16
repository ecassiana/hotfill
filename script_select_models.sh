#! /bin/bash

outnew=outfiltered

outpath=./tests/

if [[ ! -e ./tests/${outnew}/ ]]; then
    mkdir ./tests/${outnew}/
fi

rm -rf ${outpath}/${outnew}

mkdir ${outpath}/${outnew}

cp -r ${outpath}/out/adfoot_*_0* ${outpath}/${outnew}
cp -r ${outpath}/out/bkwren_*_0* ${outpath}/${outnew}
cp -r ${outpath}/out/cthang_*_90* ${outpath}/${outnew}
cp -r ${outpath}/out/dstop1_*_0* ${outpath}/${outnew}
cp -r ${outpath}/out/dstop2_*_90* ${outpath}/${outnew}
cp -r ${outpath}/out/flange_*_0* ${outpath}/${outnew}
cp -r ${outpath}/out/grille_*_0* ${outpath}/${outnew}
cp -r ${outpath}/out/hlatch_*_90* ${outpath}/${outnew}
cp -r ${outpath}/out/nkhang_*_90* ${outpath}/${outnew}
cp -r ${outpath}/out/runleg* ${outpath}/${outnew}
cp -r ${outpath}/out/tkwren_*_90* ${outpath}/${outnew}
cp -r ${outpath}/out/unwren_*_90* ${outpath}/${outnew}
