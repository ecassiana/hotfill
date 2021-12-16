#! /bin/bash
datadir=tests/results
  for m in 15; do 
    ./plot_hotfill_performance.sh \
      "max band width ${m}" \
      ${datadir}/plot_data_complete_maxband${m}.csv \
      ${datadir}/hotfill_perf_mb${m}.png
done
