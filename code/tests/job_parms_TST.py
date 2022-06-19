#! /usr/bin/python3
# Test program for module {job_parms}.

import job_parms
import sys

for tt in ( 
    ( "typical_js",   job_parms.typical_js(), ), 
    ( "typical",      job_parms.typical(), ), 
    ( "slow",         job_parms.slow(), ),
    ( "very_slow",    job_parms.very_slow(), ),
  ):
  name, parms = tt
  sys.stderr.write("%s:\n" % name)
  job_parms.write(sys.stderr, "  ", parms, None)
  
