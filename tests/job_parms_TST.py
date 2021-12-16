#! /usr/bin/python3
# Test program for module {job_parms}.
# Last edited on 2021-05-06 15:29:37 by jstolfi

import job_parms
import sys

for tt in ( 
    ( "typical_js",   job_parms.typical_js(), ), 
    ( "typical_elis", job_parms.typical_elis(), ), 
    ( "slow",         job_parms.slow(), ),
    ( "very_slow",    job_parms.very_slow(), ),
  ):
  name, parms = tt
  sys.stderr.write("%s:\n" % name)
  job_parms.write(sys.stderr, "  ", parms, None)
  
