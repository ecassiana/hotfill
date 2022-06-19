# Tools to create printer/job parameter tables.

import job_parms_IMP

# A /parameter table/ {parms} is a Python dict that specifies a bunch of
# parameters that determine the computation of move and jump times, the
# details of generated G-code, the nominal widths of traces, etc.
#
# Some parameters depend on the printer. Some on the filament used.
# Some are chosen by the user for a specific job.

def typical_js():
  # Return a parameter table with typical values for testing (JS version).
  return job_parms_IMP.typical_js()

def typical():
  # Return a parameter table with typical values for testing.
  return job_parms_IMP.typical()

def typical_paper():
  # Return a parameter table with typical values for testing (Paper version).
  return job_parms_IMP.typical_paper()

def slow():
  # A parameter table with smaller speeds, accelerations,
  # and nozzle/filament up/down time.  Handy for testing
  # timing functions.
  return job_parms_IMP.slow()
  
def very_slow():
  # A set of parameter to make debugging easier. Acceleration is {+oo}
  # and the extrusion speeds are 1 mm/s, so that the execution time of a
  # trace is equal to the distance in mm. The jump speed is 2, but there
  # is an additional penalty of 1 second at each end, so jumps are
  # better than moves only when the distance is greater than 4 mm.
  return job_parms_IMP.very_slow()

def write(wr, pref, parms, suff):
  # Writes the parameter table {parms} nicely to {wr}, one parameter per line.
  # Each line is prefixed by the string {pref} and suffixed by the string 
  # {suff} (if they are not {None}).
  job_parms_IMP.write(wr, pref, parms, suff) 
  
