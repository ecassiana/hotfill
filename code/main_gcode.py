import hotfill
import gcode_read
import paper_example_B
import raster
import raster_example
import path
import move 
import move_parms
import contact
import hacks
import palette
import txt_read
import job_parms
import os
import rn
import rmxn
import time
import pyx
import sys
from math import sqrt, sin, cos, floor, ceil, inf, nan, pi
import main_plots

##
# Some arbitrary dynamics parameters:
acc = 3000      # Acceleration/deceleration for moves and jumps (mm/s^2).
csp_trace = 40  # Cruise speed for traces (mm/s).
csp_jump = 130  # Cruise speed for jumps (mm/s).
udp = 0.05      # Trace/jump transition penalty (s).
outfolder = ""

slice_dict = {
  'adfoot': 5,
  'runleg': 17
}

# --------------------------------------------------------------- #

def do_test_build(OPHS, tag, alg, outfile):

  sys.stderr.write("--- testing {build} tag = %s ---\n" % tag)

  # Compute input raster fabtime:
  ydir = (0,1)  # Y direction.
  ytol = 1.0e-5 # Tolerance for Y coordinate variation.
  minlen = 0.1  # Minimum raster length.
  
  NrastA, TrastA, Nlink, Tlink, Njump, Tjump = raster.analyze_fabtime(OPHS, ydir, ytol, minlen)

  # Set the contact time limits:
  CTS = set()
  for oph in OPHS:
    for isd in 0,1:
      for ct in path.get_contacts(oph,isd):
        CTS.add(ct)
        contact.set_tcool_limit(ct, 1)
  CTS = list(CTS)
  
  wd_jump = 0.00; mp_jump = move_parms.make(wd_jump, acc, csp_jump, udp)

  start_time = time.time()
  fph = gcode_read.solve(OPHS, mp_jump)
  end_time = time.time()

  execution_time = end_time - start_time

  # Compute input raster fabtime:
  NrastB, TrastB, Nlink, Tlink, Njump, Tjump = raster.analyze_fabtime([fph,], ydir, ytol, minlen)

  if fph == None:
    fname = "tests/out/" + outfile + ".txt"
    sys.stderr.write("writing %s ...\n" % fname)
    wr = open(fname, "w")
    wr.write("RunTime:  %7.5f\n" % (execution_time))
    wr.close()
    sys.stderr.write("!! {hotfill.solve} returned {None}\n")
  else:
    # ??? Should validate the fullpath ???
    ncold = 1 # Number of coldest contacts to show.
    CTScold = contact.coldest(fph, CTS, ncold)
    contact.show_times(sys.stderr, [fph,], CTScold)
    contact.write_times(outfile, [fph,], CTScold, execution_time, NrastA, TrastA, Nlink, Tlink, Njump, Tjump)

    main_plots.plot_plots(fph, OPHS, OPHS, CTScold, tag, alg, outfolder)

  return

def get_slice(partname):

  for key in slice_dict.keys():
    if key in partname:
      return slice_dict[key]

  return 1

def get_angle(partname):

  angle = None

  for key in slice_dict.keys():
    if key in partname:
      angle = 1.5708 if "_0.gcode" in partname else 0

  if angle == None:
    angle = 0 if "_0.gcode" in partname else 1.5708

  return angle

# --------------------------------------------------------------- #

arguments = sys.argv

if '.gcode' not in arguments[1]:
  exit()

partname = arguments[1]
alg = 'RP3' if 'rp3' in arguments[1] else 'SLIC3R'
outfile  = arguments[2]

angle = 0 if "_0.gcode" in partname else 1.5708

islice = get_slice(partname)
shift = (0, 0)

tag = "%s_slc%04d_%s" % (partname.replace(".gcode", ""), islice, alg)

# Move parameters matching the raster spacings in the file:
wd_cont = 0.40; mp_cont = move_parms.make(wd_cont, acc, csp_trace, 0.0)
wd_fill = 0.40; mp_fill = move_parms.make(wd_fill, acc, csp_trace, 0.0)
wd_link = 0.40; mp_link = move_parms.make(wd_link, acc, csp_trace, 0.0)

infolder = "tests/in_gcodes/"
outfolder = "tests/out/" + partname.replace(".gcode", "") + "_" + str(alg) + "/"
if not os.path.exists(outfolder):
  os.mkdir(outfolder)

wd_jump = 0.00; mp_jump = move_parms.make(wd_jump, acc, csp_jump, udp)
fname_in = infolder + partname
fname_out = partname.replace(".gcode", "") + "_" + str(alg) + "/" + outfile
rd = open(fname_in, 'r')
OPHS, Z = gcode_read.read(rd, islice, angle, shift, mp_cont, mp_fill, mp_link, mp_jump)
rd.close()

do_test_build(OPHS, tag, alg, outfile)
