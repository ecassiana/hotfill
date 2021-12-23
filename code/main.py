import hotfill
#import scanline
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

# --------------------------------------------------------------- #

def do_test_build(OPHS, tag, Delta, maxband, outfile):

  sys.stderr.write("--- testing {build} tag = %s ---\n" % tag)
   
  # Set the contact time limits:
  CTS = set()
  for oph in OPHS:
    for isd in 0,1:
      for ct in path.get_contacts(oph,isd):
        CTS.add(ct)
        contact.set_tcool_limit(ct, Delta)
  CTS = list(CTS)
  
  wd_jump = 0.00; mp_jump = move_parms.make(wd_jump, acc, csp_jump,  udp)

  # Compute input raster fabtime:
  ydir = (0,1)  # Y direction.
  ytol = 1.0e-5 # Tolerance for Y coordinate variation.
  minlen = 0.1  # Minimum raster length.
  NrastA, TrastA, Nlink, Tlink, Njump, Tjump = raster.analyze_fabtime(OPHS, ydir, ytol, minlen)

  start_time = time.time()
  fph, z, i, j, BPHS, BTCV, TUK = hotfill.solve(OPHS, mp_jump, maxband, False)
  end_time = time.time()

  execution_time = end_time - start_time

  if fph == None:
    fname = "tests/out/" + outfile + ".txt"
    sys.stderr.write("writing %s ...\n" % fname)
    wr = open(fname, "w")
    wr.write("CpuTime:  %7.5f\n" % (execution_time))
    wr.close()
    sys.stderr.write("!! {hotfill.solve} returned {None}\n")
  
  else:
    # Compute input raster fabtime:
    NrastB, TrastB, Nlink, Tlink, Njump, Tjump = raster.analyze_fabtime([fph,], ydir, ytol, minlen)

    # ??? Should validate the fullpath ???
    ncold = 1 # Number of coldest contacts to show.
    CTScold = contact.coldest(fph, CTS, ncold)
    contact.show_times(sys.stderr, [fph,], CTScold)
    contact.write_times(outfile, [fph,], CTScold, execution_time, NrastA, TrastA, Nlink, Tlink, Njump, Tjump)

    main_plots.plot_plots(fph, OPHS, BPHS, CTScold, tag, 'hotfill', outfolder)
  return

# --------------------------------------------------------------- #

arguments = sys.argv

partname = arguments[1]
Delta    = float(arguments[2])
maxband  = int(arguments[3])
outfile  = arguments[4]

angle = 0 if "_0.txt" in partname else 1.5708
angle_tag = 0 if "_0.txt" in partname else 90
part = partname.split("_")[0]
shift = (0, 0)

tag = "%s-%s-D%06.0f-mu%04d" % (part, angle_tag, 1000*Delta, maxband)

# Move parameters matching the raster spacings in the file:
wd_cont = 0.40; mp_cont = move_parms.make(wd_cont, acc, csp_trace, 0.0)
wd_fill = 0.40; mp_fill = move_parms.make(wd_fill, acc, csp_trace, 0.0)
wd_link = 0.40; mp_link = move_parms.make(wd_link, acc, csp_trace, 0.0)

infolder = "tests/in/"
outfolder = "tests/out/" + partname.replace(".txt", "") + "_delta" + str(Delta) + "_maxband" + str(maxband) + "/"
if not os.path.exists(outfolder):
  os.mkdir(outfolder)

fname_in = infolder + partname
rd = open(fname_in, 'r')
smptol = 0.0
fixdir = True
OCRS, OPHS, OLKS, CTS, Z = txt_read.read(rd, mp_cont, mp_fill, mp_link, angle, shift, fixdir, smptol)
rd.close()

do_test_build(OPHS, tag, Delta, maxband, outfile)
