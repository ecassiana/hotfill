#! /usr/bin/python3

import txt_read
import job_parms
import move_parms
import move
import path
import block
import raster
import contact
import hacks
import input_data

import sys
from math import sqrt, sin, cos, log, exp, floor, ceil, inf, nan, pi

parms = job_parms.typical_paper()
mp_jump = move_parms.make_for_jumps(parms)
   
def do_show_data(title, OCRS, OPHS, OLKS, CTS, Z):
  # Prints the given contour paths {OCRS}, fill paths {OPHS}, link paths {OLKS}, and contacts {CTS},
  # as well as the nominal Z-coordinate {Z}.
  
  sys.stderr.write("  ### %s #########################################################################\n" % title)
  sys.stderr.write("Z = %.3f:\n" % Z)
  sys.stderr.write("RASTERS (%d):\n" % len(OPHS))
  path.show_list(sys.stderr, "  ", OPHS, None, True, False)
  sys.stderr.write("CONTOURS (%d):\n" % len(OCRS))
  path.show_list(sys.stderr, "  ", OCRS, None, True, True)
  sys.stderr.write("LINKS (%d):\n" % len(OLKS))
  path.show_list(sys.stderr, "  ", OLKS, None, False, False)
  sys.stderr.write("CONTACTS (%d):\n" % len(CTS))
  contact.show_list(sys.stderr, "  ", CTS, None)
  return 
  # ----------------------------------------------------------------------

def do_test_read(infolder, partname, mp_cont, mp_fill, mp_link, islice, angle, shift, debug):

  sys.stderr.write("--- testing {read} ---\n")
  sys.stderr.write("infolder = %s\n" % infolder)
  sys.stderr.write("partname =  %s\n" % partname)
  sys.stderr.write("islice =  %s\n" % islice)
  sys.stderr.write("angle =  %s\n" % angle)
  sys.stderr.write("shift =  %s\n" % str(shift))
  
  wd_fill = move_parms.width(mp_fill)
  
  tag = partname + ("%03d" % islice)
  
  fname_in = infolder + partname + ("_%03d" % islice) + ".txt"

  rd = open(fname_in, 'r')
  smptol = 0.02 * wd_fill
  fixdir = False
  OCRS, OPHS, OLKS, CTS, Z = txt_read.read(rd, mp_cont, mp_fill, mp_link, angle, shift, fixdir, smptol)
  rd.close()
   
  if debug: do_show_data("read data", OCRS, OPHS, OLKS, CTS, Z)
 
  # ??? Should validate output -- contours, links, contacts, etc. ???
  
  BCS = [ block.from_paths((oph, path.rev(oph))) for oph in OPHS ]

  outfolder = "tests/out/"
  outfile = "txt_read_TST_" + partname + ("_%03d" % islice)
  fname_plot = outfolder + outfile
  CLRS = hacks.trace_colors(len(BCS), None)
  wd_axes = 0.15*wd_fill
  deco = (len(BCS) < 100)
  links = True
  contacts = True
  input_data.plot_input(fname_plot, OCRS, BCS, None, CLRS, wd_axes,deco, links, contacts)  
  return
  # ----------------------------------------------------------------------

def test_read_chain_link_2():

  # Some arbitrary dynamics parameters:
  ac = parms['acceleration'] 
  sp = parms['max_extrusion_speed'] 

  # Move parameters matching the raster spacings in the file:
  wd_cont = 2.00; mp_cont = move_parms.make(wd_cont, ac, sp, 0.0)
  wd_fill = 3.00; mp_fill = move_parms.make(wd_fill, ac, sp, 0.0)
  wd_link = 2.00; mp_link = move_parms.make(wd_link, ac, sp, 0.0)

  partname = "chain_link_2"
  infolder = "tests/in/" + partname + "/"
  islice = 2
  angle = pi/2
  shift = (0, 0)

  do_test_read(infolder, partname, mp_cont, mp_fill, mp_link, islice, angle, shift, False)
  return
  # ----------------------------------------------------------------------

def test_read_runleg():

  # Some arbitrary dynamics parameters:
  ac = parms['acceleration'] 
  sp = parms['max_extrusion_speed'] 

  # Move parameters matching the raster spacings in the file:
  wd_cont = 0.40; mp_cont = move_parms.make(wd_cont, ac, sp, 0.0)
  wd_fill = 0.40; mp_fill = move_parms.make(wd_fill, ac, sp, 0.0)
  wd_link = 0.40; mp_link = move_parms.make(wd_link, ac, sp, 0.0)

  partname = "runleg"
  infolder = "tests/in/"
  islice = 17
  angle = 0
  shift = (0, 0)

  do_test_read(infolder, partname, mp_cont, mp_fill, mp_link, islice, angle, shift, False)
  return
  # ----------------------------------------------------------------------
  
test_read_chain_link_2()
test_read_runleg()
