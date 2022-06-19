#! /usr/bin/python3
# Test program for module {hotfill}

import hotfill
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

import rn
import rmxn

import pyx
import sys
import time
from math import sqrt, sin, cos, floor, ceil, inf, nan, pi

# Some arbitrary dynamics parameters:
acc = 3000      # Acceleration/deceleration for moves and jumps (mm/s^2).
csp_trace = 40  # Cruise speed for traces (mm/s).
csp_jump = 130  # Cruise speed for jumps (mm/s).
udp = 0.05      # Trace/jump transition penalty (s).

def plot(fph, BPHS, CTS, tag):
  # Plots the fullpath {fph} and the individual bandpaths and connectors {BPHS}.
  # Also plots the contacts in {CTS}.
  # 
  # Plot files are "tests/out/hotfill_TST_{tag}.{ext}"
  # where {ext} is "png", "jpg", "eps".
  
  LKS = list(path.get_links(fph) + path.get_links(path.rev(fph)))
  nlk = len(LKS)
  
  # Compute the bounding box of it all:
  B = None
  B = rn.box_join(B, path.bbox([fph,]))
  B = rn.box_join(B, path.bbox(LKS))
  B = hacks.round_box(B,0.5)
  
  dp = (0,0)
  
  frame = False
  grid = True
  c,szx,szy = hacks.make_canvas(B, dp, None, frame, grid, 1, 2)
  
  mp_fill = move.parameters(path.elem(BPHS[0],0))
  wd_fill = move_parms.width(mp_fill)
  wd_axes = 0.05*wd_fill

  # PLOT THE FULLPATH
  
  # Paint the material traces:

  rwd_ph = 0.80
  axes_ph = False
  dots_ph = False
  arrows_ph = False
  matter_ph = True
  path.plot_standard(c, [fph,], dp, 0, None, rwd_ph, wd_axes, axes_ph, dots_ph, arrows_ph, matter_ph)

  # Paint the nominal traces of the links:

  if nlk > 0:
    CLRS_lk = hacks.link_colors(nlk, None)
    rwd_lk = 0.60
    axes_lk = True
    dots_lk = True
    arrows_lk = False
    matter_lk = False
    path.plot_standard(c, LKS, dp, None, CLRS_lk, rwd_lk, wd_axes, axes_lk, dots_lk, arrows_lk, matter_lk)

  # Paint the nominal traces of the path:

  clr_ph = pyx.color.rgb( 0.000, 0.800, 1.000 )
  rwd_ph = 0.80
  axes_ph = False
  dots_ph = False
  arrows_ph = True
  matter_ph = False
  path.plot_standard(c, [fph,], dp, None, [clr_ph,], rwd_ph, wd_axes, axes_ph, dots_ph, arrows_ph, matter_ph)
    
  if len(CTS) > 0:
    clr_ct = pyx.color.rgb( 1.000, 0.000, 0.700 )
    wd_ct = 2.0*wd_axes
    dashpat_ct = None
    ext_ct = 0
    sz_tics_ct = 0
    arrows_ct = False
    for ct in CTS:
      plotit = True
      contact.plot_single(c, ct, dp, clr_ct, dashpat_ct, ext=ext_ct, wd=wd_ct, sz_tic=sz_tics_ct, arrow=arrows_ct)

  # PLOT THE BANDPATHS
  
  dp = rn.add(dp, (0, szy))

  # Separate the bandpaths from the connectors:
  PPHS = []
  for iph in range(len(BPHS)):
    if iph % 2 == 0:
      PPHS.append(BPHS[iph])
      
  nbp = len(PPHS) # Number of bandpaths.
  CLRS_bp = hacks.trace_colors(nbp, None)
  rwd_bp = 0.80
  axes_bp = False
  dots_bp = False
  arrows_bp = True
  matter_bp = False
  path.plot_standard(c, PPHS, dp, None, CLRS_bp, rwd_bp, wd_axes, axes_bp, dots_bp, arrows_bp, matter_bp)
  
  fname = "tests/out/hotfill_TST" "_" + tag
  hacks.write_plot(c, fname)
  return
  # ----------------------------------------------------------------------

def do_test_build(OPHS, tag, Delta, maxband, quick):

  sys.stderr.write("--- testing {build} tag = %s ---\n" % tag)
  
  debug = False
   
  # Set the contact time limits:
  CTS = set()
  for oph in OPHS:
    for isd in 0,1:
      for ct in path.get_contacts(oph,isd):
        CTS.add(ct)
        contact.set_tcool_limit(ct, Delta)
  CTS = list(CTS)
  
  wd_jump = 0.00; mp_jump = move_parms.make(wd_jump, acc, csp_jump,  udp)
  move_parms.show(sys.stderr, "mp_jump = ", mp_jump, "\n")

  hotfill.describe_and_check_input(sys.stderr, OPHS, mp_jump, maxband, quick)
  tbeg = time.clock_gettime(time.CLOCK_THREAD_CPUTIME_ID)
  fph, z, i, j, BPHS, BTCV, TUK = hotfill.solve(OPHS, mp_jump, maxband, quick)
  tend = time.clock_gettime(time.CLOCK_THREAD_CPUTIME_ID)
  Tcpu = tend - tbeg
  
  if debug:
    # Print all moves:
    for imv in range(path.nelems(fph)):
      omv = path.elem(fph, imv)
      move.show(sys.stderr, "  fph[%05d] " % imv, omv, "\n", 0)

  if fph == None:
    sys.stderr.write("!! {hotfill.solve} returned {None}\n")
    sys.stderr.write("Tcpu = %.2f s\n" % Tcpu)
  else:
    sys.stderr.write("=== results ===\n")
    sys.stderr.write("tag = %s\n" % tag)
    sys.stderr.write("maxband = %d quick = %s\n" % (maxband, quick))
    
    # Compute input raster fabtime:
    ydir = (0,1)  # Y direction.
    ytol = 1.0e-5 # Tolerance for Y coordinate variation.
    minlen = 0.1  # Minimum raster length.
    Trast_in, Tlink_in, Tjump_in = raster.analyze_fabtime(OPHS, ydir, ytol, minlen)
    assert Tlink_in == 0
    assert Tjump_in == 0
    sys.stderr.write("Trast (input) = %.6f\n" % Trast_in)
    
    # Check output raster fabtime (paranoia):
    Trast_ot, Tlink_ot, Tjump_ot = raster.analyze_fabtime([fph,], ydir, ytol, minlen)
    sys.stderr.write("Trast (output) = %.6f\n" % Trast_ot)
    sys.stderr.write("Tlink (output) = %.6f\n" % Tlink_ot)
    sys.stderr.write("Tjump (output) = %.6f\n" % Tjump_ot)

    # ??? Should validate the fullpath ???
    ncold = 3 # Number of coldest contacts to show.
    CTScold = contact.coldest(fph, CTS, ncold)
    contact.show_times(sys.stderr, [fph,], CTScold)
    sys.stderr.write("Tcpu = %.2f s\n" % Tcpu)
    assert Trast_ot >= Trast_in
    if abs(Trast_ot - Trast_in) > 1.0e-8:
      sys.stderr.write("!! warning - Trast out greater than Trast in\n")
    sys.stderr.write("===============\n")

    plot(fph, BPHS, CTScold, tag)
  return
  # ----------------------------------------------------------------------
  
def test_patch_array(cols, rows, nx, ny, Delta, maxband, quick):
  # Tests {hotfill.build} with the {raster_example.patch_array} dataset.

  # Move parameters matching the raster spacings in the file:
  wd_cont = 0.75; mp_cont = move_parms.make(wd_cont, acc, csp_trace, 0.0)
  wd_fill = 1.00; mp_fill = move_parms.make(wd_fill, acc, csp_trace, 0.0)
  wd_link = 0.50; mp_link = move_parms.make(wd_link, acc, csp_trace, 0.0)

  tag_dim = "c%02d_r%02d_nx%02d_ny%02d" % (cols, rows, nx,ny)
  tag_mbq = "mb%03d_q%d" % (maxband, int(quick))
  tag = "patch_array_%s_D%06.0f_%s" % (tag_dim, 1000*Delta, tag_mbq)

  islands = False
  OCRS, OPHS, OLKS, CTS = raster_example.patch_array(cols, rows, nx,ny, islands, mp_cont,mp_fill,mp_link) 

  do_test_build(OPHS, tag, Delta, maxband, quick)
  
  return 
  # ----------------------------------------------------------------------
  
def test_runleg(angle, Delta, maxband, quick):
  # Tests {hotfill.build} with the "runleg" data set (angle 0)

  # Move parameters matching the raster spacings in the file:
  wd_cont = 0.40; mp_cont = move_parms.make(wd_cont, acc, csp_trace, 0.0)
  wd_fill = 0.40; mp_fill = move_parms.make(wd_fill, acc, csp_trace, 0.0)
  wd_link = 0.40; mp_link = move_parms.make(wd_link, acc, csp_trace, 0.0)

  move_parms.show(sys.stderr, "mp_fill = ", mp_fill, "\n")
  move_parms.show(sys.stderr, "mp_link = ", mp_link, "\n")

  partname = "runleg"
  islice = 17

  tag_mbq = "mb%03d_q%d" % (maxband, int(quick))
  tag = "%s_slc%03d_ang%03d_D%06.0f_%s" % (partname, angle, islice, 1000*Delta, tag_mbq)

  infolder = "tests/in/"

  fname_in = infolder + partname + ("_%03d" % islice) + ".txt"
  rd = open(fname_in, 'r')

  shift = (0, 0)
  smptol = 0.0
  fixdir = True
  OCRS, OPHS, OLKS, CTS, Z = txt_read.read(rd, mp_cont, mp_fill, mp_link, angle, shift, fixdir, smptol)

  rd.close()

  do_test_build(OPHS, tag, Delta, maxband, quick)
  
  return 
  # ----------------------------------------------------------------------
   
def test_paper_example_B(Delta, maxband, quick):
  # Tests {hotfill.build} with the {paper_example_B.make_turkey} dataset.

  tag_mbq = "mb%03d_q%d" % (maxband, int(quick))
  tag = "paper_example_B_%s" % tag_mbq

  # Move parameters matching the raster spacings in the file:
  wd_cont = 0.75; mp_cont = move_parms.make(wd_cont, acc, csp_trace, 0.0)
  wd_fill = 1.00; mp_fill = move_parms.make(wd_fill, acc, csp_trace, 0.0)
  wd_link = 0.50; mp_link = move_parms.make(wd_link, acc, csp_trace, 0.0)
  wd_jump = 0.00; mp_jump = move_parms.make(wd_jump, acc, csp_jump,  udp)

  move_parms.show(sys.stderr, "mp_fill = ", mp_fill, "\n")
  move_parms.show(sys.stderr, "mp_link = ", mp_link, "\n")
  move_parms.show(sys.stderr, "mp_jump = ", mp_jump, "\n")

  wd_axes = 0.15*min(wd_fill,wd_cont)
  OCRS,OPHS,OLKS,CTS,VGS,EGS = paper_example_B.make_turkey(mp_cont, mp_fill, mp_link, mp_jump)

  do_test_build(OPHS, tag, Delta, maxband, quick)
  
  return 
  # ----------------------------------------------------------------------
  
Delta = 8
maxband = 20
angle = 0
test_runleg(angle, Delta, maxband, False)
### test_runleg(angle, Delta, maxband, True)

### Delta = 3.5
### maxband = 30
### test_paper_example_B(Delta, maxband, False)

### cols = 2*3 + 1
### rows = 1
### nx = 4
### ny = 20
### Delta = 20.0
### maxband = 20
### test_patch_array(cols, rows, nx, ny, Delta, maxband, quick)
