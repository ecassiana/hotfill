#! /usr/bin/python3
# Test program for module {bandpath}
# Last edited on 2021-11-25 12:37:20 by stolfi

import bandpath
import paper_example_B
import raster
import raster_example
import path
import move 
import move_parms
import contact
import hacks
import palette
import job_parms

import rn
import rmxn

import pyx
import sys
from math import sqrt, sin, cos, floor, ceil, inf, nan, pi

parms = job_parms.typical_js()

mp_jump = move_parms.make_for_jumps(parms)

# Some arbitrary dynamics parameters:
ac = parms['acceleration'] 
sp = parms['max_extrusion_speed'] 

# Move parameters matching the raster spacings in the file:
wd_cont = 0.75; mp_cont = move_parms.make(wd_cont, ac, sp, 0.0)
wd_fill = 1.00; mp_fill = move_parms.make(wd_fill, ac, sp, 0.0)
wd_link = 0.50; mp_link = move_parms.make(wd_link, ac, sp, 0.0)

wd_axes = 0.15*min(wd_fill,wd_cont)

parms['solid_raster_width'] = wd_fill
parms['contour_trace_width'] = wd_cont

def plot(bph, CTS_lo, CTS_hi, tag):
  # Plots the path {bph} and the contacts in the lists {CTS_lo}, {CTS_hi}.
  # Plot files are "tests/out/bandpath_regroup_TST_{tag}.{ext}"
  # where {ext} is "png", "jpg", "eps".
  
  CTS = CTS_lo + CTS_hi
  nct = len(CTS)
  
  LKS = list(path.get_links(bph) + path.get_links(path.rev(bph)))
  nlk = len(LKS)
  
  # Compute the bounding box of it all:
  B = None
  B = rn.box_join(B, path.bbox([bph,]))
  B = rn.box_join(B, path.bbox(LKS))
  B = rn.box_join(B, contact.bbox(CTS))
  B = hacks.round_box(B,0.5)
  
  dp = (0,0)
  
  frame = False
  grid = True
  c,szx,szy = hacks.make_canvas(B, dp, None, frame, grid, 1, 1)
  
  wd_axes = 0.05*wd_fill

  # Paint the material traces:

  rwd_ph = 0.80
  axes_ph = False
  dots_ph = False
  arrows_ph = False
  matter_ph = True
  path.plot_standard(c, [bph,], dp, 0, None, rwd_ph, wd_axes, axes_ph, dots_ph, arrows_ph, matter_ph)

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
  path.plot_standard(c, [bph,], dp, None, [clr_ph,], rwd_ph, wd_axes, axes_ph, dots_ph, arrows_ph, matter_ph)
    
  if nct > 0:
    clr_ct = pyx.color.rgb(1.000, 0.100, 0.000)
    wd_ct = 2*wd_axes
    sz_tics_ct = 0
    arrows_ct = False
    for ct in CTS:
      contact.plot_single(c, ct, dp, clr_ct, None, 0, wd_ct, sz_tics_ct, arrows_ct)

  fname = "tests/out/bandpath_TST" "_" + tag
  hacks.write_plot(c, fname)
  return
  # ----------------------------------------------------------------------

def do_test_build(OPHS, SCS, data_tag, Delta, z, i, j, quick):

  sys.stderr.write("--- testing {build} data_tag = %s z = %d i = %d j = %d ---\n" % (data_tag,z,i,j))
   
  tag = "%s_z%d_i%03d_j%03d" % (data_tag,z,i,j)
  
  # Set the contact time limits:
  for oph in OPHS:
    for isd in 0,1:
      for ct in path.get_contacts(oph,isd):
        contact.set_tcool_limit(ct, Delta)
  
  SCS_band = SCS[i:j]
  OPHS_band = []
  for Si in SCS_band:
    for iph in Si:
      OPHS_band.append(OPHS[iph])

  ydir = (0,1)
  ytol = 1.0e-6
  minlen = 0.1
  Trast_in, Tlink_in, Tjump_in = raster.analyze_fabtime(OPHS_band, ydir, ytol, minlen)
  sys.stderr.write("Input: Trast = %.6f  Tlink = %.6f  Tjump = %.6f\n" % (Trast_in,Tlink_in,Tjump_in))
  assert Tlink_in == 0
  assert Tjump_in == 0
  
  bph, CTS_lo, TCVS_lo, CTS_hi, TCVS_hi = bandpath.build(OPHS, SCS_band, z, mp_jump, quick)
  Tfab = path.fabtime(bph)
  sys.stderr.write("Tfab = %.6f\n" % Tfab)

  Trast_ot, Tlink_ot, Tjump_ot = raster.analyze_fabtime([bph,], ydir, ytol, minlen)
  sys.stderr.write("Output: Trast = %.6f  Tlink = %.6f  Tjump = %.6f\n" % (Trast_ot,Tlink_ot,Tjump_ot))
  
  if bph == None:
    sys.stderr.write("!! {bandpath.build} returned {None}\n")
  else:
    # ??? Should validate the bandpath ???
    # OPHS_band = [ OPHS[iph] for IPHS in SCS for iph in IPHS ]

    plot(bph, CTS_lo, CTS_hi, tag)
  return
  # ----------------------------------------------------------------------
  
def test_patch_array(cols, rows, nx, ny, quick):
  # Tests {bandpath.build} with the {raster_example.patch_array} dataset.
  # 

  tag = "patch_array_c%02d_r%02d_nx%02d_ny%02d_q%d" % (cols, rows, nx,ny,int(quick))

  islands = False
  OCRS, OPHS, OLKS, CTS = raster_example.patch_array(cols, rows, nx,ny, islands, mp_cont,mp_fill,mp_link) 
  nph = len(OPHS)

  xdir = (1,0)
  ydir = (0,1)
  ystep_exp = wd_fill
  yphase_exp = 2*wd_fill
  eps = 1.0e-6
  ystep, yphase = raster.get_spacing_and_phase(OPHS, xdir, ydir)
  sys.stderr.write("ystep = %20.16f yphase = %20.16f\n" % (ystep, yphase))
  assert abs(ystep - ystep_exp) <= 0.5*eps
  assert abs(yphase - yphase_exp) <= 0.5*eps
  
  OPHS = raster.sort_by_scanline(OPHS, xdir, ydir, ystep, yphase)
  SCS = raster.separate_by_scanline(OPHS, xdir, ydir, ystep, yphase)
  
  Delta = 60.0
  do_test_build(OPHS, SCS, tag, Delta, 0, 3, 8, quick)
  
  return 
  # ----------------------------------------------------------------------
   
def test_paper_example_B(quick):
  # Tests {bandpath.build} with the {paper_example_B.make_turkey} dataset.

  tag = "paper_example_B_q%d" % int(quick)

  OCRS,OPHS,OLKS,CTS,VGS,EGS = paper_example_B.make_turkey(mp_cont, mp_fill, mp_link, mp_jump)
  nph = len(OPHS)

  xdir = (1,0)
  ydir = (0,1)
  ystep_exp = wd_fill
  yphase_exp = wd_fill
  eps = 1.0e-6
  ystep, yphase = raster.get_spacing_and_phase(OPHS, xdir, ydir)
  sys.stderr.write("ystep = %20.16f yphase = %20.16f\n" % (ystep, yphase))
  assert abs(ystep - ystep_exp) <= 0.5*eps
  assert abs(yphase - yphase_exp) <= 0.5*eps
  
  OPHS = raster.sort_by_scanline(OPHS, xdir, ydir, ystep, yphase)
  SCS = raster.separate_by_scanline(OPHS, xdir, ydir, ystep, yphase)
  
  Delta = 60.0
  do_test_build(OPHS, SCS, tag, Delta, 0, 2, 17, quick)
  
  return 
  # ----------------------------------------------------------------------
  
### test_paper_example_B(True)
test_paper_example_B(False)

### cols = 2*3 + 1
### rows = 1
### nx = 4
### ny = 7
### test_patch_array(cols, rows, nx, ny, True)
### test_patch_array(cols, rows, nx, ny, False)
