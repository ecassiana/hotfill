
# Implementation of {input_data}.
# Last edited on 2021-10-29 17:27:25 by stolfi

import pyx
import raster
import raster_example
import paper_example_B
import path
import contact
import move
import move_parms
import block
import hacks

import rn
import rmxn

import sys
from math import sqrt, sin, cos, atan2, log, exp, floor, ceil, inf, nan, pi

def plot_input(fname, OCRS, BCS, o, CLRS, wd_axes,deco, links, contacts):

  assert len(BCS) > 0
  
  # Find the max number of choices {nch}:
  nch = 0
  for bc in BCS: nch = max(nch, block.nchoices(bc))

  # Compute the bounding box {B}:
  B = None
  if o != None: B = rn.box_include_point(B, o)
  B = rn.box_join(B, path.bbox(OCRS))
  B = rn.box_join(B, block.bbox(BCS, links, contacts))
   
  dp = (0,0)
  
  frame = True
  grid = True
  c, szx, szy = hacks.make_canvas(hacks.round_box(B,0.5), dp, None, frame, grid, 1, nch)
  
  ystep = szy # Y shift between choice plots.

  rwd = 0.80
  axes = deco
  dots = deco
  arrows = deco
  matter = True

  clr_o = pyx.color.rgb.black
  wd_o = 3*wd_axes
  
  for ich in range(nch):
    dp = (0, ich*ystep)
    plot_contours(c, OCRS, o, dp, wd_axes, matter)
    block.plot_standard(c, BCS, ich, dp, CLRS, rwd, wd_axes, axes, dots, arrows, matter, links, contacts)

  if o != None: hacks.plot_line(c, o, o, dp, clr_o, wd_o, None)

  hacks.write_plot(c, "%s_ip%02d" % (fname, ich))

  return
  # ----------------------------------------------------------------------

def plot_contours(c, OCRS, org, dp, wd_axes, matter):
  # Plots the contours paths in the list {OCRS},
  # starting at the point closest to {org}.
  
  clr_cr = pyx.color.rgb( 0.500, 0.700, 1.000 )

  rwd = 0.80
  axes = False
  dots = False
  arrows = False

  for oph in OCRS:
    if org != None:
      iorg, dorg = path.find_nearest_move_endpoint(oph, org)
      oph = path.shift_contour(oph, iorg)
      org  = path.pfin(oph)
    assert isinstance(oph, path.Path)
    path.plot_standard(
      c, [oph,], dp, None, [clr_cr,], rwd=rwd, wd_axes=wd_axes,
      axes=axes, dots=dots, arrows=arrows, matter=matter )
  return
  # ----------------------------------------------------------------------

def make_synthetic(dsname, variant, mp_cont, mp_fill, mp_link, mp_jump):
  # Return lists of contours {OCRS}, rasters {OPHS}, links {OLKS}, and contacts {CTS},
  # as well as the nominal {Z} coordinate of the slice.
  
  OCRS = []
  BCS = []
  o = None
  Z = 25.00
  
  if dsname == "patch_array":
    if variant == "1x1x7no":
      islands = False
      cols = 1*3 + 1
      rows = 1
      nx = 3
      ny = 7
    elif variant == "1x1x27no":
      islands = False
      cols = 1*3 + 1
      rows = 1
      nx = 3
      ny = 27
    elif variant == "2x1x27no":
      islands = False
      cols = 2*3 + 1
      rows = 1
      nx = 3
      ny = 27
    elif variant == "2x1x27is":
      islands = True
      cols = 2*4 + 1
      rows = 1
      nx = 3
      ny = 27
    elif variant == "3x1x27no":
      islands = False
      cols = 3*3 + 1
      rows = 1
      nx = 3
      ny = 27
    else:
      assert False, "unknown {variant}"

    OCRS, OPHS, OLKS, CTS = raster_example.patch_array(cols, rows, nx, ny, islands, mp_cont, mp_fill, mp_link)
    # Build blocks from filling rasters:
    BCS = [ block.from_paths([oph, path.rev(oph),]) for oph in OPHS ]

  elif dsname == "paper_fig_B":
    
    OCRS, OPHS, OLKS, CTS, VGS, EGS = paper_example_B.make_turkey(mp_cont, mp_fill, mp_link, mp_jump)
    
    if variant == "rasters":

      # Build blocks from filling rasters:
      BCS = [ block.from_paths([oph, path.rev(oph),]) for oph in OPHS ]
      
    elif variant == "subrivs":
      
      # Build blocks from sub-rivers:
      assert False, "not implemented"

    else:
      assert False, "unknown {variant}"
      
  else:
    assert False, "unknown {dsname}"
    
  Z = 0
  o = (0,0)
  return OCRS, BCS, o, Z
  # ----------------------------------------------------------------------

def show_check_and_plot(fname, OCRS, OPHS, OLKS, CTS, o, ydir, wd_axes,deco):
  
  if OPHS == None: OPHS = []
  assert type(OPHS) is tuple or type(OPHS) is list
  assert len(OPHS) >= 1
  
  ncr = len(OCRS)
  if ncr > 0:
    sys.stderr.write("  ... contours (%d) ...\n" % ncr)
    path.show_list(sys.stderr, "  ", OCRS, None, False, True)
    OMVS_CR = [ path.elem(ocr,imv) for ocr in OCRS for imv in range(path.nelems(ocr)) ]
    sys.stderr.write("    contour moves:\b")
    move.show_list(sys.stderr, "    ", OMVS_CR, None)
  
  nph = len(OPHS)
  if nph > 0:
    sys.stderr.write("  ... raster fill elements (%d) ...\n" % nph)
    path.show_list(sys.stderr, "  ", OPHS, None, True, False)
  
  nlk = len(OLKS)
  if nlk > 0:
    sys.stderr.write("  ... link paths (%d) ...\n" % nlk)
    path.show_list(sys.stderr, "  ", OLKS, None, False, False)
  
  nct = len(CTS)
  if nct > 0:
    sys.stderr.write("  ... contacts (%d) ...\n" % nct)
    contact.show_list(sys.stderr, "  ", CTS, None)

  sys.stderr.write("  ... validating ...\n")
  for oph in OPHS: path.validate(oph)
  path.check_links(OPHS,OLKS)
  contact.check_side_paths(CTS,OPHS,ydir)

  # ??? Should validate more ???
  
  # Make blocks from the raster traces:
  BCS = []
  for oph in OPHS:
    bc = block.from_paths([oph,])
    BCS.append(bc)
  
  links = True
  contacts = True
  CLRS = hacks.trace_colors(nph, None)
  plot_input(fname, OCRS, BCS, o, CLRS, wd_axes,deco, links, contacts)
  return
  # ----------------------------------------------------------------------
