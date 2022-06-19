#! /bin/usr/python3
# Test program for module {contact_example}

import contact_example

import move
import move_example
import move_parms
import path
import path_example
import block
import block_example
import contact
import palette
import hacks
import job_parms
import rn
import pyx
import sys
from math import sqrt, sin, cos, floor, ceil, pi, nan, inf

parms = job_parms.typical_js()
parms['solid_raster_width'] = 1.00
parms['contour_trace_width'] = 0.50

mp_jump = move_parms.make_for_jumps(parms)
mp_cont = move_parms.make_for_contours(parms)
mp_fill = move_parms.make_for_fillings(parms)

wd_fill = move_parms.width(mp_fill)
wd_cont = move_parms.width(mp_cont)

def plot_paths(tag, OPHS, CTS, deco, ct_arrows):
  # Plots to files called "tests/out/path_example_TST_{tag}" the oriented paths in the list
  # {OPHS} and the contacts in the list {CTS}. If {deco} is true also prints the matter footprint and draws
  # axes, dots, and arrowheads on the traces.  If {ct_arrows} is true, prints arrowheads 
  # on the contacts instead of tics.
  
  if OPHS == None: OPHS = []
  assert type(OPHS) is tuple or type(OPHS) is list
  assert len(OPHS) >= 1
  
  if CTS == None: CTS = []
  assert type(CTS) is tuple or type(CTS) is list

  B = path.bbox(OPHS)
  
  dp = None

  c, szx,szy = hacks.make_canvas(hacks.round_box(B,0.5), dp, None, True, True, 1, 1)

  nph = len(OPHS)
  CLRS = hacks.trace_colors(nph, None)

  # Plot the paths:
  rwd = 0.80
  wd_axes = 0.15*min(wd_fill,wd_cont) # Width of jumps and axis lines.
  axes = deco
  dots = deco
  arrows = True
  matter = deco
  path.plot_standard(c, OPHS, dp, None, CLRS, rwd, wd_axes, axes, dots, arrows, matter)

  # Plot the contacts:
  clr_ct =   pyx.color.rgb(1.000, 0.000, 0.000) # Color of contacts.
  wd_ct = 1.5*wd_axes   # Width of contact lines.
  arrow_ct = ct_arrows
  sz_tic_ct = 0 if arrow_ct else wd_ct
  for ct in CTS: 
    contact.plot_single(c, ct, dp, clr_ct, None, 0, wd=wd_ct, sz_tic=sz_tic_ct, arrow=arrow_ct)

  hacks.write_plot(c, "tests/out/contact_example_TST_" + tag)
  return
  # ----------------------------------------------------------------------

def plot_blocks(tag, BCS, CTS, deco, ct_arrows):
  # Plots to files called "tests/out/path_example_TST_{tag}" the blocks
  # in the list {BCS} and the contacts in the list {CTS}. If {deco} is
  # true also prints the matter footprint and draws axes, dots, and
  # arrowheads on the traces. If {ct_arrows} is true, prints arrowheads
  # on the contacts instead of tics.
  
  if BCS == None: BCS = []
  assert type(BCS) is tuple or type(BCS) is list
  assert len(BCS) >= 1
  
  if CTS == None: CTS = []
  assert type(CTS) is tuple or type(CTS) is list

  B = block.bbox(BCS, True, True)
  
  nbc = len(BCS)
  CLRS = hacks.trace_colors(nbc, None)
  
  # Get the max number of choices {nch} among the blocks:
  nch = 0;
  for bc in BCS: nch = max(nch, block.nchoices(bc))

  dp = (0,0)

  c, szx,szy = hacks.make_canvas(hacks.round_box(B,0.5), dp, None, True, True, 1, nch)
  ystep = szy

  # Plot the blocks:
  assert ystep != None
  rwd = 0.80
  wd_axes = 0.15*min(wd_fill,wd_cont) # Width of jumps and axis lines.
  axes = deco
  dots = deco
  arrows = True
  matter = deco
  ystep = szy
  links = True
  contacts = True
  for ich in range(nch):
    dp = (0, ich*ystep)
    block.plot_standard(c, BCS, ich, dp, CLRS, rwd, wd_axes, axes, dots, arrows, matter, links, contacts)

  # Plot the contacts (only on choice 1 for now):
  clr_ct =   pyx.color.rgb(1.000, 0.000, 0.000) # Color of contacts.
  wd_ct = 1.5*wd_axes   # Width of contact lines.
  arrow_ct = ct_arrows
  sz_tic_ct = 0 if arrow_ct else wd_ct
  for ct in CTS: 
    contact.plot_single(c, ct, dp, clr_ct, None, 0, wd=wd_ct, sz_tic=sz_tic_ct, arrow=arrow_ct)

  hacks.write_plot(c, "tests/out/contact_example_TST_" + tag)
  return
  # ----------------------------------------------------------------------

def test_raster_raster_contact():
  
  tag = "raster_raster_contact"
  
  mp_trace = mp_fill
  wd = move_parms.width(mp_trace)
  
  # Build two snake paths {OPHS[0..1]}, one above the other, with a contact in between 
  n = [4, 3]       # Number of rasters in each path.
  plo = [None,None] # Lower left corner of bbox each path.
  plo[0] = ( 2, 1 )
  plo[1] = ( 1, plo[0][1] + n[0]*wd )
  szx = [3, 6] # X extent of traces of each path.
  OPHS = [None, None]
  for isd in range(2):
    # Build a snake block but use only the first choice of each block:
    axis = 0
    alt = True
    PHSi,TRSi,LJS0i,LJS1i = path_example.raster_rectangle(plo[isd], 0, n[isd], alt, szx[isd], wd, mp_trace,None)
    OPHS[isd] = PHSi[0]
  
  # Create the contact {ct} between them:
  ct = contact_example.raster_raster_contact(OPHS[0], OPHS[1])
  
  # It must be between the last trace of {OPHS[0]} and the first one of {OPHS[1]}:
  nmv0 = path.nelems(OPHS[0])
  mv0,dr0 = move.unpack(path.elem(OPHS[0],nmv0-1))
  nmv1 = path.nelems(OPHS[1])
  mv1,dr1 = move.unpack(path.elem(OPHS[1],0))
  assert contact.side_move(ct, 0) == mv0
  assert contact.side_move(ct, 1) == mv1
  
  # Check the path-contact links:
  ph = [None,None]
  dr = [None,None]
  imv = [None,None]
  ph[0], dr[0] = path.unpack(OPHS[0]); imv[0] = path.find_move(OPHS[0],mv0)
  ph[1], dr[1] = path.unpack(OPHS[1]); imv[1] = path.find_move(OPHS[1],mv1)
  assert imv[0] == path.nelems(OPHS[0]) - 1
  assert imv[1] == 0
  for isd in range(2):
    tr_exp = (ph[isd],dr[isd],imv[isd])  # Triplet expected from {get_side_paths}
    assert contact.get_side_paths(ct, isd) == set((tr_exp,))
    assert path.get_contacts(OPHS[isd], isd) == set((ct,))
    assert path.get_contacts(path.rev(OPHS[isd]), isd) == set((ct,))
  
  deco = True
  ct_arrows = True
  plot_paths(tag, OPHS, [ct,], deco, ct_arrows)
  return
  # ----------------------------------------------------------------------

def test_two_roads_and_islands(nmv, nmg, nis):
  sys.stderr.write("--- testing {two_roads_and_islands} nmv = %d nmg = %d nis = %d ---\n" % (nmv,nmg,nis))
  tag = "two_roads_and_islands_nrd%02d_nbc%02d_nmv%03d" % (nmv,nmg,nis)
  BCS, CTS = contact_example.two_roads_and_islands(nmv, nmg, nis, mp_cont, mp_fill, mp_jump)
  block.show_list(sys.stderr, "    ", BCS, None, True)
  contact.show_list(sys.stderr, "    ", CTS, None)
  # ??? seam.print_contact_table(sys.stderr, BCS, CTS, None)
  o = (1,1)
  deco = True
  ct_arrows = False
  plot_blocks(tag, BCS, CTS, deco, ct_arrows)
  return
  # ----------------------------------------------------------------------
 
def test_misc_A():
  sys.stderr.write("--- testing {plot,misc_A} ---\n")
  tag = "misc_A"
  ct = contact_example.misc_A(mp_fill)
  tr0 = contact.side_move(ct,0)
  tr1 = contact.side_move(ct,1)
  CTS = [ct,]
  OPHS = [ path.from_moves([tr0,]), path.from_moves([tr1,]), ]
  deco = True
  ct_arrows = False
  plot_paths(tag, OPHS, [ct,], deco, ct_arrows)
  return
  # ----------------------------------------------------------------------

def test_misc_B(ct_arrows):
  sys.stderr.write("--- testing {plot,misc_B} with arrow = %s ---\n" % str(ct_arrows))
  tag = "misc_B"
  CTS, OPHS, TRS = contact_example.misc_B(mp_fill, mp_jump)
  deco = True
  plot_paths(tag, OPHS, CTS, deco, ct_arrows)
  return
  # ----------------------------------------------------------------------

def test_misc_C():
  sys.stderr.write("--- testing {plot,misc_C} ---\n")
  tag = "misc_C"
  CTS, OPHS = contact_example.misc_C(mp_fill)
  deco = True
  ct_arrows = True
  plot_paths(tag, OPHS, CTS, deco, ct_arrows)
  return
  # ----------------------------------------------------------------------

def test_misc_F(alt):
  sys.stderr.write("--- testing {plot,misc_F} alt = %d ---\n" % int(alt))
  tag = "misc_F_alt%d" % int(alt)
  ct, ph = contact_example.misc_F(alt, mp_fill, mp_jump)
  tr0 = contact.side_move(ct,0)
  tr1 = contact.side_move(ct,1)
  CTS = [ ct, ]
  OPHS = [ ph, ]
  deco = True
  ct_arrows = False
  plot_paths(tag, OPHS, [ct,], deco, ct_arrows)
  return
  # ----------------------------------------------------------------------

def test_misc_K():
  sys.stderr.write("--- testing {plot,misc_K} ---\n")
  tag = "misc_K"
  CTS, OPHS = contact_example.misc_K(mp_fill, mp_jump)
  deco = True
  ct_arrows = True
  plot_paths(tag, OPHS, CTS, deco, ct_arrows)
  return
  # ----------------------------------------------------------------------

test_misc_K()
test_two_roads_and_islands(8,1,3)
test_two_roads_and_islands(19,3,2)

test_misc_F(False)
test_misc_F(True)
test_raster_raster_contact()
test_misc_A()
test_misc_B(False)
test_misc_B(True)
test_misc_C()

