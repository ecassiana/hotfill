#! /usr/bin/python3
# Test program for module {block_example}.
# Last edited on 2021-10-16 04:00:13 by stolfi

import block_example
import block
import path
import move
import move_parms
import contact
import hacks
import job_parms
import color
import rn
import pyx
import sys
from math import sqrt, sin, cos, floor, ceil, inf, nan, pi

parms = job_parms.typical_js()
parms['solid_raster_width'] = 1.00
parms['contour_trace_width'] = 0.50

mp_jump = move_parms.make_for_jumps(parms)
mp_cont = move_parms.make_for_contours(parms)
mp_fill = move_parms.make_for_fillings(parms)

wdf = move_parms.width(mp_fill)
wdc = move_parms.width(mp_cont)

def show_and_plot(tag, BCS, deco):
  # Prints the blocks {BCS} to {stderr}, validates them, and plots them.

  block.show_list(sys.stderr, "    ", BCS, None, True)

  # Collect and show all choices of all blocks, ignoring reversals:
  PHS = set()
  for bc in BCS:
    for ich in range(block.nchoices(bc)):
      phi, dri = path.unpack(block.choice(bc,ich))
      PHS.add(phi)
  path.show_list(sys.stderr, "    ", list(PHS), None, True, True)

  for bc in BCS: block.validate(bc)

  fname = "tests/out/block_example_TST_" + tag
  rwd = 0.80
  wd_axes = 0.05*wdf
  axes = deco
  dots = deco
  arrows = deco
  matter = True
  CLRS = hacks.trace_colors(len(BCS), None)
  block.plot_to_files(fname, BCS, CLRS, rwd, wd_axes, matter=True, links=True, contacts=True)
  return
  # ----------------------------------------------------------------------

def test_single_raster():
  sys.stderr.write("--- testing {single_raster} ---\n")
  tag = "single_raster"
  xlo = 2
  xhi = 4
  y = 2
  bc = block_example.single_raster(xlo, xhi, y, mp_fill)
  show_and_plot(tag, [bc,], True)
  return
  # ----------------------------------------------------------------------

def test_spiral_rectangle():
  sys.stderr.write("--- testing {spiral_rectangle} ---\n")
  tag = "spiral_rectangle"
  plo = (2,1)
  szx = 7*wdf; szy = 6*wdf
  bc = block_example.spiral_rectangle(plo, szx, szy, mp_fill)
  show_and_plot(tag, [bc,], True)
  return
  # ----------------------------------------------------------------------

def test_raster_rectangle(hor, ver, alt):
  sys.stderr.write("--- testing {rater_rectangle} hor = %s ver = %s ---\n" % (str(hor), str(ver)))
  tag = ("raster_rectangle_h%d_v%d_a%d" % (int(hor),int(ver),int(alt)))
  plo = (2,2)
  nx = 4
  ny = 6
  bc = block_example.raster_rectangle(plo, nx, ny, hor, ver, alt, mp_fill, mp_jump)
  show_and_plot(tag, [bc,], True)
  return
  # ----------------------------------------------------------------------

def test_onion():
  sys.stderr.write("--- testing {onion} ---\n")
  tag = "onion"
  nch = 6
  ctr = (4,4)
  Rc = 4
  Rf = 1
  bc = block_example.onion(nch, ctr, Rc, mp_cont, Rf, mp_fill, mp_jump)
  show_and_plot(tag, [bc,], True)
  return
  # ----------------------------------------------------------------------

def test_misc_A():
  sys.stderr.write("--- testing {misc_A} ---\n")
  tag = "misc_A"
  BCS, OPHS, TRS, JMS = block_example.misc_A(mp_fill, mp_jump)
  for bc in BCS: block.validate(bc)
  show_and_plot(tag, BCS, True)
  path.show_list(sys.stderr, "    ", OPHS, None, True, True)
  move.show_list(sys.stderr, "    ", TRS, None)
  move.show_list(sys.stderr, "    ", JMS, None)
  return
  # ----------------------------------------------------------------------

def test_misc_B():
  sys.stderr.write("--- testing {misc_B} ---\n")
  tag = "misc_B"
  BCS, OPHS, TRS, JMS = block_example.misc_B(mp_fill, mp_jump)
  show_and_plot(tag, BCS, True)
  path.show_list(sys.stderr, "    ", OPHS, None, True, True)
  move.show_list(sys.stderr, "    ", TRS, None)
  move.show_list(sys.stderr, "    ", JMS, None)
  return
  # ----------------------------------------------------------------------

def test_misc_C():
  sys.stderr.write("--- testing {misc_C} ---\n")
  tag = "misc_C"
  BCS = block_example.misc_C(mp_fill)
  show_and_plot(tag, BCS, True)
  return
  # ----------------------------------------------------------------------

def test_misc_D():
  sys.stderr.write("--- testing {misc_D} ---\n")
  tag = "misc_D"
  BCS, CTS = block_example.misc_D(mp_fill, True)
  show_and_plot(tag, BCS, True)
  return
  # ----------------------------------------------------------------------

def test_misc_E():
  sys.stderr.write("--- testing {misc_E} ---\n")
  tag = "misc_E"
  BCS = block_example.misc_E(mp_fill, mp_jump)
  show_and_plot(tag, BCS, True)
  return
  # ----------------------------------------------------------------------

def test_misc_G():
  sys.stderr.write("--- testing {misc_G} ---\n")
  tag = "misc_G"
  BCS, PHS, TRS0, TRS1 = block_example.misc_G(mp_cont, mp_fill, mp_jump)
  show_and_plot(tag, BCS, True)
  path.show_list(sys.stderr, "    ", PHS, None, True, True)
  move.show_list(sys.stderr, "    ", TRS0, None)
  move.show_list(sys.stderr, "    ", TRS1, None)
  return
  # ----------------------------------------------------------------------

test_single_raster()
test_raster_rectangle(True,False,False)
test_raster_rectangle(True,False,True)
test_raster_rectangle(False,True,True)
test_raster_rectangle(True,True,False)
test_spiral_rectangle()
test_onion()
test_misc_A()
test_misc_B()
test_misc_C()
test_misc_D()
test_misc_E()

test_misc_G()
