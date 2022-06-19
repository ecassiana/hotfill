#! /usr/bin/python3
# Test program for module {move_example}

import move_example
import move
import move_parms
import hacks
import job_parms
import rn
import pyx
import sys
from math import sqrt, sin, cos, floor, ceil, inf, nan, pi

parms = job_parms.slow()

mp_fill = move_parms.make_for_fillings(parms)
mp_cont = move_parms.make_for_contours(parms)
mp_jump = move_parms.make_for_jumps(parms)

wd_fill = move_parms.width(mp_fill)
wd_cont = move_parms.width(mp_cont)

def plot(tag, OMVS):
  # Plots the moves in the list {OMVS} to files.
  nmv = len(OMVS)
  CLRS = hacks.trace_colors(nmv, 0.300)
  fname = ("tests/out/move_example_TST_%s" % tag)
  rwd = 0.80
  wd_axes = 0.15*min(wd_fill,wd_cont) # Reference line width.
  move.plot_to_files(fname, OMVS, CLRS, rwd, wd_axes)
  return 
  # ----------------------------------------------------------------------

def test_rectangle_rasters(axis):
  sys.stderr.write("--- testing {rectangle_raster,rectangle_links,rectangle_jumps} axis = %d ---\n" % int(axis))
  plo = (1,2)
  step = wd_fill
  nx = 3; ny = 5
  if axis == 0:
    sz = (nx-1)*step; n = ny
  else:
    sz = (ny-1)*step; n = nx
  TRS = move_example.rectangle_rasters(plo, axis, n, sz, step, mp_fill)
  LKS0, LKS1 = move_example.rectangle_links(plo, axis, n, sz, step, mp_fill)
  JMS0, JMS1 = move_example.rectangle_jumps(plo, axis, n, sz, step, mp_jump)
  
  move.show_list(sys.stderr, "    ", TRS, None)
  move.show_list(sys.stderr, "    ", LKS0, None)
  move.show_list(sys.stderr, "    ", LKS1, None)
  move.show_list(sys.stderr, "    ", JMS0, None)
  move.show_list(sys.stderr, "    ", JMS1, None)

  tag = ("rectangle_raster_ax%d" % axis)
  plot(tag, TRS+LKS0+LKS1+JMS0+JMS1)
  return
  # ----------------------------------------------------------------------

def test_misc_A():
  sys.stderr.write("--- testing {misc_A} ---\n")
  tag = "misc_A"
  MVS = move_example.misc_A(mp_fill, mp_cont, mp_jump)

  move.show_list(sys.stderr, "    ", MVS, None)

  plot(tag, MVS)
  return
  # ----------------------------------------------------------------------

def test_misc_B():
  sys.stderr.write("--- testing {misc_B} ---\n")
  tag = "misc_B"
  TRS, JMS = move_example.misc_B(mp_fill, mp_jump)
  
  move.show_list(sys.stderr, "    ", TRS, None)
  move.show_list(sys.stderr, "    ", JMS, None)
  
  plot(tag, TRS+JMS)
  return
  # ----------------------------------------------------------------------
 
def test_misc_C():
  sys.stderr.write("--- testing {misc_C} ---\n")
  tag = "misc_C"
  TRS, JMS = move_example.misc_C(mp_fill, mp_jump)
  
  move.show_list(sys.stderr, "    ", TRS, None)
  move.show_list(sys.stderr, "    ", JMS, None)
  
  plot(tag, TRS+JMS)
  return
  # ----------------------------------------------------------------------
  
# Run the tests:

test_rectangle_rasters(0)
test_rectangle_rasters(1)
test_misc_A()
test_misc_B()
test_misc_C()
 
