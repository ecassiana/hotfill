#! /usr/bin/python3
# Test program for module {path}
# Last edited on 2021-10-17 17:26:08 by stolfi

import raster
import path
import rootray_shape
import raster_example
import path_example
import block_example
import move 
import move_parms
import input_data
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

def do_plot(tag, OCRS, BCS, o, links, contacts):

  CLRS = hacks.trace_colors(len(BCS), None)

  outname = "tests/out/input_data_TST_" + tag
  deco = True
  links = True
  contacts = True
  input_data.plot_input(outname, OCRS, BCS, o, CLRS, wd_axes,deco, links, contacts)
  return
  # ----------------------------------------------------------------------


def test_plot_input():
  sys.stderr.write("--- testing {plot_input}---\n")
  
  tag = "plot_input"

  OCRS = [] # ??? For now.  We need contours.
  BCS = block_example.misc_C(mp_fill)
  o = (0,0)
  
  do_plot(tag, OCRS, BCS, o, True, True)
  return
  # ----------------------------------------------------------------------

def test_make_synthetic(dsname, variant):
  sys.stderr.write("--- testing {make_synthetic}---\n")
  
  tag = "make_synthetic_" + dsname + "_" + variant
  
  OCRS, BCS, o, Z =  input_data.make_synthetic(dsname, variant, mp_cont, mp_fill, mp_link, mp_jump)

  do_plot(tag, OCRS, BCS, o, True, True)
  return
  # ----------------------------------------------------------------------

test_plot_input()

test_make_synthetic("paper_fig_B", "rasters")

test_make_synthetic("patch_array", "1x1x7no")

test_make_synthetic("patch_array", "1x1x27no")
test_make_synthetic("patch_array", "2x1x27no")
test_make_synthetic("patch_array", "3x1x27no")

test_make_synthetic("patch_array", "2x1x27is")
