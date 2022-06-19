# /usr/bin/python3
# Test program for module {raster_example}

import raster_example
import move
import move_parms
import path
import contact
import block
import input_data
import job_parms
import hacks
import pyx
import rn
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

def test_patch_array(cols, rows, nx, ny, islands):

  sys.stderr.write("--- testing {patch_array} ---\n")
  sys.stderr.write("cols =  %s\n" % cols)
  sys.stderr.write("rows = %s\n" % rows)
  sys.stderr.write("nx =  %s\n" % nx)
  sys.stderr.write("ny =  %s\n" % ny)
  sys.stderr.write("islands =  %s\n" % islands)
  
  split = True
  max_lines = 5
  
  tag = \
    "patch_array" + \
    ("_c%02d" % cols) + \
    ("_r%02d" % rows) + \
    ("_nx%02d" % nx) + \
    ("_ny%02d" % nx) + \
    ("_is%s" % "FT"[islands])
  
  mp_link = mp_fill
  OCRS, OPHS, OLKS, CTS = raster_example.patch_array(cols, rows, nx,ny, islands, mp_cont, mp_fill, mp_link)
  o = None

  # Plot:
  deco = True
  ydir = (0,1)
  fname = "tests/out/raster_example_TST_" + tag
  input_data.show_check_and_plot(fname, OCRS, OPHS, OLKS, CTS, o, ydir, wd_axes, deco)
  return
  # ----------------------------------------------------------------------

def test_rasters_A():
  sys.stderr.write("--- testing {rasters_A} ---\n")
  tag = ("rasters_A")

  xdir = ( cos(pi/6), sin(pi/6) )
  ydir = ( -xdir[1], +xdir[0] )
  yphase = 0.45
  eps = 0.15*wd_fill # Magnitude of perturbations.
  
  OCRS = []
  OPHS = raster_example.rasters_A(mp_fill, xdir, ydir, yphase, eps)
  OLKS = []
  CTS = []
  o = None

  deco = True
  ydir = (0,1)
  fname = "tests/out/raster_example_TST_" + tag
  input_data.show_check_and_plot(fname, OCRS, OPHS, OLKS, CTS, o, ydir, wd_axes, deco)
  # ----------------------------------------------------------------------

def test_rasters_B(nph):
  sys.stderr.write("--- testing {rasters_B} ---\n")
  tag = ("rasters_B")

  xdir = ( cos(pi/6), sin(pi/6) )
  ydir = ( -xdir[1], +xdir[0] )
  yphase = 0.45
  eps = 0.15*wd_fill # Magnitude of perturbations.
  
  OCRS = []
  TRS, PHS, LKS, CTS = raster_example.rasters_B(nph, xdir, mp_fill, mp_link)
  o = None

  deco = True
  ydir = (0,1)
  fname = "tests/out/raster_example_TST_" + tag
  input_data.show_check_and_plot(fname, OCRS, PHS, LKS, CTS, o, ydir, wd_axes, deco)
  # ----------------------------------------------------------------------

test_rasters_B(3)
test_rasters_A()

nx = 3
ny = 20
test_patch_array( 5, 2, nx, ny, True) 
test_patch_array( 7, 1, nx, ny, False) 
test_patch_array(10, 1, nx, ny, False) 
test_patch_array(13, 1, nx, ny, False) 

