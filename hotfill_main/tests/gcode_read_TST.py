#! /usr/bin/python3
# Test program for module {gcode_read}
# Last edited on 2021-10-14 12:13:40 by stolfi

import gcode_read
import path
import move 
import move_parms
import hacks
import job_parms
import rn
import pyx
import sys
from math import sqrt, sin, cos, floor, ceil, inf, nan, pi

parms = job_parms.typical_js()

mp_jump = move_parms.make_for_jumps(parms)
mp_cont = move_parms.make_for_contours(parms)
mp_fill = move_parms.make_for_fillings(parms)

wdf = move_parms.width(mp_fill)
wdc = move_parms.width(mp_cont)

def do_plot(name, ph):
  # Plot the path {ph}. Writes files 
  # "tests/out/gcode_read_TST_{name}.{ext}" where {ext} is "eps", "png", and "jpg".
  #
  # If {axes} is true, plots trace axes too. 
  #
  # If {matter} is true, shows the estimated actual area covered
  # by the extruded material of traces.

  ctrace =  pyx.color.rgb(0.050, 0.850, 0.000) # Color of nominal trace area.
  rwd = 0.80
  wd_axes =   0.05*wdf; 
  grid = True
  deco = False # ??? Should depend on size of slice ???
  fname = ("tests/out/gcode_read_TST_" + name)
  path.plot_to_files(fname, [ph,], [ctrace,], rwd, wd_axes, grid, deco)
  return 
  # ----------------------------------------------------------------------

def do_test_read(dir, name, fdiam, zstep, xyzpos, epos, eret, fan,ntemp,btemp):
  fname = "tests/in/" + dir + "/" + name + ".gcode"
  grd = open(fname, 'r', encoding="ISO-8859-1")
  ac = 300      # Acceleration.
  ud = 0.20     # Trace/jump transition penalty time.
  state = gcode_read.make_state(fdiam, zstep, ac, ud, fname)
  unit = 1
  pabs = True
  eabs = True
  sp = 40
  lnum = 0
  gcode_read.set_state(state, unit,pabs,eabs, xyzpos, epos, eret, sp, fan,ntemp,btemp, lnum)
  gcode_read.show_state(sys.stderr,state)
  ph = gcode_read.slice(grd, state)
  gcode_read.show_state(sys.stderr,state)
  grd.close()
  
  nsh = 5
  sys.stderr.write("first and last %d moves:\n" % nsh)
  nmv = path.nelems(ph)
  for kmv in range(nmv):
    if kmv < nsh or kmv >= nmv - nsh:
      omvk = path.elem(ph,kmv)
      move.show(sys.stderr, ("%6d " % kmv), omvk, "\n", 6)
  
  do_plot(name, ph)
  return 
  # ----------------------------------------------------------------------

def test_read_1():
  name = "springy_pen_holder_000068"
  dir = "g-code-samples"
  fan = 255
  ntemp = 200
  btemp = 60
  fdiam = 1.75
  zstep = 0.175
  xyzpos = (60.658, 49.019, 12.206)
  epos = 2060.30324
  eret = 4.5 
  do_test_read(dir, name, fdiam, zstep, xyzpos, epos, eret, fan,ntemp,btemp)
  return 
  # ----------------------------------------------------------------------

def test_read_2():
  name = "daisy_knife_holder_000014"
  dir = "g-code-samples"
  fan = 255
  ntemp = 200
  btemp = 60
  fdiam = 1.75
  zstep = 0.175
  xyzpos = (60.658, 49.019, 12.206)
  epos = 1941.01665
  eret = 1947.51665 - epos
  do_test_read(dir, name, fdiam, zstep, xyzpos, epos, eret, fan,ntemp,btemp)
  return 
  # ----------------------------------------------------------------------

def test_read_3(layer, xpos, ypos):
  name = ("rueda_escornabot_%06d" % layer)
  dir = "g-code-samples"
  fan = 0
  ntemp = 210
  btemp = 60
  fdiam = 1.75
  zstep = 0.200
  xyzpos = (xpos, ypos, (layer+1)*zstep)
  epos = -2.000
  eret = 2.000
  do_test_read(dir, name, fdiam, zstep, xyzpos, epos, eret, fan,ntemp,btemp)
  return 
  # ----------------------------------------------------------------------

def test_read_4(layer, xpos, ypos):
  name = ("c-part_%06d" % layer)
  dir = "g-code-samples"
  fan = 0
  ntemp = 210
  btemp = 60
  fdiam = 1.75
  zstep = 0.200
  xyzpos = (xpos, ypos, (layer+1)*zstep)
  epos = -2.000
  eret = 2.000
  do_test_read(dir, name, fdiam, zstep, xyzpos, epos, eret, fan,ntemp,btemp)
  return 
  # ----------------------------------------------------------------------

test_read_3(40, 247.651, 47.307)
test_read_3(50, 80.762, 61.812)
test_read_3(30, 247.547, 47.149)
test_read_3(5, 106.005, 60.569)
test_read_2()
test_read_1()
