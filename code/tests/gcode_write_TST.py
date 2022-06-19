#! /usr/bin/python3
# Test program for module {gcode_write}

import gcode_write
import path
import move 
import move_parms
import hacks
import path_example
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

zstep = parms['slice_thickness']
zspeed = parms['job_z_speed']
fdiam = parms['filament_diameter']
ex_temp = parms['nozzle_temperature']

pyx.unit.set(uscale=0.5, wscale=0.5, vscale=0.5)

def do_plot_slice(tag, islice, ph):
  # Plot the path {ph}. Writes files 
  # "tests/out/gcode_write_write_TST_{tag}.{ext}" where {ext} is "eps", "png", and "jpg".
  #
  # If {axes} is true, plots trace axes too. 
  #
  # If {matter} is true, shows the estimated actual area covered
  # by the extruded material of traces.

  ctrace =  pyx.color.rgb(0.050, 0.850, 0.000) # Color of nominal trace area.
  rwd = 0.80
  wd_axes =   0.05*wdf; 
  grid = True
  deco = True # ??? Should depend on size of slice ???
  fname = ("tests/out/gcode_write_TST_%s_s%04d" % (tag, islice))
  path.plot_to_files(fname, [ph,], [ctrace,], rwd, wd_axes, grid, deco)
  return 
  # ----------------------------------------------------------------------

def do_test_write_thing(tag, PHS):
  # Tests {gcode_write.thing} with the list {PHS} of slice tool-paths.
  # Writes files "tests/out/gcode_write_TST_{tag}_s{islice}.{ext}"
  # where {islice} is printed with 4 digits and {ext} is "eps", "png",
  # "jpg", and "gcode".

  # Plot the path for reference:
  for islice in range(len(PHS)):
    ph = PHS[islice]
    do_plot_slice(tag, islice, ph)

  # Generate the G-code:
  fname = ("tests/out/gcode_write_TST_%s.gcode" % tag)
  gwr = open(fname, "w", encoding="ISO-8859-1")
  gcode_write.thing(gwr, PHS, ex_temp, zstep, zspeed, fdiam, mp_jump)
  gwr.close()
  return 
  # ----------------------------------------------------------------------

def test_misc_F():
  # Generates G-code for a single-layer thing with a simple path.
  tag = "scfill"
  alt = True # Alternating diections.
  ph = path_example.misc_F(alt, mp_fill, mp_jump) 
  do_test_write_thing(tag, [ph,])
  return 
  # ----------------------------------------------------------------------

def test_gearloose():
  # Generates G-code for a single-layer thing with a medium complexity path.
  tag = "gearloose"
  R = 10
  ph = path_example.gearloose(R, True, mp_cont, mp_fill, mp_jump)
  do_test_write_thing(tag, [ph,])
  return 
  # ----------------------------------------------------------------------
 
test_misc_F()
test_gearloose()
