# Implementation of the module {gcode_write}.
# Last edited on 2021-03-21 13:22:40 by jstolfi

import gcode_write
import path
import move
import move_parms
import rn
import re
import sys
from math import sqrt, sin, cos, floor, ceil, pi, inf, nan

def thing(wr, OPHS, ex_temp, zstep, zspeed, fdiam, mp_jump):
  file_preamble(wr, ex_temp)
  for islice in range(len(OPHS)):
    oph = OPHS[islice]
    zslice = islice * zstep # {Z} of BOTTOM of slice.
    slice(wr, oph, islice, zslice, zstep, zspeed, fdiam, mp_jump)
  file_postamble(wr)
  return
  # ----------------------------------------------------------------------
  
def file_preamble(wr, ex_temp):

  wr.write("M104 S%.0f ; set temperature, no wait" % ex_temp)
  wr.write("M107 ; print cooling fan off\n")
  wr.write("G21 ; set units to millimeters\n")
  wr.write("G28 ; move to home position\n")
  wr.write("G92 X0 Y0 E0 ; set coordinate origin\n")
  wr.write("G90 ; use absolute coordinates\n")
  wr.write("M83 ; use relative coordinates for the filament\n")
  wr.write("M109 S%.0f ; wait for temperature" % ex_temp)

  # Extrude 20 mm of filament at 100 mm/min without moving, to prime the nozzle:
  wr.write("G1 E20.0 F100.0 ; extrude\n")
  wr.write("G1 EPOS-2.00000 F24000 ; retract filament\n")
  wr.write("\n") 
  wr.flush()
  return
  # ----------------------------------------------------------------------
  
def slice(wr, oph, islice, zslice, zstep, zspeed, fdiam, mp_jump):
  slice_preamble(wr, islice, zslice, zstep, zspeed)
  tool_path(wr, oph, zstep, fdiam, mp_jump)
  slice_postamble(wr, islice)
  return 
  # ----------------------------------------------------------------------
 
def file_postamble(wr):
  wr.write("M107 ; print cooling fan off\n")
  wr.write("M104 S0 ; nozzle heater off\n")
  wr.write("G28 ; home all axes\n")
  wr.flush()
  return
  # ----------------------------------------------------------------------
  
def slice_preamble(wr, islice, zslice, zstep, zspeed):
  wr.write("(begin layer %d)\n" % islice)
  zex = zslice + zstep # Z cordinate of nozzle while extruding.
  # G-code speeds are in mm/min not mm/s
  wr.write("G1 Z%.3f F%.0f ; move to layer's Z\n" % (zex, 60*zspeed))
  wr.write("\n") 

def slice_postamble(wr, islice):
  # Writes the G-code to finalize the current slice.
  wr.write("(end layer %d)\n" % islice)
  wr.write("\n") 
  wr.write("\n") 
  wr.flush()

def tool_path(wr, oph, zstep, fdiam, mp_jump):
  # ??? Should condense jumps ???
  # ??? Should add the start drop at each jump-trace transitions ???
  
  retracted = True # False if filament is at the nozzle, true if it is retracted by 2 mm.
  p = path.pini(oph)
  # Jump to initial point of path.
  jump(wr, p, mp_jump)  

  pant = p
  for imv in range(path.nelems(oph)):
    # wr.write("(Move %d)\n" % imv)
    omvk = path.elem(oph, imv)
    mpk = move.parameters(omvk)
    p, q = move.endpoints(omvk)
    assert p == pant, "discontinuous path"
    if move.is_jump(omvk):
      if not retracted:  wr.write("G1 E-2 F2400\n")
      jump(wr, q, mpk)
      retracted = True
    else:
      if retracted: wr.write("G1 E2 F2400\n")
      trace(wr, p, q, mpk, zstep, fdiam)
      retracted = False
    pant = q
  
  # Make sure that the filament is retracted at the end:
  if not retracted: wr.write("G1 E-2.00000 F24000 ; retracts 2mm of material\n")
    
def jump(wr, q, mp):
  ac, sp, ud = move_parms.dynamics(mp)
  # Speed in G-code is mm/min not mm/s
  wr.write("G0 E0 X%.6f Y%.6f F%d\n" % (q[0], q[1], 60*sp))

def trace(wr, p, q, mp, zstep, fdiam):
  wd = move_parms.width(mp)
  ac, sp, ud = move_parms.dynamics(mp)
  tfeed = compute_feed_length(p, q, wd, zstep, fdiam)
  # Speed in G-code is mm/min not mm/s
  wr.write("G1 X%.6f Y%.6f E%.3f F%d\n" % (q[0], q[1], tfeed, 60*sp))
    
def compute_feed_length(p, q, wd, zstep, fdiam):
  dpq = rn.dist(p, q)

  # Compute area {tarea} of cross-section of extruded material, assuming that it is 
  # a rectangle with height {zstep} and width {wdt-zstep} with semicircles
  # of radius {zstep/2} on each side:
  tarea = zstep * wd

  volume = dpq*tarea # Volume of material needed
  
  # Compute area of cross-section of filament:
  farea = pi*fdiam*fdiam/4
  
  # Length of filament
  flength = volume/farea

  return flength

  
  
