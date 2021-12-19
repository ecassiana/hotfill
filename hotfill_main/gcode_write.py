# Tools to generate G-code files from tool-paths.
# Last edited on 2021-10-03 21:35:02 by stolfi

import gcode_write_IMP
import path
import move
import move_parms
import sys

# The following parameters are common to several functions:
#
#   {wr}      the G-code file, that should be opened
#             (with {open}) and closed (with {wr.close()}) by the caller.
#
#   {islice}  the slice index, counting from 0 at the bottom-most slice.
# 
#   {zslice}  the {Z} coordinate of the BOTTOM of the slice.
#
#   {zstep}   the slice thickness, also the height of the extruder nozzle
#             above the bottom of the layer while extruding. Also used
#             to compute the amount of filament needed for each trace.
#
#   {zspeed}  the maximum (cruise) speed to use when moving in {Z} between layers.
# 
#   {fdiam}   the filament's diameter.
#
#   {mp_jump} a {Move_Parms} object that defines the acceleration and cruise speed
#             to use for the initial jump to the start of each slice's path.
# 
# All dimensions are in millimeters. 
#
# Every path should be well-formed, meaning that the starting point of
# each element (trace or jump) is the final point of the previous
# element (if any). These procedures will not generate any jumps, except
# possibly at the start and end of each slice.

# FULL THING AT ONCE

def thing(wr, OPHS, ex_temp, zstep, zspeed, fdiam, mp_jump):
  # Writes to file {wr} the G-code to fabricate the thing whose slices
  # are the oriented paths in the list {OPHS}. Writes the file preamble
  # and postamble. The paths must include any desired raft, scaffolding,
  # etc.
  #
  # The preamble will set the extruder temperature to {ex_temp} and wait
  # for it to be reached. Assumes that the slices all have the same
  # thickness {zstep}, so that the BOTTOM of slice {islice} is at {Z =
  # islice*step}.
  #
  # The structure of the file will be:
  #
  #   file preamble
  #     slice 0 preamble
  #       tool-path of slice 0
  #     slice 0 postamble
  #     slice 1 preamble
  #       tool-path of slice 1
  #     slice 1 postamble
  #     ...
  #     slice {n-1} preamble
  #       tool-path of slice {n-1}
  #     slice {n-1} postamble
  #   file postamble
  # 
  gcode_write_IMP.thing(wr, OPHS, ex_temp, zstep, zspeed, fdiam, mp_jump)

# SLICE BY SLICE

def file_preamble(wr, ex_temp):
  # Writes to {wr} the preamble of a G-code file, before the first 
  # slice.
  # 
  # The preamble sets the G-code unit of length to mm and specifies
  # absolute coordinates for all three axes, except for the filament
  # position which is set to relative. Sends the head and platform to
  # the home position and sets the {X,Y} origins there.
  #
  # Also turns the print cooling fan off. Sets the nozzle temperature to
  # {ex_temp} and waits for it to be reached.
  #
  # Extrudes some filament in order to prime the nozzle, then retracts it
  # by 2mm.
  gcode_write_IMP.file_preamble(wr, ex_temp)

def slice(wr, oph, islice, zslice, zstep, zspeed, fdiam, mp_jump):
  # Writes to {wr} the G-code to fabricate a single slice of the thing,
  # described by the oriented path {oph}. It calls {slice_preamble},
  # then {tool_path}, then {slice_postamble}.
  #
  # The procedure will move the nozzle to {Z}-coordinate {zslice + zstep}
  # with max (cruise) {Z} speed {zspeed} before starting to fabricate the path{oph}.
  # 
  # The procedure assumes that the printer parameters like the unit of length and
  # the absolute/relative options have been set as specified in {thing} above.
  # It should normally be called only after {slice},
  # {file_preamble}, or {slice_postamble}, or the equivalent.
  gcode_write_IMP.slice(wr, oph, islice, zslice, zstep, zspeed, fdiam, mp_jump)

def file_postamble(wr):
  # Writes to file {wr} the G-code file postamble, ending the print job 
  # after the last slice.
  #
  # Assumes that the filament is retracted by 2 mm and leaves it there.
  # Moves the nozzle and platform to the home position and turns heating 
  # and fan off.
  gcode_write_IMP.file_postamble(wr)
  
# BY PARTS OF SLICES

def slice_preamble(wr, islice, zslice, zstep, zspeed):
  # Writes to file {wr} the G-code preamble for a new slice with index
  # {islice}.
  # 
  # This procedure generates commands to position nozzle at the the
  # proper {Z} height {zslice+zstep} with max speed {zspeed}, as in
  # {slice} above.
  # 
  # The preamble does not change the {X,Y} coordinates. It assumes that
  # the filament has been retracted by 2 mm and leaves it in that
  # position.
  #
  # The first line of the preamble will be a comment "(start of slice {islice})".
  # 
  # This procedure should normally be called only after {slice}, {file_preamble}, or
  # {slice_postamble}.
  gcode_write_IMP.slice_preamble(wr, islice, zslice, zstep)

def tool_path(wr, oph, zstep, fdiam, mp_jump):
  # Writes to {wr} the G-code to fabricate the tooi-path {oph}.
  #
  # Assumes that the nozzle is at some arbitrary {X,Y} coordinates and
  # at the correct {Z} to extrude the layer, with the filament retracted
  # by 2mm. At the end, the filament will be retracted by 2 mm.
  #
  # Should normally be called only after {slice_preamble}. The caller is
  # expected to write the slice postamble after the procedure returns.
  gcode_write_IMP.tool_path(wr, oph, zstep, fdiam)
 
def slice_postamble(wr, islice):
  # Writes to file {wr} the G-code pstamble to finalize the execution of a
  # slice with index {islice}. 
  #
  # Assumes that the filament is retracted by 2 mm and leaves it there.
  # Leaves the {X,Y,Z} coordinates unchanged.   The postamble will end  
  # with a comment "(end of slice {islice}).
  # 
  # This procedure should normally be called only after {tool-path}.
  gcode_write_IMP.slice_postamble(wr, islice)

# MOVE BY MOVE:

def jump(wr, q, mp):
  # Writes the G-code to move from the current position to point {q}
  # without extruding, with the nominal width and dynamics specified by
  # the {Move_Parms} object {mp}. Assumes that the filament is retracted
  # by 2mm.
  gcode_write_IMP.jump(wr, q, mpa)

def trace(wr, p, q, mp, zstep, fdiam):
  # Writes the G-code to extrude a trace from {p} (assumed to be the
  # current position) to point {q}, with the nominal width and dynamics
  # specified by the {Move_Parms} object {mp}. Assumes that the filament
  # is at the nozzle (not retracted), and leaves it so at the end.
  gcode_write_IMP.trace(wr, q, mp, zstep, fdiam)

# AUXILIARY PROCS

def compute_feed_length(p, q, wd, zstep, fdiam):
  # Returns the length of filament needed to extrude a trace of nominal width 
  # {wd} from {p} to {q}.
  return gcode_write_IMP.compute_feed_length(p, q, wd, zstep, fdiam)

