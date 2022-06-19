# Implementation of the module {hotfill}.

import scanline
import bandpath
import raster
import path
import contact
import move
import move_parms
import hacks
import rn
import sys
import time
from math import sqrt, sin, cos, log, exp, floor, ceil, inf, nan, pi

def solve(OPHS, mp_jump, alt):

  # Describe and check the input data:
  describe_and_check_input(sys.stderr, OPHS, mp_jump)

  nph = len(OPHS) # Number of raster elements.
  assert nph >= 1

  # Separate the rasters by scan-line:
  #xdir = (1,0)
  #ydir = (0,1)
  #ystep, yphase = raster.get_spacing_and_phase(OPHS, xdir, ydir)
  #SCS = raster.separate_by_scanline(OPHS, xdir, ydir, ystep, yphase)

  if alt:

    # Separate the rasters by scan-line:
    xdir = (1,0)
    ydir = (0,1)
    ystep, yphase = raster.get_spacing_and_phase(OPHS, xdir, ydir)
    SCS = raster.separate_by_scanline(OPHS, xdir, ydir, ystep, yphase)

    # Reverse order and orientation of raster elements of {OPHS} on alternate scanlines:
    PFSnew = [] # Reaster paths, reversed if needed.
    rev = False # True if next scan-line is to be reversed.
    for Si in SCS:
      # Now {Si} is a list of indices of rasters in one scanline.
      if rev:
        PFSi = [ path.rev(OPHS[jrs]) for jrs in Si ]
        PFSi.reverse()
      else:
        PFSi = [ OPHS[jrs] for jrs in Si ]
      rev = not rev
      PFSnew += PFSi
    OPHS = PFSnew

  # Assemble the path:
  use_links = alt
  ph = path.concat(OPHS, use_links, mp_jump)
  return ph, OPHS

def describe_and_check_input(wr, OPHS, mp_jump):
  wr.write("!! Not implemented !!\n")
  return
  # ----------------------------------------------------------------------

