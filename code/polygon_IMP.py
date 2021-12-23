# Implementation of module {polygon}
# Last edited on 2021-06-04 03:41:10 by jstolfi

import polygon

import path

import rn

class Polygon_IMP:
  def __init__(self, PTS):
    self.PTS = tuple(PTS) # List of vertices of the polygon, implicitly closed.

def from_points(PTS):
  return polygon.Polygon(PTS)
  # ----------------------------------------------------------------------
  
def from_point_lists(PTSS):
  PGS = []
  for PTS in PTSS:
    pg = from_points(PTS)   
    PGS.append(pg)
  return PGS
  # ----------------------------------------------------------------------
  
def make_path(pg, mp_cont):
  PTS = list(pg.PTS).copy()
  PTS.append(PTS[0])
  ph = path.from_points(PTS, mp_cont, None)
  return ph
  # ----------------------------------------------------------------------

def get_points(pg):
  return pg.PTS
  # ----------------------------------------------------------------------

def bbox(PGS):
  B = None
  for pg in PGS:
    for pt in pg.PTS:
      B = rn.box_include_point(B, pt)
  return B
  # ----------------------------------------------------------------------
