#! /usr/bin/python3

import polygon
import path
import move
import move_parms
import job_parms
import rn
import sys

parms = job_parms.typical_js()
parms['solid_raster_width'] = 1.0
parms['contour_trace_width'] = 0.5

mp_jump = move_parms.make_for_jumps(parms)
mp_cont = move_parms.make_for_contours(parms)
mp_fill = move_parms.make_for_fillings(parms)

def test_basic():
  sys.stderr.write("--- testing {from_points,get_points,make_path,bbox} ---\n")
  
  p10 = ( 1, 1)
  p11 = ( 5, 2)
  p12 = ( 4, 4)
  p13 = ( 2, 3)
  PTS1 = (p10,p11,p12,p13,)
  npt1 = len(PTS1)
  
  pg1 = polygon.from_points(PTS1)
  assert isinstance(pg1, polygon.Polygon)
  
  PTSx = polygon.get_points(pg1)
  assert tuple(PTS1) == tuple(PTSx)
  
  # Chech {make_path}:
  ph = polygon.make_path(pg1, mp_cont)
  assert isinstance(ph, path.Path)
  assert path.nelems(ph) == len(PTS1) 
  assert path.pini(ph) == PTS1[0]
  assert path.pfin(ph) == PTS1[0]
  for ipt in range(npt1):
    pt1 = PTS1[ipt % npt1]
    ptx = move.pini(path.elem(ph, ipt))
    assert pt1 == ptx

  # Another polygon:
  p20 = ( 6, 3)
  p21 = (11, 4)
  p22 = (10, 6)
  p23 = ( 7, 5)
  PTS2 = (p20,p21,p22,p23,)
  pg2 = polygon.from_points(PTS2)
  
  PGS = [ pg1,pg2, ]
  
  # Check {bbox}:
  bb1 = polygon.bbox(PGS)
  bbx = None
  for pt in PTS1 + PTS2: bbx = rn.box_include_point(bbx, pt)
  assert bb1 == bbx
  
  return
 
test_basic()
