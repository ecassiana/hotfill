#! /usr/bin/python3
# Test program for module {rootray_shape}

import rootray_shape
import rootray
import hacks
import pyx
import rn
import sys
from math import sqrt, sin, cos, floor, ceil, inf, nan, pi

def closenough(f1, f2):
  s1, T1 = f1
  s2, T2 = f2
  if s1 != s2: return False
  if len(T1) != len(T2): return False
  if rn.dist(T1, T2) > 1.0e-10: return False
  return True
  # ----------------------------------------------------------------------

def do_test_shape(S, tag):

  ylo= 0.125
  yhi = 15.0
  ystep = 0.25
  xlo = 0.0 
  xhi = 15.0
  
  B = ((xlo-3,ylo-3), (xhi+3,yhi+3))
  dp = (0, 0)
  frame = False
  grid = True
  c, szx, szy = hacks.make_canvas(hacks.round_box(B,0.5), dp, None, frame, grid, 1, 1)
  
  cline = pyx.color.rgb( 1.000, 0.200, 0.000 )
  
  y = ylo-2
  while y < yhi+2:
    p = (xlo, y - 2.0)
    q = (xhi, y + 2.0)
    fraw = rootray_shape.trace(p, q, S)
    f = rootray.intersection(fraw, (+1, [-0.125, +1.125]))
    fs, ft = f
    assert fs == +1
    nt = len(ft)
    assert nt % 2 == 0
    PT = [ rn.mix(1-t,p, t, q) for t in ft ]
    for kt in range(nt//2):
      pk = PT[2*kt]
      qk = PT[2*kt + 1]
      hacks.plot_line(c, pk, qk, dp, cline, 0.5*ystep, None)
    y = y + ystep
      
  hacks.write_plot(c, "tests/out/rootray_shape_TST_" + tag);
  return
  # ----------------------------------------------------------------------

def test_basics():
  sys.stderr.write("--- testing {prism,quadric,box,affine,union,intersection,complement} ---\n")

  # Builds a {RootRay_Shape} object with multiple parts.
  
  b = 11
  c = 12
  d = 13
  e = 14
  
  P1 = rootray_shape.prism([ (0,1), (9,1), (9,9), (0,9), (0,3), (7,3), (7,7), (8,7), (8,2), (0,2), (0,1), ])
  P2 = rootray_shape.prism([ (1,4), (1,7), (6,7), (6,4), (1,4), ])
  P3 = rootray_shape.prism([ (b,7), (e,7), (e,9), (b,9), (b,7), ])
  P4 = rootray_shape.prism([ (2,5), (3,5), (3,6), (2,6), (2,5), ])
  P5 = rootray_shape.prism([ (4,5), (5,5), (5,6), (4,6), (4,5), ])
  
  PA = rootray_shape.union([P2, P4, P5])
  PB = rootray_shape.intersection([P1, PA])
  PC = rootray_shape.union([PB, P3])
  
  C = rootray_shape.quadric(2, 2, -1) # Unit circle {X^2 + Y^2 - 1}.
  Adir1 = ((2,1), (1,2));               bdir1 = (+12, +3)
  Ainv1 = ((2/3,-1/3), (-1/3,2/3));     binv1 = ( -7, +2)
  E1 = rootray_shape.affine(C, Ainv1, binv1)  # Elliptic island at {bdir1}.
  
  Adir4 = ((  +7,  00), ( 00,   +7));     bdir4 = ( +7, +7)
  Ainv4 = ((+1/7,  00), ( 00, +1/7));     binv4 = ( -1, -1)
  E4 = rootray_shape.affine(C, Ainv4, binv4)  # Big circle at {bdir4}.
  
  BX = rootray_shape.box((c,6), (d,8))
  
  CC = rootray_shape.complement(C)
  Adir3 = (( +2/5,  -1/5), ( -1/5,  +2/5)); bdir3 = (+12,  +3)
  Ainv3 = ((+10/3,  +5/3), ( +5/3, +10/3)); binv3 = (-45, -30)
  E2 = rootray_shape.affine(CC, Ainv3, binv3)  # Elliptic hole at {bdir3}.
  
  EX = rootray_shape.intersection([E1, E2])

  H = rootray_shape.hspace(0)
  Adir2 = ((-2, -1), (+1,-2));          bdir2 = (+2, +1)
  Ainv2 = ((-2/5, +1/5), (-1/5, -2/5)); binv2 = (+1, 00)
  H1 = rootray_shape.affine(H, Ainv2, binv2) 

  U = rootray_shape.union([PC, BX, EX])
  
  # S = CC
  # S = E2
  # S = H1
  # S = H
  # S = E4
  # S = rootray_shape.intersection([E4, H])
  S = rootray_shape.intersection([U, H1])
  
  do_test_shape(S, "basics")
  return
  # ----------------------------------------------------------------------
  
def test_from_polygons():
  sys.stderr.write("--- testing {from_polygons} ---\n")

  # Builds a {RootRay_Shape} object from nested polygons.
  
  b = 11
  c = 12
  d = 13
  e = 14
  
  P1 = [ (0,1), (9,1), (9,9), (0,9), (0,3), (7,3), (7,7), (8,7), (8,2), (0,2), (0,1), ] # CCW
  P2 = [ (1,4), (1,7), (6,7), (6,4), (1,4), ] # CW
  P3 = [ (b,7), (e,7), (e,9), (b,9), (b,7), ] # CCW
  P4 = [ (2,5), (3,5), (3,6), (2,6), (2,5), ] # CCW
  P5 = [ (4,5), (5,5), (5,6), (4,6), (4,5), ] # CCW
  
  PTCS = [ P1, P2, P3, P4, P5 ]
  
  S = rootray_shape.from_polygons(PTCS, +1)
  
  do_test_shape(S, "from_polygons")
  return
  # ----------------------------------------------------------------------
  
test_basics()
test_from_polygons()
