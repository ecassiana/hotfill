#! /usr/bin/python3
# Test program for module {rootray_cart}

import rootray_cart
import rootray
import hacks
import rn
import sys
from math import sqrt, sin, cos, floor, ceil, inf, nan, pi

def checksol(name, f_cmp, f_exp):
  sys.stderr.write("%s_cmp = %s\n" % (name, str(f_cmp)))
  sys.stderr.write("%s_exp = %s\n" % (name, str(f_exp)))
  assert closenough (f_cmp, f_exp)
  return
  # ----------------------------------------------------------------------

def closenough(f1, f2):
  s1, T1 = f1
  s2, T2 = f2
  if s1 != s2: return False
  if len(T1) != len(T2): return False
  if rn.dist(T1, T2) > 1.0e-10: return False
  return True
  # ----------------------------------------------------------------------

def test_halfspace():
  sys.stderr.write("--- testing {halfspace} ---\n")

  p1 = (-1, -2, -3, -4)
  q1 = (+4, +5, +6, +7)
  f1_cmp = rootray_cart.halfspace(p1, q1, 2)
  f1_exp = (-1, [ 1.0/3.0 ])
  checksol("f1", f1_cmp, f1_exp)
  return
  # ----------------------------------------------------------------------

def test_slab():
  sys.stderr.write("--- testing {slab} ---\n")

  p1 = (-1, -2, -3, -4)
  q1 = (+4, +5, +6, +7)
  f1_cmp = rootray_cart.slab(p1, q1, 2)
  f1_exp = (+1, [ 2.0/9.0, 4.0/9.0 ])
  checksol("f1", f1_cmp, f1_exp)
  return
  # ----------------------------------------------------------------------

def test_cube():
  sys.stderr.write("--- testing {cube} ---\n")

  p1 = (-1, -2, -3, -4)
  q1 = (+4, +5, +6, +7)
  f1_cmp = rootray_cart.cube(p1, q1)
  f1_exp = (+1, [ 3.0/11.0, 2.0/5.0 ])
  checksol("f1", f1_cmp, f1_exp)
  return
  # ----------------------------------------------------------------------

def test_box():
  sys.stderr.write("--- testing {box} ---\n")

  p1 = (-1, -2, -3, -4)
  q1 = (+4, +5, +6, +7)
  ca = (00, 00, -2, -3)
  cb = (+3, +4, 00, +6)
  f1_cmp = rootray_cart.box(p1, q1, ca, cb)
  f1_exp = (+1, [ 2.0/7.0, 3.0/9.0 ])
  checksol("f1", f1_cmp, f1_exp)
  return
  # ----------------------------------------------------------------------

def test_quadric():
  sys.stderr.write("--- testing {quadric} ---\n")

  # Hyperbola:
  p1 = (-4, -1, +2)
  q1 = (+4, +1, +7)
  n1 = 2
  m1 = 1
  v1 = -4
  sys.stderr.write("n = %d m= %d bias v = %.8f\n" % (n1, m1, v1))
  f1_cmp = rootray_cart.quadric(p1, q1, n1, m1, v1)
  t1a = 0.5*1 - 1/sqrt(15); t1b = 0.5*1 + 1/sqrt(15)
  f1_exp = (+1, [ t1a, t1b ])
  checksol("f1", f1_cmp, f1_exp)

  # Ball:
  p2 = (-4, 00, +2)
  q2 = (+4, +1, +7)
  n2 = 2
  m2 = 2
  v2 = -4
  sys.stderr.write("n = %d m= %d bias v = %.8f\n" % (n2, m2, v2))
  f2_cmp = rootray_cart.quadric(p2, q2, n2, m2, v2)
  sdel = 2*sqrt(61)
  t2a = (32 - sdel)/65; t2b = (32 + sdel)/65
  f2_exp = (+1, [ t2a, t2b ])
  checksol("f2", f2_cmp, f2_exp)
 
  # Ball hole:
  p3 = (-4, 00, +2)
  q3 = (+4, +1, +7)
  n3 = 2
  m3 = 0
  v3 = +4
  sys.stderr.write("n = %d m= %d bias v = %.8f\n" % (n3, m3, v3))
  f3_cmp = rootray_cart.quadric(p3, q3, n3, m3, v3)
  sdel = 2*sqrt(61)
  t3a = (32 - sdel)/65; t3b = (32 + sdel)/65
  f3_exp = (-1, [ t3a, t3b ])
  checksol("f3", f3_cmp, f3_exp)
 
  return
  # ----------------------------------------------------------------------

def test_ball():
  sys.stderr.write("--- testing {ball} ---\n")

  p1 = (00, 00, 00)
  q1 = (+2, +2, +1)
  f1_cmp = rootray_cart.ball(p1, q1)
  f1_exp = (+1, [ -1.0/3.0, +1.0/3.0 ])
  checksol("f1", f1_cmp, f1_exp)

  p2 = (-3.00, 00.00, +0.50)
  q2 = (+6.00, 00.00, +0.50)
  f2_cmp = rootray_cart.ball(p2, q2)
  cc2a = cos(pi/6)
  t2a = (-cc2a - p2[0])/(q2[0] - p2[0])
  t2b = (+cc2a - p2[0])/(q2[0] - p2[0])
  f2_exp = (+1, [ t2a, t2b ])
  checksol("f2", f2_cmp, f2_exp)
  return
  # ----------------------------------------------------------------------

def test_cone():
  sys.stderr.write("--- testing {cone} ---\n")

  p1 = (+0.5, 00.0, +1.0)
  q1 = (00.0, +0.5, -1.0)
  f1_cmp = rootray_cart.cone(p1, q1, 2)
  t1a = 0.5*(1 - 1/sqrt(7)); t1b = 0.5*(1 + 1/sqrt(7))
  f1_exp = (-1, [ t1a, t1b ])
  checksol("f1", f1_cmp, f1_exp)

  p2 = (+1.00, 00.00, +0.50)
  q2 = (00.00, +0.50, 00.00)
  f2_cmp = rootray_cart.cone(p2, q2, 1)
  t2a = 0.5*(3 - sqrt(3)); t2b = 0.5*(3 + sqrt(3))
  f2_exp = (+1, [ t2a, t2b ])
  checksol("f2", f2_cmp, f2_exp)
  return
  # ----------------------------------------------------------------------

def test_prism():
  sys.stderr.write("--- testing {prism} ---\n") 

  PTS = [ (00, -3), (+5, -3), (+5, -1), (+1, -1), (+1, +2), (00, +1), (00, -3), ]

  p1 = (-0.50, +0.50, +1.00)
  q1 = (+1.50, +0.50, -2.00)
  f1_cmp = rootray_cart.prism(p1, q1, PTS)
  t1a = +0.25; t1b = +0.75
  f1_exp = (+1, [ t1a, t1b ])
  checksol("f1", f1_cmp, f1_exp)

  p2 = (+0.5, +1.0)
  q2 = (00.0, +2.0)
  f2_cmp = rootray_cart.prism(p2, q2, PTS)
  t2a = -4.0; t2b = -2.0; t2c = -1.0; t2d = +1/3
  f2_exp = (+1, [ t2a, t2b, t2c, t2d ])
  checksol("f2", f2_cmp, f2_exp)
  return
  # ----------------------------------------------------------------------

test_halfspace()
test_slab()
test_cube()
test_box()
test_quadric()
test_ball()
test_cone()
test_prism()
