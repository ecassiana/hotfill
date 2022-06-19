# Implementation of module {rootray_cart}.

import rootray
import hacks
import rn
import sys
from math import sqrt, sin, cos, inf, nan, pi

def halfspace(p, q, k):
  d = len(p)
  assert len(q) == d, "incompatible ponint lengths"
  assert 0 <= k and k < d
  # Equation of the
  # Coefficients of equation of the hyperplane {r(t)[k] = A t + B}:
  pk = p[k]
  qk = q[k]
  A =  qk - pk
  B = pk
  if A == 0 or abs(B/A) > 1.0e+100:
    # Line is nearly parallel to the hyperplane:
    if B > 0:
      # Pretend whole line is outside:
      f = rootray.vacuum()
    elif B < 0:
      # Pretend whole line is inside:
      f = rootray.plenum()
    else:
      # Line practically lies on hyperplane. Pretend it enters at {t = 0}:
      f = (+1, [ 0 ])
  else:
    # Compute intersection
    t = -B/A
    # Decide whether {r(-oo)} is inside or outside:
    if A < 0:
      # {r(-oo)} is outside:
      f = (+1, [t])
    else:
      f = (-1, [t])
  return f

def slab(p, q, k):
  d = len(p)
  assert len(q) == d, "incompatible ponint lengths"
  assert 0 <= k and k < d
  pk = p[k]
  qk = q[k]
  if pk == qk:
    # Seg is parallel to slab:
    f = rootray.vacuum() if abs(pk) > 1 else rootray.plenum()
  else:
    f1 = halfspace((+pk - 1,), (+qk - 1,), 0)
    f2 = halfspace((-pk - 1,), (-qk - 1,), 0)
    f = rootray.intersection(f1, f2)
  return f

def cube(p, q):
  d = len(p)
  assert len(q) == d, "incompatible ponint lengths"
  f = rootray.plenum()
  for k in range(d):
    fk = slab(p, q, k)
    f = rootray.intersection(f, fk)
  return f

def box(p, q, ca, cb):
  d = len(p)
  assert len(q) == d, "incompatible point lengths"
  assert len(ca) == d, "incompatible corner lengths"
  assert len(cb) == d, "incompatible corner lengths"
  # sys.stderr.write("box p = %s q = %s\n" % (str(p),str(q)))
  f = rootray.plenum()
  for i in range(d):
    # Remap {p[i],q[i]} by the map that takes box to the signed unit cube:
    cmdi = (ca[i] + cb[i])/2
    magi = 2/(cb[i]-ca[i])
    p1i = magi*(p[i] - cmdi)
    q1i = magi*(q[i] - cmdi)
    # sys.stderr.write("box p1i = %s q1i = %s\n" % (str(p1i),str(q1i)))
    fi = slab((p1i,), (q1i,), 0)
    f = rootray.intersection(f, fi)
  return f
  
def quadric(p, q, n, m, bias):
  d = len(p)
  assert len(q) == d, "incompatible ponint lengths"
  assert p != q, "points must be distinct"
  n = min(n, len(p))
  m = min(m, len(p))
  assert 0 <= m and m <= n
  if m == 0 and bias < 0:
    f = rootray.plenum()
  elif m == n and bias > 0:
    f = rootray.vacuum()
  else:
    # Let {r(t)} be {t*p + (1-t)*q}. 
    # Determine the coefficients of the quadratic formula
    # {A t^2 + B t + C} for the characteristic function
    # {sum1 - sum2 + bias} where {sum1} is the sum of squares
    # of the first {m} coordinates and {sum2} is the sum of suqares
    # of the followinng {n-m} coordinates.
    v = rn.sub(q, p)
    
    vm = v[0:m]; vn = v[m:n]
    # sys.stderr.write("vm = %s vn = %s\n" % (str(vm), str(vn)))
    pm = p[0:m]; pn = p[m:n]
    A = rn.dot(vm,vm) - rn.dot(vn,vn)
    B = 2*(rn.dot(vm,pm) - rn.dot(vn,pn))
    C = rn.dot(pm,pm) - rn.dot(pn,pn) + bias
    # sys.stderr.write("equation = %.6f t^2 + %.6f t + %.6f = 0\n" % (A, B, C))
    if A == 0 or abs(A) < 1.0e-16:
      # Line is on the surface?
      if B == 0 or abs(B) < 1.0e-16:
        if C > 0:
          f = rootray.vacuum()
        elif C < 0:
          f = rootray.plenum()
        else:
          f = (+1, [0])
      else:
        s = -1 if B > 0 else +1
        f = (s, [ -C/B ])
    else:
      t0, t1 = hacks.real_quadratic_roots(A, B, C)
      # sys.stderr.write("t0 = %.6f t1 = %.6f\n" % (t0, t1))
      # val0 = A*t0*t0 + B*t0 + C 
      # val1 = A*t1*t1 + B*t1 + C 
      # sys.stderr.write("val0 = %23.15e t1 = %23.15e\n" % (val0, val1))
      if t0 == None and t1 == None:
        # Line does not cross the quadric:
        f = rootray.vacuum()
      else:
        assert t0 != None and t1 != None
        assert t0 < t1
        s = -1 if A < 0 else +1
        if t0 == -inf and t1 == +inf:
          # Roots are very far apart:
          f = (-s, [])
        elif t0 == -inf:
          # Line entrance is very far in the past:
          f = (-s, [t1])
        elif t1 == +inf:
          # Line exit is very far in the future:
          f = (+s, [t0])
        else: 
          f = (+s, [t0, t1])
    return f
  
def ball(p, q):
  d = len(p)
  f = quadric(p, q, d, d, -1)
  return f
  
def cylinder(p, q, m):
  d = len(p)
  assert 0 <= m and m <= d
  f= quadric(p, q, d, m, -1)
  return f 

def cone(p, q, m):
  d = len(p)
  assert 0 <= m and m <= d
  f= quadric(p, q, d, m, 0)
  return f 

def prism(p, q, PTS):
  p0 = p[0:2]
  q0 = q[0:2]
  # Collect the crossings {TS} in the order they occur on the polygon.
  # Each element of {TS} is a pair {(t,sgn)} where {t} is the 
  # time parameter along {p--q} and {sgn} is {+1} if the segment enters
  # the polygon, {-1} if it exits.

  eps = 1.0e-6 # Tolerance to detect ray-vertex and popy-poly intersections.

  sgn_inf = hacks.poly_orientation(PTS) 
  
  npt = len(PTS)
  TSS = []   # Crossing times and signs
  for ipt in range(npt-1):
    p1 = PTS[ipt]
    q1 = PTS[ipt+1]
    t0, t1, sgn = line_seg_intersect(p0,q0, p1,q1)
    if t0 != None:
      assert t1 != None and t1 >= 0 and t1 <= 1
      assert sgn == +1 or sgn == -1
      if abs(t1) < eps or abs(1 - t1) < eps:
        sys.stderr.write("!! warning - crossing very close to vertex\n")
      TSS.append((t0,sgn))
  
  # Sort the crossing times:
  TSS.sort(key = (lambda ts: ts[0]))
  
  # Check consistency, remove duplicates:
  nt = 0    # Number of crossings kept in {TSS}
  sgn_prev = None
  t_prev = -inf
  for t, sgn in TSS:
    # Check alternation of signs:
    t_prev, sgn_prev = (-inf, sgn_inf) if nt == 0 else TSS[nt-1]
    if abs(t - t_prev) < eps:
      sys.syderr.write("!! warning - crossings are very close to each other\n")
    if sgn == sgn_prev:
      sys.stderr.write("!! warning - crossing sign %+d inconsistent\n" % sgn)
    if t == t_prev:
      # Duplicate crossing: cancel them.
      assert nt > 0
      nt = nt-1
    else:
      # Append the crossing to the kept ones:
      TSS[nt] = (t, sgn); nt += 1
  # 
  TS = ( sgn_inf, [ TSS[i][0] for i in range(nt) ])
  return TS
  # ----------------------------------------------------------------------
  
def line_seg_intersect(p0,q0, p1,q1):
  # Computes the intersection of the straight line {p0--q0} and the segment {p1--q1}.
  #
  # The segment is assumed to be closed at the first point and open 
  # at the second point.
  #
  # If the segment crosses the line transversely, returns {t0,t1,sgn}
  # where {t0} and {t1} are the relative positions of the intersection
  # points along the two segments (floats in {[0 _ 1]}), and {sgn} is
  # {+1} if {p0--q0} crosses {p1--q1} in from the "out" side to the "in"
  # side, and {-1} otherwise.
  #
  # If the segment does not cross the line, or lies entirely on it
  # (within float roundoff tolerance) then returns {None,None,None}.
  
  dmin = 1.0e-5; # Bomb if points are closer than this.

  # Unit vectors:
  u0, d0 = rn.dir(rn.sub(q0,p0)); v0 = (-u0[1],u0[0])
  u1, d1 = rn.dir(rn.sub(q1,p1)); v1 = (-u1[1],u1[0])
  
  assert d0 > dmin, "points {p0,q0} are too close"
  assert d1 > dmin, "points {p1,q1} are too close"

  # Test {p1--q1} against line {p0--q0}:
  p1s0 = rn.dot(v0, rn.sub(p1,p0))
  q1s0 = rn.dot(v0, rn.sub(q1,p0))
  if (p1s0 == 0 and q1s0 == 0):
    return None, None, None
  if (p1s0 > 0 and q1s0 >= 0) or (p1s0 < 0 and q1s0 <= 0):
    return None, None, None
  t1 = p1s0/(p1s0 - q1s0)
  
  # Compute position of crossing on line {p0--q0}:  
  p0s1 = rn.dot(v1, rn.sub(p0,p1))
  q0s1 = rn.dot(v1, rn.sub(q0,p1))
  if (p0s1 == 0 and q0s1 == 0):
    return None, None, None
  assert p0s1 != q0s1
  t0 = p0s1/(p0s1 - q0s1)

  if u0[0]*u1[1]-u0[1]*u1[0] < 0:
    sgn = -1
  else:
    sgn = +1
  
  return t0, t1, sgn
  # ----------------------------------------------------------------------
