# Python3 procedures for linear algebra and Cartesian geometry.

MODULE_NAME = "rn"
MODULE_DESC = "Linear algebra operations on numeric vectors"
MODULE_VERS = "1.0"

MODULE_COPYRIGHT = "Copyright © 2009 State University of Campinas"

MODULE_INFO = \
  "A library module to perform linear algebra operations on numeric vectors.\n" \
  "\n" \
  "  Input vectors can be tuples or lists.  Output vectors will be tuples.\n"

import sys
import copy
from math import sqrt,sin,cos

# ------------------------------------------------------------------------
# GENERAL OPERATIONS ON POINTS AND VECTORS

def add(x,y) :
  # Vector sum of {x+y}.
  #
  n = len(x);
  assert len(y) == n, "incompatible {x,y} lenghts";
  r = [None]*n;
  for i in range(n) :
    r[i] = x[i] + y[i];
  return tuple(r);
  # ----------------------------------------------------------------------

def sub(x,y) :
  # Vector difference {x-y}.
  #
  n = len(x);
  assert len(y) == n, "incompatible {x,y} lenghts";
  r = [None]*n;
  for i in range(n) :
    r[i] = x[i] - y[i];
  return tuple(r);
  # ----------------------------------------------------------------------

def scale(s,x) :
  # Scales the vector {x} by {s}, which may be a float or a vector.
  #
  n = len(x);
  r = [None]*n;
  if type(s) is tuple or type(s) is list:
    assert len(s) == n, "incompatible {x,s} lengths"
    for i in range(n) :
      r[i] = s[i]*x[i];
  elif type(s) is int or type(s) is float:
    for i in range(n) :
      r[i] = s*x[i];
  else:
    assert False, "invalid scale {s}"
  return tuple(r);
  # ----------------------------------------------------------------------

def mix(s,x,t,y) :
  #
  # Returns {s*x+t*y}.
  n = len(x);
  assert len(y) == n, "incompatible {x,y} lenghts";
  r = [None]*n;
  for i in range(n) :
    r[i] = s*x[i] + t*y[i];
  return tuple(r);
  # ----------------------------------------------------------------------

def mix3(s,x,t,y,u,z) :
  # Returns {s*x+t*y+u*z}.
  #
  n = len(x);
  assert (len(y) == n) and (len(z) == n), "incompatible {x,y,z} lenghts";
  r = [None]*n;
  for i in range(n) :
    r[i] = s*x[i] + t*y[i] + u*z[i];
  return tuple(r);
  # ----------------------------------------------------------------------

def mix4(s,x,t,y,u,z,v,o) :
  # Returns {s*x+t*y+u*z+v*o}.
  #
  n = len(x);
  assert (len(y) == n) and (len(z) == n) and (len(o) == n), "incompatible {x,y,z,o} lenghts";
  r = [None]*n;
  for i in range(n) :
    r[i] = s*x[i] + t*y[i] + u*z[i] + v*o[i];
  return tuple(r);
  # ----------------------------------------------------------------------

def dir(x) :
  # Vector {x} normalized to unit Euclidean length. Also returns the original norm.
  #
  n = len(x);
  e = norm(x) + 1.0e-290;
  r = [None]*n;
  for i in range(n) :
    r[i] = x[i]/e;
  return tuple(r), e;
  # ----------------------------------------------------------------------

def dot(x,y) :
  # Scalar product of {x} by {y}.
  #
  n = len(x);
  assert len(y) == n, "incompatible {x,y} lenghts";
  s = 0;
  for i in range(n) :
    s += x[i] * y[i];
  return s;
  # ----------------------------------------------------------------------

def wdot(x,y,w):
  # Computes the weighted scalar product of vectors {x} and {y}
  # namely {\sum x[i]*y[i]*w[i]}.
  # However, if {w} is {None}, assumes {w[i]==1} for all {i}.
  n = len(x)
  s = 0
  for i in range(n):
    xyi = x[i] * y[i]
    if w != None: xyi = xyi * w[i]
    s += xyi
  return s
  # ----------------------------------------------------------------------

def norm_sqr(x) :
  # Square of Euclidean norm of {x}.
  #
  n = len(x);
  s = 0;
  for i in range(n) :
    xi=x[i]; s += xi*xi;
  return s;
  # ----------------------------------------------------------------------

def norm(x) :
  # Euclidean norm of {x}.
  #
  return sqrt(norm_sqr(x));
  # ----------------------------------------------------------------------

def dist(x,y) :
  # Euclidean distance between {x} and {y}.
  return norm(sub(x,y));
  # ----------------------------------------------------------------------

def dist_sqr(x,y) :
  # Euclidean distance squared between {x} and {y}.
  return norm_sqr(sub(x,y));
  # ----------------------------------------------------------------------

def pos_on_line(x,y,z):
  # Parameter {r} such that {mix(1-r,x,r,y) is closest to {z}.
  #
  vxz = sub(z, x)
  vxy = sub(y, x)
  dxy2 = norm_sqr(vxy)
  r = dot(vxz,vxy)/dxy2  # Relative position of {z}
  return r
  # ----------------------------------------------------------------------

def cross2d(x,y) :
  # Cross product of two vectors in R^2 (a real number).
  #
  assert len(x) == 2, "{x} must be a point of R^2";
  assert len(y) == 2, "{y} must be a point of R^2";
  return x[0]*y[1]-x[1]*y[0];
  # ----------------------------------------------------------------------

def cross3d(x,y) :
  # Cross product of two vectors in R^3 (a vector of R^3).
  #
  assert len(x) == 3, "{x} must be a point of R^3";
  assert len(y) == 3, "{y} must be a point of R^3";
  return ( x[1]*y[2]-x[2]*y[1], x[2]*y[0]-x[0]*y[2], x[0]*y[1]-x[1]*y[0] );
  # ----------------------------------------------------------------------

def rotate2(x,ang) :
  # Rotates the first two coords of {x} by {ang} radians around the origin
  #
  assert len(x) >= 2, "{x} must have at least 2 coords";
  c = cos(ang);
  s = sin(ang);
  return ( c*x[0] - s*x[1], s*x[0] + c*x[1] ) + tuple(x[2:]);
  # ----------------------------------------------------------------------
 
def write(wr, pref, x, fmt, sep, suff):
  # Writes to {wr} the string {pref}, then the coordinates of {x}, then the string {suff}.
  # Each coordinate is printed with float format {fmt}, and coords are separated by 
  # the string {sep}.
  if sep == None: sep = ' '
  if pref != None: wr.write(pref)
  for i in range(len(x)):
    if i > 0: wr.write(sep)
    wr.write(fmt % x[i])
  if suff != None: wr.write(suff)
  return
  # ----------------------------------------------------------------------
  
def poly_area(PS):
  # The argument {PS} should be a list of of vertices of a simple
  # polygon of {\RR^2}, in cyclic order. Returns the area of that
  # polygon.
  # 
  # A side from the last vertex to the first is implicitly assumed. The
  # sign of the result is {+1} if the vertices are in {CCW} order, {-1}
  # if in {CW} order. Result is undefined if the polygon has
  # self-intersections.
  
  np = len(PS)
  assert np > 0
  if np < 3: return 0
  
  A = 0
  o = PS[np-1]
  for i in range(np-1):
    p = PS[i]
    q = PS[(i+1) % np]
    if p != q:
      vp = sub(p, o)
      vq = sub(q, o)
      Ai = vp[0]*vq[1] - vp[1]*vq[0]
      A += Ai
  return A/2
  # ----------------------------------------------------------------------

# ------------------------------------------------------------------------
# SPECIAL-PURPOSE OPERATIONS

def pt_seg_dist(p, q0, q1):
  # Returns the Euclidean distance from point {p} to the segment {q0--q1}.
  # Also returns the point {q} on that segment that is closest to {p}
  u, dseg = dir(sub(q1, q0))
  if dseg == 0:
    # Segment has zero length:
    q = q0
    dpq = dist(p, q0)
  else:
    # Project {p} along {u}:
    r = dot(sub(p, q0), u);
    if r < 0:
      # Point is off the {q0} end:
      q = q0
      dpq = dist(p, q0)
    elif r > dseg:
      # Point is off the {q1} end:
      q = q1
      dpq = dist(p, q1)
    else:
      # Nearest point is interior to segment:
      q = mix(1.0, q0, r, u)
      dpq = dist(p, q)
  return dpq, q
  # ----------------------------------------------------------------------

def seg_seg_overlap(p0,p1, q0,q1, tol):
  # Checks whether the segments {p0--p1} and {q0--q1} are longer than
  # {tol} and parallel or antiparallel within tolerance {tol}, and their
  # projections on their mean direction overlap by a positive length.
  #
  # If those conditions are satisfied, returns four points {a0,a1} on
  # {p0--p1} and {b0,b1} on {q0--q1} which delimit the maximal subsets of
  # the two segments that, within tolerance {tol}, realize the minimum
  # distance between them.
  #
  # If those conditions are not satisfied, returns
  # {None,None,None,None}.
  
  # Displacement vecors of the two segments:
  p01 = sub(p1, p0)
  q01 = sub(q1, q0)

  # Check the lengths {dp,dq} of the two segments:
  dp = norm(p01)
  dq = norm(q01)
  if dp < tol or dq < tol: 
    # sys.stderr.write("segments too short\n")
    return None, None, None, None
    
  # Segments long enough; check the angle {ang} between them:
  c = dot(p01,q01)/(dp*dq)    # Cos(ang).
  s = sqrt(max(0, 1 - c**2))  # Abs(sin(ang)).
  # sys.stderr.write("p01 = ( %+.9f %+.9f ) q01 = ( %+.9f %+.9f )" % (p01[0],p01[1],q01[0],q01[1]))
  # sys.stderr.write(" sin(ang) = %.9f\n" % s)
  if s*min(dp, dq) > tol: 
    # sys.stderr.write("do not seem parallel/antiparallel\n")
    return None, None, None, None

  # Pretend that the segments are parallel or antiparallel.
  # sys.stderr.write("seem parallel or antiparallel\n")
  if c < 0:
    # Antipparallel -- flip {q0,q1} to make them parallel:
    q0,q1 = q1,q0; q01 = sub(q1, q0); c = -c;
  # Compute the unit vector of the mean direction:
  w, wm = dir(add(p01,q01))
  # Compute the segment lengths in that direction:
  dwp = dot(p01,w)
  dwq = dot(q01,w)
  # Compute the endpoint coordinates on the common direction, clipped:
  ra0 = min(dwp, max(0, dot(sub(q0,p0), w)))
  ra1 = min(dwp, max(0, dot(sub(q1,p0), w)))
  rb0 = min(dwq, max(0, dot(sub(p0,q0), w)))
  rb1 = min(dwq, max(0, dot(sub(p1,q0), w)))
  # sys.stderr.write("ra = [ %.9f _ %.9f ] rb =  [ %.9f _ %.9f ]\n" % (ra0/dwp,ra1/dwp,rb0/dwq,rb1/dwq))
  assert ra0 <= ra1 and rb0 <= rb1
  assert ra1 == dwp or rb1 == dwq
  if ra1 - ra0 < tol or rb1 - rb0 < tol:
    # sys.stderr.write("projections do not really overlap\n")
    return None, None, None,None

  # Non-trivial overlap. Get the low points:
  if ra0 == 0:
    a0 = p0; b0 = mix(1, q0, rb0/dwq, q01)
  elif rb0 == 0:
    a0 = mix(1, p0, ra0/dwp, p01); b0 = q0
  else:
    assert False
  # Get the high points:
  if ra1 == dwp:
    a1 = p1; b1 = mix(1, q0, rb1/dwq, q01)
  elif rb1 == dwq:
    a1 = mix(1, p0, ra1/dwp, p01); b1 = q1
  else:
    assert False
  return a0,a1, b0,b1
  # ----------------------------------------------------------------------

# ----------------------------------------------------------------------
# BOXES

# An {n}-dimensional /box/ is a subset of {R^n} that is the Cartesian
# product of {n} intervals. It is represented here as a pair (2-tuple)
# of points {(plo,phi)}, whose coordinates are respectively the low and
# high ends of those intervals.  A valid box must have {plo[i] <= phi[i]}
# for all{i}.

def box_from_point(x):
  # Retuns the box that contains the single point {x}.
  #
  return (tuple(x), tuple(x))
  # ----------------------------------------------------------------------
  
def box_inside(x, B):
  # {True} if the box {B} is not {None} and the point {x} is inside or on the boundary of {B}.
  # {False} otherwise.
  if B == None: return False
  n = len(x);
  xlo = B[0]; xhi = B[1]
  assert (len(xlo) == n) and (len(xhi) == n), "incompatible {B[i],x} lenghts";
  for i in range(n):
    if x[i] < xlo[i] or x[i] > xhi[i]: return False
  return True
  # ----------------------------------------------------------------------
  
def box_include_point(B,x):
  # Returns the smallest box that includes the box {B} and the point {x}.
  # However, if {x} is {None}, returns {B}; else, if {B} is {None}, returns 
  # a box with the single point {x}.
  #
  if x == None:
    return B
  elif B == None:
    return box_from_point(x)
  else:
    n = len(x);
    assert (len(B[0]) == n) and (len(B[1]) == n), "incompatible {B[i],x} lenghts";
    plo = [None]*n;
    phi = [None]*n
    for i in range(n):
      plo[i] = min(B[0][i], x[i])
      phi[i] = max(B[1][i], x[i])
    return (tuple(plo), tuple(phi),)
  # ----------------------------------------------------------------------
  
def box_join(A,B):
  # Returns the smallest box that includes the two boxes {A} and {B}.
  # However, if either box is {None}, returns the other box.
  #
  if A == None:
    return B
  elif B == None:
    return A
  else:
    n = len(A[0]);
    assert (len(A[1]) == n), "incompatible {A[0],A[1]} lenghts";
    assert (len(B[0]) == n) and (len(B[1]) == n), "incompatible {A[i],B[i]} lenghts";
    plo = [None]*n;
    phi = [None]*n
    for i in range(n):
      plo[i] = min(A[0][i], B[0][i])
      phi[i] = max(A[1][i], B[1][i])
    return (tuple(plo), tuple(phi))
  # ----------------------------------------------------------------------
 
def box_size(B):
  # Returns a tuple where element {i} is the extent of the box {B} along coordinate axis {i}.
  #
  n = len(B[0]);
  assert (len(B[1]) == n), "incompatible {B[0],B[0]} lenghts";
  s = [None]*n
  for i in range(n):
    s[i] = B[1][i] - B[0][i]
  return tuple(s)
  # ----------------------------------------------------------------------

def box_expand(B,dlo,dhi):
  # The parametr {B} should be a box in some space {\RR^d}, and the
  # parameters {dlo,dhi} should be vectors of {\RR^d}. Returns a copy of
  # {B} with the lower and upper corners displaced by those vectors --
  # outwards if positive, inwards if negative.
  #
  # The procedure fails if the displacements would result in
  # an invalid box. If {B} is {None}, returns {None}.
  #
  if B == None:
    return B
  else:
    n = len(B[0]);
    assert (len(B[1]) == n), "incompatible {B[0],B[1]} lenghts";
    assert (len(dlo) == n) and (len(dhi) == n), "incompatible {dlo,dhi} lenghts";
    plo = [None]*n;
    phi = [None]*n
    for i in range(n):
      plo[i] = B[0][i] - dlo[i]
      phi[i] = B[1][i] + dhi[i]
      assert plo[i] <= phi[i], "invalid displacements for this box"
    return (tuple(plo), tuple(phi))
  # ----------------------------------------------------------------------

# ----------------------------------------------------------------------
