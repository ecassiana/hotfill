# Tools for tracing straight rays through simple solids in Cartesian space {\RR^d}.

import rootray_cart_IMP
import rootray

# See {rootray.py} for an explanation of root-rays.

# RAYTRACING PRIMITIVES

# In the procedure of this section, the trajectory {r} is defined by two
# points {p,q} of some Cartesian space {\RR^d}, as {r(t) = t p +
# (1-t)q}. Note that this trajectory spans the whole infinite straight
# line defined by the two points. Each procedure returns a root ray for
# {r} relative to some simple subset {F} of {\RR^n}. 
#
# The two points must be distinct tuples or lists of floats with the same length {d}.
# These shapes assume the Cartesian model of Euclidean space.
  
def halfspace(p, q, k):
  # The figure {F} is the halfspace whose interior is all points whose coordinate {k} 
  # is negative. The index{k} should be in {0..d-1}.
  return rootray_cart_IMP.halfspace(p, q, k)
  
def slab(p, q, k):
  # The figure {F} consists of the points whose coordinate {k} is in the range {(-1 _ +1)}.
  return rootray_cart_IMP.slab(p, q, k)
  
def cube(p, q):
  # The figure {F} is the axis-aligned {d}-dimensional hypercube of side 
  # 2 centered at the origin; that is, {[-1 _+1]^d}.
  return rootray_cart_IMP.cube(p, q)
  
def box(p, q, ca, cb):
  # The figure {F} is the axis-aligned {d}-dimensional box with opposite corners 
  # {ca,cb}.  The box must have non-zero size along all axes.
  return rootray_cart_IMP.box(p, q, ca, cb)

def quadric(p, q, n, m, bias):
  # The figure {F} is the quadratic form of {\RR^n} which is the sum of
  # squares of the first {m} coordinates, minus the sum of squares of
  # the following {n-m} coordinates, plus the given {bias}. Any
  # additional coordinates are ignored. If {p} and {q} have less than n
  # coordinates, they are implicitly extended with zeros.
  #
  # If {m = n} and {bias} is not negative, the quadric is equivalent to
  # {rootray.vacuum()} . If {m} is zero and {bias} is non-positive, it
  # is equivalent to {rootray.plenum()}.
 return rootray_cart_IMP.quadric(p, q, n, m, bias)

def cylinder(p, q, m):
  # The figure {F} is a cylinder centered at the origin
  # that is a ball of unit radius on the first {m} coordinates,
  # and unrestrited in the last {d-m}. For example, if {d=3} and {m=2}, then {F}
  # is the infinite cylinder in {\RR^3} radius 1 whose axis is the {Z} axis.
  # It is equivalent to {quadric(p, q, m, m, -1)}.
  #
  # If {m == 0}, this is equivalent to {rootray.plenum()}. If {m == 1},
  # this is equivalent to {slab(p,q,0)}. If {m == d}, this is equivalent
  # to {ball(p,q)}.
  return rootray_cart_IMP.cylinder(p, q, m)
  
def cone(p, q, m):
  # The figure {F} is the quadratic form of {\RR^d} which is the sum of squares
  # of the first {m} coordinates minus the sum of squares of the next {d-m}.
  # It is the same as {quadric(p, q, d, m, 0)}.
  #
  # For example, if {d=3} and {m=2}, {F} consists of the two oppositely directed infinite cones
  # with vertex at the origin and revolving around the {Z} axis, with walls sloping at 45 degrees.
  # If {d=2} and {m=1}, {F} is the 
  #
  # If {m == 0}, this is essentially equivalent to {plenum()}. If {m == d}, this is essentially
  # equivalent {vacuum}.
  return rootray_cart_IMP.cone(p, q, m)

def ball(p, q):
  # The figure {F} is the {d}-dimensional ball of unit radius centered 
  # at the origin.  It is equivalent to {quadric(p,q,d,d,-1)}.
  return rootray_cart_IMP.ball(p, q)

def prism(p, q, PTS):
  # The figure {F} is a prism consisting of all points of {\RR^d} whose 
  # first two coordinates lie inside the simple polygon of {\RR^2} 
  # whose vertices are the elements of the list {PTS}.
  #
  # The polygon is assumed to be oriented so that the ineterior is on the
  # left margin.  Namely, if the vertices are listed in CCW order, the 
  # interior is the bounded side; if listed in CW order, the interior
  # is the unbounded side.
  return rootray_cart_IMP.prism(p, q, PTS)
