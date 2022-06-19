# Tools for ray-tracing csg models in arbitrary spaces with arbitrary probing trajectories.

import rootray_shape_IMP; from rootray_shape_IMP import RootRay_Shape_IMP

class RootRay_Shape(RootRay_Shape_IMP):
  # A {RootRay_Shape} object represents a complex shape obtained from 
  # primitive objects by boolean operations and affine transformations.
  # 
  # A {RootRay_Shape} object has an attribute {op}, a string that 
  # determines the type of obejct.
  pass
  
# PRIMITIVES

def hspace(m):
  # Creates an object with {op} ="hspace" and an attribute {m}. The
  # characteristic function is coordinate {m} of the point.
  return rootray_shape_IMP.hspace(m)

def box(ca, cb):
  # Creates an object with {op} ="box" and an attributes {ca,cb}.  The points {ca,cb} must have 
  # the same number {m} of coordinates.  The
  # characteristic function is positive if the first {m} coordinaes
  # of the point are inside the box of {\RR^m} which has {ca} and {cb}
  # as opposite corners.
  return rootray_shape_IMP.box(ca, cb)
  
def quadric(n, m, bias):
  # Creates an object with {op} = "quadric" and attributes {n}, {m}, and
  # {bias}. Its characteristic function is the sum of the squares of the
  # first {m} coordinates minus the squares of the next {n-m}
  # coordinates, plus the {bias}.
  return rootray_shape_IMP.quadric(n, m, bias)

def prism(PTS):
  # Creates an object with {op} = "prism" object has a list of points of
  # {\RR^2}. See {rootray_cart.prism} for the semantics.
  return rootray_shape_IMP.prism(PTS)

# TRANSFORM

def affine(arg, Ainv, binv):
  # Creates an object with {op} = "affine". The attributes are: a child {RootRay_Shape} object {arg}, a {d} by
  # {d} non-singular array {Ainv}, and a {d}-vector {binv}.
  #
  # It represents the shape represented by {arg}, mapped by a linerar map
  # {Adir} and then translated by {bdir}. The attributes {Ainv} and {binv}
  # are the analogous parameters of the inverse affine map; that is,
  # {Ainv=Adir^{-1}} and {binv = bdir*Ainv}.
  return rootray_shape_IMP.affine(arg, Ainv, binv)
  
# BOOLEAN

def union(ARGS):
  # Creates an object with {op} = "union" object, whose attribute is a
  # list {ARGS} (possibly empty) of other {RootRay_Shape} objects. It
  # represents the union of the shapes defined by those objects. If the
  # list is empty, it describes the vacum shape.
  return rootray_shape_IMP.union(ARGS)
  
def intersection(ARGS):
  # Creates an object with {op} = "intersection" object, whose attribute is a
  # list {ARGS} (possibly empty) of other {RootRay_Shape} objects. It
  # represents the intersection of the shapes defined by those objects. If the
  # list is empty, it describes the plenum shape.
  return rootray_shape_IMP.intersection(ARGS)

def complement(arg):
  # Creates an object with {op} = "complement" object with one attribute
  # {arg}, a {RootRay_Shape} object. It represents the set-complement
  # with respect to {\RR^n} of the shape described by {arg}.
  return rootray_shape_IMP.complement(arg)

# TRACING

def trace(p, q, obj):
  # Traces the straight line through the points {p,q} and returns 
  # a root-ray {(s,T)} that describes the intersections of that line 
  # with the spahe defined by the {RootRay_Shape} object {obj}.
  # See the {rootray.py} interface for more details.
  return rootray_shape_IMP.trace(p, q, obj)
  
# SPECIAL
  
def from_polygons(PTCS, sgn_inf):
  # Assumes that {PTCS} is a list of lists, where each element is a
  # list of points that describe a simple closed polygon, with
  # the last point repeating the first point.  The boundaries of the polygons should be
  # pairwise disjoint, and nested and oriented so that the nearest enclosing polygon
  # of each polygon should be oriented in opposite directions.  Specifically,
  # island vertices should be in CCW order, and hole vertices should be in CW order.
  #
  # The parameter {sgn_inf} should be {+1} the outermost polygons are islands,
  # and {-1} if they are holes.
  #
  # Returns a {RootRay_Shape} object that describes the interior of that
  # collection of polygons.
  return rootray_shape_IMP.from_polygons(PTCS, sgn_inf)
