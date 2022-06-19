# Implementation of module {rootray_shape}

import rootray_shape
import rootray_cart
import rootray
import hacks
import rn
import rmxn
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
import pyx 
from math import nan, inf, sqrt
import sys

class RootRay_Shape_IMP:
  
  def __init__(self, op):
    # Read-only fields:
    self.op = op
    self.n = None # Total relevant coords, for "quadric"
    self.m = None # Positive coord, for "quadric", and coord index, for "hplane".
    self.bias = None # Constant term of quadric.
    self.ARGS = None # Arguments for boolean and "affine", points for "prism".
    self.Ainv = None # Inverse linear matrix for "affine".
    self.binv = None # Inverse displacement for "affine". 
    self.ca = None   # Corner of box.
    self.cb = None   # Corner of box.

def hspace(m):
  obj = rootray_shape.RootRay_Shape("hspace")
  obj.m = m
  return obj
  # ----------------------------------------------------------------------
     
def box(ca, cb):
  assert len(ca) == len(cb), "wrong dimensions of {ca,cb}"
  obj = rootray_shape.RootRay_Shape("box")
  obj.ca = ca
  obj.cb = cb
  return obj
  # ----------------------------------------------------------------------
 
def quadric(n, m, bias):
  obj = rootray_shape.RootRay_Shape("quadric")
  obj.n = n
  obj.m = m
  obj.bias = bias
  return obj
  # ----------------------------------------------------------------------

def prism(PTS):
  obj = rootray_shape.RootRay_Shape("prism")
  assert len(PTS) >= 3, "too few points in polygon"
  obj.ARGS = tuple(PTS)
  return obj
  # ----------------------------------------------------------------------

def affine(arg, Ainv, binv):
  obj = rootray_shape.RootRay_Shape("affine")
  obj.ARGS = (arg,)
  obj.Ainv = tuple(Ainv)
  obj.binv = tuple(binv)
  return obj
  # ----------------------------------------------------------------------
  
def union(ARGS):
  obj = rootray_shape.RootRay_Shape("union")
  obj.ARGS = tuple(ARGS)
  return obj
  # ----------------------------------------------------------------------

def intersection(ARGS):
  obj = rootray_shape.RootRay_Shape("intersection")
  obj.ARGS = tuple(ARGS)
  return obj
  # ----------------------------------------------------------------------

def complement(arg):
  obj = rootray_shape.RootRay_Shape("complement")
  obj.ARGS = (arg,)
  return obj
  # ----------------------------------------------------------------------

def trace(p, q, obj):
  if obj.op == "hspace":
    f = rootray_cart.halfspace(p, q, obj.m)
  elif obj.op == "box":
    f = rootray_cart.box(p, q, obj.ca, obj.cb)
  elif obj.op == "quadric":
    f = rootray_cart.quadric(p, q, obj.n, obj.m, obj.bias)
  elif obj.op == "prism":
    f = rootray_cart.prism(p, q, obj.ARGS)
  elif obj.op == "affine":
    d = len(p)
    assert len(q) == d, "inconsistent point dimensions"
    p1 = rn.add(rmxn.map_row(p, obj.Ainv), obj.binv)
    q1 = rn.add(rmxn.map_row(q, obj.Ainv), obj.binv)
    f = trace(p1, q1, obj.ARGS[0])
  elif obj.op == "union":
    f = rootray.vacuum()
    for arg in obj.ARGS:
      g = trace(p, q, arg)
      f = rootray.union(f, g)
  elif obj.op == "intersection":
    f = rootray.plenum()
    for arg in obj.ARGS:
      g = trace(p, q, arg)
      f = rootray.intersection(f, g)
  elif obj.op == "complement":
    g = trace(p, q, obj.ARGS[0])
    f = rootray.complement(g) 
  else:
    assert False, "invalid op"
  return f
  # ----------------------------------------------------------------------

def from_polygons(POLYS, sgn_inf):
 
  npo = len(POLYS)
  SHPOLYS = [None]*npo   # Contour paths converted to {shapely} polygons.
  
  # Convert all point lists to {shapely} polygons:
  for ipo in range(npo):
    PTS = POLYS[ipo]
    assert PTS[0] == PTS[-1]
    SHPOLYS[ipo] = Polygon(PTS)
  
  # Find the parent of each contour:
  # This is quadratic on {npo} and probably proportional to the number of points
  PARENT = [None]*npo # {POLYS[PARENT[ipo]]} is the parent (nearest enclosing polygon) of {POLYS[ipo]}.
  for ipo1 in range(npo):
    # Find the index {ipar} of the parent of {POLYS[ipo1]}:
    ipar = None
    pg1 = SHPOLYS[ipo1]
    for ipo2 in range(npo):
      pg2 = SHPOLYS[ipo2]
      if ipo1 != ipo2:
        # See if {ipo2} can be the parent of {ipo1}:
        if pg2.contains(pg1):
          if ipar == None:
            # First container of {ipo1}, keep it for now:
            ipar = ipo2
          else:
            # We already have a container {ipar}. Check which is closer:
            pg3 =  SHPOLYS[ipar] 
            if pg3.contains(pg2):
              # {ipo2} is a closer container:
              ipar = ipo2
            else:
              assert pg2.contains(pg3), "contours are not properly nested"
    if ipar != None:
      PARENT[ipo1] = ipar
    
  # Now define the children (nearest contained polys) of each simple
  # poly. The children of {POLYS[ipo]} are {POLYS[kch]} for {kch} in
  # {CHILDREN[ipo]}.
  CHILDREN = [None]*npo 
  for ipo in range(npo): CHILDREN[ipo] = []
  
  for ipo in range(npo):
    ipar =  PARENT[ipo]
    if ipar != None:
      CHILDREN[ipar].append(ipo)

  def from_poly_forest(ROOTS, sgn_inf, ind):
    # Builds a shape from a collection of polyon trees whose root polys are pairwise un-nested.
    # Expects each root poly to be CCW (island) if {sgn_inf = +1}, CW (hole) if {sgn_inf = -1}.
    sys.stderr.write("%*s> forest roots = %s sgn_inf = %+d\n" % (2*ind,"",str(ROOTS),sgn_inf))
    ARGS = []
    for ipo in ROOTS:
      S1 = from_poly_tree(ipo, sgn_inf, ind+2)
      ARGS.append(S1)
    if sgn_inf == +1:
      # Trees are islands:
      S = rootray_shape.union(ARGS)
    else:
      # Trees are holes:
      S = rootray_shape.intersection(ARGS)
    sys.stderr.write("%*s< forest S = %s\n" % (2*ind,"",str(S)))
    return S
    # ....................................................................

  def from_poly_tree(ipo, sgn_inf, ind):
    # Builds a shape from the polygon tree whose root is {POLYS[ipo]}. 
    # Expects that poly to be CCW (island) if {sgn_inf = +1}, CW (hole) if {sgn_inf = -1}.
    PTS = POLYS[ipo]
    sys.stderr.write("%*s>tree ipo = %d sgn_inf = %+d\n" % (2*ind,"",ipo,sgn_inf))
    assert sgn_inf == hacks.poly_orientation(PTS), "inconsistent polygon orientation"
    R = rootray_shape.prism(PTS) # Root island or hole.
    ARGS = [ R ]
    sys.stderr.write("%*sR = %s\n" % (2*ind+2,"",str(R)))
    CHILDS = CHILDREN[ipo]
    if len(CHILDS) != 0:
      # Build shape from nested polygons:
      C = from_poly_forest(CHILDS, -sgn_inf, ind+2)
      sys.stderr.write("%*sC = %s\n" % (2*ind+2,"",str(C)))
      ARGS.append(C)
    if sgn_inf == +1:
      # Root is an island:
      S = rootray_shape.intersection(ARGS)
    else:
      # Root is a hole:
      S = rootray_shape.union(ARGS)
    sys.stderr.write("%*s<tree S = %s\n" % (2*ind,"",str(S)))
    return S
    # ----------------------------------------------------------------------

  # Collect the indices of the root contour :
  ROOTS = [ ipo for ipo in range(npo) if PARENT[ipo] == None ]


  # Convert to shape:
  S = from_poly_forest(ROOTS, sgn_inf, 0)

  return S
  # ----------------------------------------------------------------------
