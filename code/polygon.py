# Tools and types for representing the outline of a slice of the object.
# Last edited on 2021-06-04 03:37:27 by jstolfi

import polygon_IMP; from polygon_IMP import Polygon_IMP

class Polygon(Polygon_IMP):
  # A {Polygon} object represents a simple polygon on the plane,
  # such as one component of the outline of a slice of the object. It is stored as a list
  # of points, the vertices. The list is implicitly closed, with an edge
  # from the last vertex to the first one.  The order of the vertices
  # is clockwise for internal contours (holes) and counterclockwise
  # for external ones (islands).
  pass

def from_points(PTS):
  # Creates a {Polygon} object from a list {PTS} of points.
  return polygon_IMP.from_points(PTS)

def get_points(pg):
  # Returns the list of vertices (points) of the {Polygon} object {pg}.
  return polygon_IMP.get_points(pg)

def make_path(pg, mp):
  # Creates a closed {Path} object from the point list of the {Polygon} object {pg},
  # using {path.from_points} with the {Move_Parms} object {mp} for the moves.
  # The path will include a closing move from the last point to the first one.
  return polygon_IMP.make_path(pg, mp)

def bbox(PGS):
  # Given a list of {Polygon} objects, returns the bounding box
  # of all their points.  If the list is empty, returns {None}.
  return polygon_IMP.bbox(PGS)

def from_point_lists(PTSS):
  # The argument {PTSS} should be a list of lists of points that
  # describes the outline of a slice of the object. Each element of
  # {PTSS} should be the list of vertices of a polygon that is the
  # boundary of island or hole in a slice of the object.
  #
  # The procedure converts each element of {PTSS} into a {Polygon}
  # object, using {from_points}, and returns the list of those polygons.
  return polygon_IMP.from_point_lists(PTSS)
