# Miscellaneous utility functions for HotFill/HotPath.

import color
import hacks_IMP
import pyx

# MATH HACKS:

def real_quadratic_roots(A, B, C):
  # Returns the roots {t0} and {t1} of the equation {A t^2 + B t + C =
  # 0}, with {t0 < t1}, if they are real and distinct. Otherwise returns
  # {None,None}. Requires {A != 0}. May return {-inf,+inf}.
  return hacks_IMP.real_quadratic_roots(A, B, C)

# GEOMETRIC HACKS

def is_point(p):
  # Returns {True} if {p} is a pair (2-tuple or 2-list) of floats,
  # {False} otherwise.
  return hacks_IMP.is_point(p)

def same_point(p, q, tol):
  # Returns {True} if the coordinates of {p} and {q} differ by at most {tol} 
  return hacks_IMP.same_point(p, q, tol)

def poly_cleanup(PTS, dmin):
  # Returns a copy of the point list {PTS} with points that are closer than {dmin}
  # replaced by a single point.  Result may be funny if there are too many close 
  # points.
  return hacks_IMP.poly_cleanup(PTS, dmin)

def poly_orientation(PTS):
  # The argument {PS} should be a list of of vertices of a simple polygon.
  # Returns {+1} if the vertices are in {CCW} order, {-1} if in {CW}
  # order, 0 if cannot figure it out.
  return hacks_IMP.poly_orientation(PTS)

def filet_center(p, q, ro, rt):
  # Computes the center {t} of the smaller circle of radius {rt} that
  # is tangent to the circles of radius {ro} with centers {p} and {q}.
  # Returns the solution that lies to right of the line from {p} to {q}.
  return hacks_IMP.filet_center(p, q, ro, rt)

def circle_x(y, yctr, rad):
  # positive abscissa of point on circle with center {(0,yctr)} and radius {rad}
  # with ordinate {y}.  Returns 0 if {|y-yctr| >= rad}.
  return hacks_IMP.circle_x(y, yctr, rad)

def cut_poly_with_band(CS, ydir, y0, y1):
  # Given the list {CS} of vertices of a simple polygon {D}, and two
  # floats {y0,y1}, returns a list {LS} of the pieces of polygon that
  # cross the strip of the plane between the lines perpendicular to
  # {ydir} whose coordinates in that direction are {y0,y1}.
  #
  # Specifically, {LS[jpo][kpt]} will be point {kpt} of the piece with index
  # {jpo}. Only considers pieces that completely cross the strip.
  # The pieces will be listed in no particular order.
  return hacks_IMP.cut_poly_with_band(CS, ydir, y0, y1)
  
def clip_seg_to_strip(p,q, ydir, y0,y1):
  # Clips the segment from {p} to {q} to the interior of the strip of the plane
  # perpendicular to {ydir} between coordinates {y0} and {y1} in that direction.
  # Requires {y0 < y1}.
  #
  # If the segment intersects the interior of the strip, returns four
  # values {ps,ploc,qs,qloc} where {ps,qs} are the endpoints of the
  # clipped segment, At most one of them will be on the line {y0}, or on
  # the line {y1}. The value {ploc} is 0 if {ps} is on the {y0} line, 1
  # if {ps} is on the {y1} line, and {None} otherwise (that is, {ps} is
  # the original {p}). The value {qloc} has the same meaning for {qs}
  # and {q}.
  # 
  # Returns {None,None,None,None} if the segment is entirely outside or on the boundary of the strip.
  return hacks_IMP.clip_seg_to_strip(p,q, ydir, y0,y1)
  
def find_point_on_poly(CS,idr,g):
  # Given the list {CS} or vertices of a polygonal line, Returns an
  # index {i_first} such that {CS[i_first]} is the first vertex that is
  # at distance greater than {g} away from the endpoint {idr} (0 = beg,
  # 1 = end) along the polygonal.
  #
  # Also returns the point {q_first} that is at eactly that distance
  # from the endpoint.
  #
  # Fails if {g} is greater than or equal to the computed length of the
  # polygonal (which may be affected by roundoff errors) or the polygonal
  # has fewer than 2 vertices.
  return hacks_IMP.find_point_on_poly(CS,idr,g)

def trim_poly(CS, g0, g1):
  # Given the list {CS} or vertices of a polygonal line,
  # trims a length {g0} from the beginning and a length {g1} from the end,
  # and returns the list of vertices of the resulting polygonal line.
  # The first and last vertices are usually new; all the others are 
  # copied from {CS}.  Fails if the length of the given line is 
  # less than {g0+g1}.
  return hacks_IMP.trim_poly(CS, g0, g1)
  
def fbomb(die):
  # Returns {False} if die is {False}, else fails. Handy in procedures that 
  # verify things.
  hacks_IMP.fbomb(die)
  
# PLOTTING HACKS

def make_canvas(B, dp, scale, frame, grid, nx, ny):
  # Creates a {Pyx} canvas {c} for an array of plots, each with bounding
  # box {B}. The array will have {nx} columns and {ny} rows. Optionally
  # draws a frame and/or a grid on each plot, as requested by the
  # booleans {frame} and {grid}.
  #
  # The procedure returns the new canvas {c} and the sizes {szx,szy} of
  # each array slot.
  #
  # The frame will be drawn around each slot of the array 
  # but slightly displaced inwards from its edges. The grid
  # will span the slot minus a slightly wider margin, so it will
  # almost but not quite touch the frame.
  #
  # If {dp} is not {None}, it must be a 2-vector (pair of floats).
  # Then all frames and grids will be displaced by {dp}.
  #
  # If {scale} is not {None}, the procedure also sets the {pyx} scale
  # factor to that value. If {scale} is {None}, it sets the scale so
  # that the total plot will have a reasonable size in pixels.
  #
  return hacks_IMP.make_canvas(B, dp, scale, frame, grid, nx, ny)

def round_box(B, mrg):
  # Returns the given {box} rounded outwards to edges with integer
  # coordinates, at least {mrg} away from the original ones.
  return hacks_IMP.round_box(B, mrg)
  
def plot_line(c, p, q, dp, clr, wd, dashpat):
  # Draws a line segment from {p} to {q}, displaced by the vector {dp},
  # on the {pyx} context {c}, with color {clr} and width {wd}.
  # 
  # If {dashpat} is not {None}, it should be a pair {(dlen,dgap)}; the
  # line is dashed with dashes of length {dlen} separated by gaps of
  # length {dgap}. These dimensions are absolute (not relative to {wd}).
  # 
  # Nothing is drawn if {wd} is {None} or non-positive.  
  # If {clr} is {None}, it defaults to black.
  hacks_IMP.plot_line(c, p, q, dp, clr, wd, dashpat)

def plot_box(c, B, dp, clr):
  # Paints the interior of the box {B} displaced by the vector {dp} on the {pyx} context {c},
  # with color {clr}. 
  #
  # If {clr} is {None}, uses black. If {dp} is {None}, it defaults to {(0,0)}.
  hacks_IMP.plot_box(c, B, dp, clr)
  
def plot_text(c, tx, p, dp, fsize, texclr):
  # Writes the string {txt} with font size {fsize} (in points) and color {texclr} at point {p}
  # shifted by {dp}.
  # 
  # The {texclr} must be a string naming one of the predefined colors in
  # the {color.sty} LaTeX package. If {None}, the default color (black,
  # normally) is used. If {dp} is {None}, it defaults to {(0,0)}.
  hacks_IMP.plot_text(c, tx, p, dp, fsize, texclr)

def plot_frame(c, B, dp, clr, wd, dashpat, inset):
  # Plots on the {pyx} context {c} a frame around the box {B}, displaced
  # by the vector {dp}, with color {clr}, line width {wd}.
  #
  # If {inset} is not {None}, the frame will be displaced inwards by
  # {inset} (or outwards by {-inset}, if negative). The lines will be
  # dashed with the dash pattern {dashpat} if it isnot {None}.
  #
  # If {dp} is {None}, it defaults to {(0,0)}. If {clr} is None,
  # defaults to black. If {wd} is {None} or non-positive, does nothing.
  hacks_IMP.plot_frame(c, B, dp, clr, wd, dashpat, inset)

def plot_grid(c, B, dp, clr, wd, dashpat, inset, xstep, ystep):
  # Plots on the {pyx} context {c} a grid spanning the box {B}. The grid
  # lines will be at coordinates that are integer multiples of {xstep}
  # and {ystep}.
  #
  # If {inset} is not {None}, the grid edger will be displaced inwards
  # by {inset} (or outwards by {-inset}, if negative) relative to the
  # box {B}. If {dp} is not {None}, the whole grid will be displaced by
  # {dp}. The lines will be dashed with the dash pattern {dashpat} if it
  # is not {None}.
  #
  # If {dp} is {None}, it defaults to {(0,0)}. If {clr} is {None},
  # defaults to a light gray. If {wd} is {None} or non-positive, does
  # nothing.
  hacks_IMP.plot_grid(c, B, dp, clr, wd, dashpat, inset, xstep, ystep)
  
def write_plot(c, name, scale_n = 4):
  # Writes the {pyx} canvas {c} to an EPS file "{name}.eps", 
  # a PNG file named "{name}.png", and a JPEG file "{name}.jpg".
  hacks_IMP.write_plot(c, name, scale_n)

def adjust_dash_pattern(rdist, rlen, rgap):
  # Adjusts a dash pattern so that a dashed line will begin and end
  # with full dashes.
  #
  # The parameter {rdist} should be the distance to be traced. The
  # procedure assumes that the ideal dash pattern should have dashes of
  # length {rlen} separated bt gaps of lenth {rgap}. (These measurements
  # may be absolute, or relative to some common unit, e. g. the line
  # width.)
  #
  # If the procedure succeeds, it returns the adjusted dash pattern
  # {(a*rlen, a*rgap)} where {a} is a suitable adjustment factor, near
  # 1. The factor {a} will be such that {rdist} will be exactly
  # {a*((n+1)*rlen + n*rgap)}.
  #
  # The procedure fails if the distance {rdist} is too small (compared
  # to {rlen} and {rgap}) to use a dashed line. In that case it returns
  # {None}.
  #
  # ??? Generalize to allow the first and last dashes be some given multiple of {rlen} ??? 
  return hacks_IMP.adjust_dash_pattern(rdist, rlen, rgap)

def choose_grid_steps(sz):
  # Returns the major and minor grid steps appropriate for a picture
  # that measures {sz} millimeters in the largest dimension.
  # They will be powers of 10 times some simple factor like 1, 2, 5.
  # The minor grid step will be a sub-multiple of the major one.
  return hacks_IMP.choose_grid_steps(sz)

def trace_colors(nc, Y):
  # Returns an array {C} of {nc} colors that can be used to paint
  # different paths, etc. Each element is a {pyx.color.rgb} triple.
  # 
  # The colors will have mean brightness {Y}. If {Y} is {None}, a
  # default value will be used.
  return hacks_IMP.trace_colors(nc, Y)

def link_colors(nc, Y):
  # Returns an array {C} of {nc} colors that can be used to paint different links.
  # Each element is a {pyx.color.rgb} triple.
  # 
  # The colors will have mean brightness {Y}. If {Y} is {None}, a
  # default value will be used.
  return hacks_IMP.link_colors(nc, Y)

def matter_color(Y):
  # Returns the standard {pyx.color.rgb} value for the matter traces
  # to be used by procedures like {move.plot_standard} that do not
  # have that as a parameter.  
  #
  # If {Y} is not {None}, itmust be a number between 0 and 1, and the color will have that 
  # brightness.
  return hacks_IMP.matter_color(Y)

# IMAGE HACKS

def convert_eps_to_png(name):
  # Reads the file "{name}.eps" and converts it to a PNG image file "{name}.png".
  hacks_IMP.convert_eps_to_png(name)

def convert_eps_to_jpg(name):
  # Reads the file "{name}.eps" and converts it to a JPEG image file "{name}.jpg".
  hacks_IMP.convert_eps_to_jpg(name)

