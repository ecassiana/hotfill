# Tools for raster fillings.

import raster_IMP
import path

# A /raster fill element/ is a path consisting of a single move which is
# a trace. A /raster filling/ is a collection of raster fill elements of
# the same width {wd} that are all parallel to some direction {xdir},
# and whose positions along the perpendicular direction {ydir} are all
# approximately integer multiples of {wd} plus some arbitrary shift, the
# /scanline phase/.

def get_spacing_and_phase(OPHS, xdir, ydir):
  # Given a list {OPHS} of filler elements, in arbtrary order, returns
  # estimtes of the average scanline spacing {ystep} (which should be
  # approximately the common trace width) and scanline phase {yphase}. 
  #
  # The phase should be the approximate {ydir} projection of the
  # lowest scanline in direction {ydir}.
  # 
  # Fails if the list is empty, or if the elements are not a proper raster
  # fill.
  return raster_IMP.get_spacing_and_phase(OPHS, xdir, ydir)

def point_scanline_index(p, xdir, ydir, ystep, yphase):
  # Returns the index of the scanline that contains the point {p}.
  # Assumes that the scanline  projections in the {ydir} direction are
  # approximately {yphase + isc*ystep}. 
  return raster_IMP.point_scanline_index(p, xdir, ydir, ystep, yphase)

def path_scanline_index(oph, xdir, ydir, ystep, yphase):
  # Returns the index of the scanline that contains some approximate midpoint of the path {oph}.
  # Assumes that the scanline  projections in the {ydir} direction are
  # approximately {yphase + isc*ystep}. 
  return raster_IMP.path_scanline_index(oph, xdir, ydir, ystep, yphase)

def sort_by_scanline(OPHS, xdir, ydir, ystep, yphase):
  # Given a list {OPHS} of filler elements, in arbtrary order, returns a
  # copy of the list with the same rasters sorted by increasing scanline
  # order.  Assumes that the rasters are parallel to the unit vector {xdir},
  # and that the projecttion of scanline with index {isc} on the perpendicular 
  # direction {ydir} is approximately {yphase + isc*ystep}.
  return raster_IMP.sort_by_scanline(OPHS, xdir, ydir, ystep, yphase)

def separate_by_scanline(OPHS, xdir, ydir, ystep, yphase):
  # Given a list {OPHS} of filler elements, in scanline order, returns a
  # list {SCS} of sub-lists, wehere each sub-list {SCS[isc]} contains
  # the indics into {OPHS} of the fillers on the scan-line with index
  # {isc}, also in scan-line order. 
  #
  # Assumes that the rasters are parallel to the unit vector {xdir};
  # that the projection of scanline with index {isc} on the
  # perpendicular direction {ydir} is approximately {yphase +
  # isc*ystep}; and that the scanline indices are all non-negative.
  return raster_IMP.separate_by_scanline(OPHS, xdir, ydir, ystep, yphase)

def make_raster_raster_link(oph0, oph1, xdir, dmax, dxmid, mp_link):
  # Creates a link path betwen the initial points of the raster fill
  # elements {oph0} and {oph1}, if appropriate. The link path will
  # consist of one or two traces with parameters {mp_link}. The
  # prcoedure returns the link path, or {None} if it is not created.
  #
  # The parameter {xdir} must be a unit vector (2-tuple of floats). The
  # link is created only if the two traces have the same general
  # orientation in the direction {xdir}, and their initial points are at
  # most {dmax} apart. Otherwise the procedure does nothing.
  #
  # If the link is created, it is attached the link to the lists of
  # links of the two paths with {path.add_link}.
  #
  # If {dxmid} is {None}, the link path will consist of a single trace
  # with parameters {mp_link}.
  #
  # If {dxmid} is not {None}, it must be a displacement (in mm). The
  # link path will then consist of two traces with parameters {mp_link},
  # and the middle point will be displaced by {dxmid} in the direction
  # opposite to that of the two traces.
  #
  # The names of the links and paths will be left undefined.
 return raster_IMP.make_raster_raster_link(oph0, oph1, xdir, dmax, dxmid, mp_link)
  
def make_raster_raster_contact(oph0,oph1, xdir, szmin, rszmin, tol):
  # The paths {oph0,oph1} should have a single move each. Attempts to
  # create a contact {ct} with {contact.from_moves} with their moves
  # {omv0,omv1} and the given parameters {xdir,szmin,rszmin,tol}.
  #
  # If it succeeds, attaches the contact {ct} to the two paths with
  # {path.add_contact} and {contact.add_side_path}, and returns {ct}.
  # Otherwise does nothing and returns {None}.
  return raster_IMP.make_raster_raster_contact(oph0,oph1, xdir, szmin, rszmin, tol)
  
def create_all_raster_raster_links_and_contacts(OPHS, xdir, ydir, dmax, dxmid, mp_link):
  # Given a list {OPHS} of raster fill elements, in arbitrary order,
  # creates contacts between adjacent rasters in consecutive 
  # scanlines, and link paths between their endpoints, whenever appropriate.
  #
  # Returns a list {LKS} of the {Path} objects of the created link
  # paths, and a list {CTS} of the created contacts.
  # Also attaches those items to the paths in {OPHS} with
  # {path.add_link}, {path.add_contact}, and {contact.add_side_path}.
  #
  # Specifically, the procedure calls {make_raster_raster_contact}, for
  # any two rasters {oph0,oph1} in consecutive scanlines; and
  # {make_raster_raster_link}, with the same {xdir,dmax,dxmid,mp_link}
  # parameters, for their four oriented versions.
  #
  # The name of eac contact will be set to
  # "C[{isc0}.{irs0}:{isc1}.{irs1}]" where {isc0,isc1} are the indices
  # of the scanlines and {irs0,irs1} are the indices of the rasters
  # within in each scanline.
  #
  # The name of each link path will be set to
  # "L[{isc0}:{irs0}{dr0}{isc1}:{irs1}{dr1}]" where {dr0,dr1} are "a" or
  # "z" to indicate whether the link connects to the {pini} (0) or
  # {pfin} (1) of the respective raster element. Its traces will be
  # named "{name}.{kmv}" where {name} is the above string and {kmv} is
  # the trace index in the path (0 or 1).
  # 
  # The ".{irs0}" and ".{irs1}" parts of these names are omitted if 
  # there is only one raster on that scan-line.
  #
  return raster_IMP.create_all_raster_raster_links_and_contacts(OPHS, xdir, ydir, dmax, dxmid, mp_link)

def endpoints_from_shape(S, B, xdir,ydir, ystep,yphase):
  # Computes the endpoints of traces that make upa solid raster fill,
  # from a description {S,B} of the region {D} to be filled.
  #
  # The argument {S} must be a two-dimensional {RootRay_Shape} object. It is implicitly 
  # intersected with the box {B} to yield the region {D}.
  #
  # The result is a list {PTRS} such that {PTRS[isc][jrs][kep]} is endpoint
  # {kep} (0=left, 1=right) of raster {jrs} on scan-line {isc} from bottom.
  # These endpoints will be on the boundary of {D}. Note that the traces
  # will extend outside of {D}.
  #
  # Rasters will be generated on every scan-line that crosses the box {B}.
  #
  # The scan lines will be parallel to the direction {xdir} and will be
  # spaced in the direction {ydir} by {ystep}. The axis of one of the
  # scan-lines passes through the point {yphase*ydir}. All scan-lines
  # that cross the box {B} will be considered.
  #
  # The returned rasters will be oriented the direction {xdir}. They will be
  # sorted by projction on {ydir} and ties brokem by projection on
  # {xdir}
  # 
  return raster_IMP.endpoints_from_shape(S, B, xdir,ydir, ystep,yphase)

def endpoints_from_polys(PTCS, xdir,ydir, ystep,yphase):
  # Computes the endpoints of the traces that make up a solid raster fill,
  # of a polygonal region {D}, given the vertex lists {PTCS} of {D}.
  # 
  # This procedure is equivalent to {endpoints_from_shape}, except that
  # the region {D} is a polygon defined by its vertices, instead of a
  # {RootRay_Shape} object and a box {B}.
  #
  # The parameter {PTCS} must be a list of lists of points. Each element
  # {PTS} of {PTCS} must be the vertices of a simple polygon that is a
  # connected component of the boundary of {D}. The last point of {PTS}
  # must be the same as the first point. The vertices must be in CCW
  # order if the polygon is the border of an "island", including the
  # outermost border; and CW if it is the border of a "hole". 
  #
  # Rasters will be generated on every scan-line that crosses the region {D}.
  #
  return raster_IMP.endpoints_from_polys(PTCS, xdir,ydir, ystep,yphase)

def link_points_from_raster_endpoints(PTCS, PTRS, ydir):
  # Given the polygons {PTCS} that define the perimeter of the area to be filled, 
  # and the endpoints {PTRS} of a collection of raster filling elements,
  # Returns the points that comprise the link paths between those raster elements.
  #
  # The parameter {PTCS} must be a list of lists of points, which are
  # the vertices of a polygonal region {D} as in {endpoints_from_polys}.
  #
  # The parameter {PTRS} must decribe the endpoints of a set of raster
  # fill elements, as would be produced by {endpoints_from_polys}.
  # Namely, {PTRS[isc][jrs][s]} should be the endpoint {s} (0 = left, 1 =
  # right) of raster {jrs} on scan-line {isc}. These raster endpoints must
  # lie on the border of {D}, apart from rounding error. Scanlines
  # should be perpendicular to {ydir}, sorted by increasing projection
  # in that direction.
  #
  # The result is a list {PTLS} such that {PTLS[isc][jrs][kep]} is point {kep} 
  # on the link path of index {jlk} between scan-lines {isc} and {isc+1}.
  # These ponts too lie on the boundary of the filling region {D}
  # and are sorted in inceasing Y coordinate.
  return raster_IMP.link_points_from_raster_endpoints(PTCS, PTRS, ydir)

def from_endpoints(PTRS, mp):
  # Creates a list {OPHS} of paths, each a single raster line, from the endpoints in {PTRS}.
  #
  # The parameter {PTRS} should be a list such that {PTRS[isc][jrs][kep]} is endpoint
  # {kep} (0=left, 1=right) of raster {jrs} on scan-line {isc} from bottom
  # 
  # The result {OPHS} of paths, each a single raster line with the given endpoinst
  # and {MoveParms} parameter record {mp}.
  #
  # Also calls {path.clear_links} for each raster element, in both senses.
  return raster_IMP.from_endpoints(PTRS, mp)

def analyze_fabtime(OPHS, ydir, ytol, minlen):
  # Analyzes the fabtimes of the paths in the list {OPHS}, classifying
  # each move of each path as raster, link, or jump. Returns the total
  # fabtimes {Trast,Tlink,Tjump} for each class.
  #
  # A trace of a path is considered a raster if the coordinates of its
  # endpoints on the direction {ydir} differ by at most {ytol}.
  # Otherwise the trace is assumed to be a link. In particular, if
  # {ytol} is zero, the move is considered a raster only if the
  # projections are equal.
  #
  # The trace-to-jump transition penalty {move.ud_penalty(ojp)} of a
  # jump {ojp} in a path {oph} is counted in {Tjump} only if {ojp} is
  # not the first move of {oph}, and the previous move was a trace.
  # Symmetrically, the jump-to-trace penalty is counted only if {ojp} is
  # not the last move of {oph}, and the next move is a trace.
  # 
  # The procedure also flags any move {omv} in those paths that is
  # classified as raster by the above criterion but has length (not
  # fabtime) less than {minlen}.
  #
  # The procedure also examines the links associated to each path {oph},
  # namely {path.get_links(oph)} and {path.get_links(rev(oph))}. It
  # flags any move in those links that is a jump or that would be
  # considered a raster by the above criterion. The latter would be
  # confused with a raster if {oph} is incorporated in a longer path.
  # However, the fabtimes of these moves are NOT included in {Trast} or
  # {Tlink}.
  return raster_IMP.analyze_fabtime(OPHS, ydir, ytol, minlen)
