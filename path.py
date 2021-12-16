# Tools for handling tool-paths.
# Last edited on 2021-11-09 13:59:25 by stolfi

import path_IMP; from path_IMP import Path_IMP
import move
import move_parms
import pyx

class Path(Path_IMP):
  # An object of the {Path} class represents a /tool-path/, or /path/
  # for short: the traveling of the nozzle along a polygonal path --
  # sometimes down, extruding material, sometimes raised, not extruding.
  # It is a sequence of zero or more oriented moves (traces or jumps),
  # such that the final point of each one is the initial point of the
  # next one.  The {Move} objects of all those moves must be distinct.
  #
  # Every {Path} object {ph} has an /initial point/ and a /final point/
  # the initial and final position of the nozzle as it fabricates the
  # path. These points are defined even if for a /trivial/ path with zero
  # moves, in which case they are the same point. If the path is not
  # trivial, they are the initial point of the first move and the last
  # point of the last move.
  #
  # The path may begin and/or end with a jump. It makes no sense for a
  # path to have two consecutive jumps, since those can be condensed
  # into a single jump with significant time saving. It also does not
  # make sense to have jumps of zero length. On the other hand, a trace
  # of zero length makes sense, provided that it is not followed or
  # preceded by another trace: it deposits a dot of material.
  #
  # The /execution time/ of path is the time needed to fabricate
  # all its moves, plus one /transition pernalty time/ 
  # for each internal transition from a trace to a jump or 
  # vice-versa.  This penalty is assumed to be part of the 
  # execution time of the jump.
  #
  # A {Path} object also has a mutable /name/ attribute, that is used
  # for debugging and documentation. It must be either a string or {None}.
  #
  # An /oriented path/ is a pair {(ph,dr)} where {ph} is a {path} object
  # and {dr} is 0 or 1 to indicate a direction of motion. If {dr} is zero
  # (the /native/ orientation), moves are assumed to be fabricated in the
  # order and orientations specified by the {Path} object. If {dr} is 1,
  # the moves are assumed to be fabricated in the reversed order, each in
  # the opposite direction. Note that the initial and final points of the
  # path are swapped in this case. A {Path} object {ph} by itself is
  # generally treated as the oriented path {(ph,0)}.
  pass
   
# ATTRIBUTES

def nelems(oph):
  # Returns the number of moves in the oriented path {oph}.
  return path_IMP.nelems(oph)

def elem(oph, imv):
  # Returns the oriented move with index {imv} in the oriented path {oph=(ph,dr)},
  # counting in the order implied by the orientation.  The index {imv}
  # must be in {0..nmv-1} where {nmv = nelems(oph)}.
  return path_IMP.elem(oph, imv)

def find_move(oph, omv):
  # Given an oriented path {oph} and an oriented or unoriented move {omv},
  # returns the index {imv} such that {elem(oph, imv)} is {omv} or its 
  # reverse.  If the move does not occur in {oph}, returns {None}.
  return path_IMP.find_move(oph, omv)

# GEOMETRY

def pini(oph):
  # Returns the initial point of the oriented path {oph}, taking the
  # orientation bit into account.
  return path_IMP.pini(oph)

def pfin(oph):
  # Returns the final point of the oriented path {oph}, taking the
  # orientation bit into account.
  return path_IMP.pfin(oph)

def endpoints(oph):
  # Returns a tuple {(pini(oph),pfin(oph))}.
  return path_IMP.endpoints(oph)

def points(oph):
  # Returns a list with all endpoints of all moves in the path {oph}.
  return path_IMP.points(oph)

def find_nearest_move_endpoint(oph, p):
  # Returns the index {imv} of the move endpoint in the path {oph} that
  # is nearest to the point {p} and the distance between the two.
  #
  # If the path is a single point, returns 0 and the distance from {p}
  # to {pini(oph)==pfin(oph)}. If the returned index {imv} is in
  # {0..nmv-1} where {nmv=nelems(oph)} (always the case if the path is
  # non-trivial and closed), the nearest point will be
  # {move.pini(elem(oph, imv)}.
  #
  # If {oph} is trivial or not closed, the index {imv} may be equal to
  # {nmv}, and then the nearest point will be {pfin(oph)}.
  return path_IMP.find_nearest_move_endpoint(oph, p)
   
def find_nearest_midline_point(oph, p):
  # Finds the point {q} on the midline of {oph} that is clsoest to the given point {p}.
  # 
  # Returns the closest point {q}, the index {imv} of a move of {oph} that contains it,
  # the total length {t} of the midline from {pini(oph)} to {q},
  # and the distance {d} from {q} to {p}. Beware that the computed length {t} may be
  # affected by roundoff errors.
  #
  # The meaning of {imv} is the same as in {find_nearest_move_endpoint}, except that
  # the valu emay be different, since the closest point may not be in the move that 
  # has the closest endpoint.
  return path_IMP.find_nearest_midline_point(oph, p)

def find_midline_point_by_length(oph, t):
  # Returns the index {imv} of the move of {oph} that contains the point {q} that is at distance {t}
  # from {pini(oph)} along the midline of {oph}. Also returns {q}.
  #
  # If {t} is zero, returns {imv=0} and {q=pini(oph)}. if {t} is
  # negative, or greater than or equal to the computed length of the
  # path, returns {None,None}. Beware that the computed length may be
  # affected by roundoff errors.
  return path_IMP.find_midline_point_by_length(oph,t)

def extract_section(oph, imv0,q0, imv1,q1):
  # Creates a path {sec} that is a copy of path {oph} between point
  # {q0} on the move {omv0} with index {imv0} and point {q1} on move
  # {omv1} with index {imv1}.
  #
  # The procedure assumes that all moves of {oph} are traces with
  # non-zero length and the same {MoveParms} record {mp}. Otherwise the
  # results may be inconsistent. The result {sec} may have up to two
  # moves with new {Mone} objects, also with parameters {mp}; but
  # otherwise its moves will share the {Move} objects in {oph}.
  #
  # If either point is {pini(oph)}, the corresponding index can be {-1},
  # which will be taken as 0. If either point is {pfin(oph)}, the
  # corresponding index can be {n} where {n = nelems(oph)}, and will be
  # taken as {n-1}.
  #
  # Apart from those special cases, {imv0} and {imv1} must be in {0..nmv-1}.
  #
  # If {imv0} is greater than {imv1}, the result will be {None}.
  #
  # Now assume {imv0 <= imv1}. Let {omv0} be {elem(oph,imv0)}, and let
  # {p0a,p0z} be its endpoints. Let {omv1,p1a,p1z} be the same for
  # {imv1}. As per above, we will assume that {p0a != p0z} and {p1a != p1z}.
  #
  # Consider the case {inv0==imv1}. Then {q0} and {q1} must be points of the segment 
  # {p0a--p0z}.  If {q0==q1}, the result {sec} will be the trivial path at 
  # that point.  If {q0} is ahead of {q1} onthe segment {p0a--p1z},
  # the result is {None}. if {q0==p0a} and {q1==p0z}, then {sec} will be 
  # a path that has {omv0} as its single move.  Othwewise {sec}  
  # will be a path with a new single move with endpoints {q0} and {q1}.
  #
  # Consider now the case {imv0 < imv1}. If {q0==p0a}, the returned path {sec} will start with
  # move {omv0}. If {q0==p0z}, the path {sec} will not have
  # any part of {omv0}. Otherwise the first move of {sec} will be a new
  # {Move} object with endpoints {q0} and at {p0z}.
  #
  # The situation for {q1,omv1,p1z,p1a} is symmetrical.
  #
  # If {oph} has a non-null name {pname}, any new {Move} objects created by the procedure
  # will be given names "{pname}.X{imv}" where {imv} is the index of the move of {oph} that was 
  # truncated to make this one. The new path itself will not have a name.
  return path_IMP.extract_section(oph, imv0,q0, imv1,q1)

def trim(oph, d0, d1):
  # Returns a path {otr} that is {oph} with a sections of length {d0}
  # chopped off from the beginning and a section of length {d1} chopped
  # off from the end.
  # 
  # The first and last moves of {otr} are usually new; any others are
  # shared with {oph}.  
  #
  # May return a trivial path (a single point). Returns {None} if the
  # computed length of the midline of {oph} is less than {d0 + d1}.
  # Beware that the computed length may be affected by roundoff errors.
  return path_IMP.trim(oph, d0, d1)

def bbox(OPHS):
  # The parameter {OPHS} must be a list or tuple of oriented paths.
  # Returns a bounding box of all moves of all those paths,
  # as a pair of points {(plo, phi)} that its lower left and upper right
  # corners.  If {OPHS} is empty, returns {None}
  #
  # The bounding box is the smallest axis-aligned rectangle that
  # contains all the endpoints of all moves and jumps in all those
  # paths, as well as the starting/ending point of any trivial paths. Note
  # that the nominal "sausages" of traces, as well as decoration such as
  # dots and arrowheads, may extend outside the box.
  return path_IMP.bbox(OPHS)

def mean_projections(oph, xdir, ydir):
  # Obtains the pair of mean coordinates of the endpoints of {oph}
  # along the directions of the unit vectors {xdir} and {ydir}.
  # If either direction is {None}, the coordinate is {None}.
  return path_IMP.mean_projections(oph, xdir, ydir)

# CREATION

# Unless said otherwise, the paths returned by these procedures
# have their link and contact lists set to empty.
 
def make_trivial(p):
  # Creates and returns an (unoriented) trivial {Path} object, with zero
  # moves, that begins and ends at the point {p}.
  return path_IMP.make_trivial(p)
   
def from_move(p, q, mp):
  # Creates and returns an (unoriented) {Path} object, with exactly
  # one move from {p} to {q}.  The width and dynamics will be as defined
  # by the {Move_Parms} object {mp}
  return path_IMP.from_move(p, q, mp)
  
def from_moves(OMVS):
  # The argument {OMVS} must be a list of one or more moves,
  # such that each move ends where the next one begins.
  # Creates and returns a path consisting of those moves.
  #
  # The same {Move} object must not appear more than once in this list.
  # That is, {OMVS[i]} must not be {OMVS[j]} or its reverse for any
  # distinct {i,j} However, the proceure does not check this condition,
  # and subtle errors may result later if it is not satisfied.
  # 
  return path_IMP.from_moves(OMVS)

def from_points(PTS, mp_trace, mp_jump):
  # The argument {PTS} must be a list whose elements are points or lists of points.
  #
  # If {n} consecutive elements of {PTS} are points, the procedure creates
  # a sub-path of {n-1} traces with parameters {mp_trace}, connecting those points
  #
  # An element of {PTS} may be itself a sub-list of points. In that
  # case, those points are connected by traces, as above, and that
  # sub-path is connected by jumps to the previous and following path
  # elements.
  #
  # The argument {mp_jump} is used only if the resulting path has jumps, as 
  # per above.  Otherwise it may be {None}.
  return path_IMP.from_points(PTS, mp_trace, mp_jump)
 
def concat(OPHS, use_links, mp_jump):
  # The argument {OPHS} must be a list of one or more oriented paths.
  # Returns a new {Path} object {ph} that is the concatenation of all
  # those paths. 
  #
  # If the final point of a path {oph1} in the list is not not the
  # initial point of the next path {oph2}, a /connector/ is inserted
  # between them. If {use_links} is true, the procedure looks for a path
  # {olk} in {get_links(oph2)} that starts at {pfin(oph1)}. If there is
  # such a path, the procedure inserts {olk} between {oph1} and {oph2}.
  # If there is no such path, or {use_links} is false, the connector is
  # a jump with parameters {mp_jump},
  #
  # The resulting path {ph} will share all the {move.Move} objects of
  # the given paths and any added links. Note that, depending on the
  # inputs, it may have multiple consecutive jumps; they are not
  # automatically condensed.
  #
  # The same {Move} object must not be used in more than one of those
  # paths. However, the procedure does not check this condition, and
  # subtle errors may result later if it is not satisfied.
  #
  # The input links of {ph}, as per {get_links}, will be copied from the
  # aplicable input links of the paths of {OPHS}. Note that {OPHS}
  # starts with {k} trivial paths, the links of all paths {OPHS[0..k]}
  # may be used. However, any link {olk} that starts on a points of {ph}
  # itself (namely, that bridges two paths of {OPHS}) will be ignored.
  # The links of {rev(ph}) are simlarly set from the links of the
  # reverses of paths in {OPHS}, in reverse order.
  #
  # The set {get_contacts(ph,isd)}, for {isd} in {0..1}, is the union of
  # the contact sets {get_contacts(OPHS[iph],isd)} of the given paths.
  # These sets must be pairwise disjoint, since the {Move} objects are.
  # However, the path {ph} may close some contacts that were
  # partially closed by two distinct input paths.
  #
  # This operation does not modify the given paths. It takes computing time
  # proportional to the number of paths plus the total number of moves
  # in those paths.
  #
  # The computation cost of this operation is normally proportional to
  # the total number of moves in the given paths, that is, to {n =
  # nelems(ph)}. The cost includes a factor proportional to the number
  # of input and output links of {ph}, but this number should normally
  # limited to to a small constant independent of {n}.
  return path_IMP.concat(OPHS, use_links, mp_jump)
 
def displace(oph, ang, v):
  # Returns a copy of the oriented path {oph} rotated by {ang} radians
  # around the origin and translated by the vector {v},
  # in that order.
  #
  # All moves will be new {move.Move} objects, even if {ang}and {v} are
  # zero. The parameter object {Move_Parms} of each move will be
  # preserved, however all the timing data will be recomputed (even if the 
  # length of moves is not supposed to change).
  #
  # The resulting path will have no links and no contacts. 
  return path_IMP.displace(oph, ang, v)
  
def split_at_jumps(oph):
  # Breaks the path {oph} into one or more trace-only paths at the jumps (if any).
  # Returns the list {OPHS} of the trace-only paths, and a list {OJMS} of the jumps.
  return path_IMP.split_at_jumps(oph)
  
def simplify(oph, tol, mp_jump):
  # Builds a new path that has the same endpoints as {oph}, but where
  # some sequences of two or more traces may be replaced by fewer
  # traces, maybe just one.
  #
  # The replacement only happens if all the traces in each replaced
  # sequence have the same {MoveParms} record, and the endpoints of
  # those original traces are at most {tol} millimeters away from the
  # new nozzle trajectory.
  #
  # The {Move} objects of the traces of the old path that are not merged
  # with other traces are reused in the new one. Each new trace will use
  # that common {MoveParms} record of the replaced traces, and will be named
  # "T{imv}" where {imv} is the index in {oph} of the last of those traces.
  #
  # As a side effect, the procedure also replaces any sequences of two
  # or more successive jumps in {oph}, with any {MoveParms} by a single
  # jump with parameters {mp_jump}. If {oph} contains no jumps,
  # {mp_jump} may be {None}.
  #
  # The name of the path {oph}, as per {get_name}, is copied to the new
  # path. The group index, as per {get_group}, is copied too. Any links
  # associated with {oph} (and its reversal) are copied to the new path
  # (and its reversal).
  #
  # The contact-path associations, as per {get_contacts} and
  # {contact.get_side_paths}, are NOT transferred to the new path, since
  # they depend on the original moves.
  return path_IMP.simplify(oph, tol, mp_jump)

# ORIENTATION
  
def rev(oph):
  # Returns the reversal of the oriented path {oph}.
  return path_IMP.rev(oph)

def spin(oph, dr):
  # Applies {rev} to the oriented path {oph}, {dr} times.
  # Thus the result is {oph} if {dr} is even, and {rev(oph)}
  # if {dr} is odd.
  #
  # If the {oph} argument is an unoriented {Path} object {ph}, the
  # result is {(ph,dr)}.
  return path_IMP.spin(oph, dr)

def unpack(oph):
  # Checks whether {oph} is an oriented path.  Returns the underlying 
  # {Path} object and the direction bit as two separate results.
  # If {oph} is already an (undirected) {Path} object, returns {oph,0}.
  return path_IMP.unpack(oph)
  
# TIMING

def fabtime(oph):
  # Returns the total time to fabricate the path {oph}, including all jumps
  # in it.  The orientation of {oph} is irrelevant.  The time
  # includes the penalty times for all internal jump/trace transitions,
  # as specified in the {Move_Parms} records of the jumps.
  return path_IMP.fabtime(oph)

def tini(oph, imv):
  # Returns the total nozzle travel time from the initial point of the oriented path {oph}
  # to the INITIAL point of the oriented move {elem(oph,imv)}, including all acceleration
  # and deceleration times. 
  #
  # This time includes all the trace/jump transition penalties that
  # occur in moves {0..imv-1}, which are counted as adding to the
  # execution time of the jumps. Thus, the result includes a penalty for
  # the transition from move {imv-1} to {imv} if the former is a jump and
  # the latter is a trace -- but not vice-versa.
  #
  # For convenience, if {imv == nelems(oph)}, returns {fabtime(oph)}.
  return path_IMP.tini(oph, imv)
  
def tfin(oph, imv):
  # Returns the total nozzle travel time from the initial point of the
  # oriented path {oph} to the FINAL point of the oriented move
  # {elem(oph,imv)}, including all acceleration and deceleration times.
  #
  # This time includes all the trace/jump transition penalties that
  # occur in moves {0..imv}, which are counted as adding to the execution
  # time of the jumps. Thus, the result includes a penalty for the
  # transition from move {imv} to {imv+1} if the former is a jump and the
  # latter is a trace -- not vice-versa.
  #
  # For convenience, if {imv == -1}, returns zero.
  return path_IMP.tfin(oph, imv)
 
# INCREMENTAL PATH BUILDING

def move_to(MVS, p, q, mp, tag):
  # If the point {p} is not {None},
  # appends to {MVS} a trace from {p} to {q} with parms {mp}.
  # The moves are named "M{tag}{kmv}" where {kmv} is sequentially incremented.
  # In any case, returns the point {q}.
  return path_IMP.move_to(MVS, p, q, mp, tag)

def finish(PHS, MVS, tag):
  # Makes the moves in {MVS} into a path, appends it to {PHS},
  # and resets {MVS} to empty
  # The path name is set to "P{tag}".
  return path_IMP.finish(PHS, MVS, tag)

# LINK PATHS  

# Each path {oph} can have a set of precomputed link paths that can be
# used to connect it to other paths. These links are associated to the
# endpoints of the underlying {Path} object of {oph} by {clear_links}
# and {add_link} below. The orientation of the path {oph} is used only
# to determine which endpoint the links are supposed to connect to.
#
# The links of {oph} are paths that end at {pini(oph)}.
#
# The links of {rev(oph)} are paths than end at {pfin(oph)},
# thus their reverses are paths that start at {pfin(oph)}.

def clear_links(oph):
  # Sets the list of links of the path {oph} (those that end at
  # {pini(oph)} to empty. Does not affect the links of {rev(oph)}
  path_IMP.clear_links(oph)
  
def add_link(oph, olk):
  # Appends the oriented link path {olk} to the list of 
  # links associated with {oph} object.
  #
  # The procedure fails if the link {olk} does not end at {pini(oph)}, or 
  # if its is a closed loop (with {pini(olk) == pfin(olk)}. 
  #
  # The link is NOT added if it starts at any endpoint of any move of
  # {oph}. Apart from that, the procedure does not check whether the
  # link intersects of gets too close to {oph}. It also does not check
  # for duplicate, overlapping, or redundant links.
  path_IMP.add_link(oph, olk)
   
def set_links(oph, OLKS):
  # Saves in the {Path} record of {oph} a copy of the list {OLKS}, as
  # the list of all links that connect to {oph}. They all must end at
  # {pini(oph)}. The links of {path.rev(oph}) are not affected. It is
  # equivalent to {clear_links(oph)} and then adding each link with
  # {add_link}.
  return path_IMP.set_links(oph, OLKS)
 
def get_links(oph):
  # Returns a list of all links stored in the {Path} record of {oph}
  # that end at {pini(oph)}.
  return path_IMP.get_links(oph)

def get_connecting_link(oph0, oph1):
  # Returns the link path that connects the end of oriented
  # path {oph0} to the start of the oriented path {oph1},
  # if such a link has been stored in the respective {Path}
  # objects; or {None} if there is no such link.
  #
  # There must be at most one link path satisfyng those conditions.
  # That path must have been associated with both {oph0} and 
  # {oph1} through {add_link}.
  return path_IMP.get_connecting_link(oph0, oph1)

def get_connector(oph0, oph1, use_links, mp_jump):
  # Returns a link path or jump that bridges the gap between paths {oph0} to {oph1},
  # if any. 
  # 
  # Namely, if {pfin(oph0) == pini(oph1)}, returns {None}. Othwerwise,
  # if {use_links} is true and there is a link path that connects the
  # two paths, returns that link path. Otherwise, the result is a {Move}
  # object is a jump across the gap,] with parameters {mp_jump}.
  return path_IMP.get_connector(oph0, oph1, use_links, mp_jump)

def get_all_connecting_links(oph0, oph1):
  # Returns a list of all link paths that connect either endpoint of
  # oriented path {oph0} to either endpoint of the oriented path {oph1}.
  # 
  # The list may be empty, and it will have at most one link,
  # obtained through {get_connecting_link}, for each pair of endpoints.
  return path_IMP.get_all_connecting_links(oph0, oph1)


def connection_time(oph0, oph1, use_links, mp_jump):
  # Estimates the extra fabrication time of {concat([oph0, oph1,], use_links, mp_jump)}
  # compared to the fabtimes of paths {oph0} and {oph1}.
  # 
  # If {pfin(oph0)} is the same as {pini(oph1)}, the result is zero.
  # Otherwise, {use_links} is true and {get_connecting_link(oph0, oph1)} is not {None},
  # returns the fabtime of that link.  Otherwise returns the time of a jump with parameters
  # {mp_jump}.  Includes the penalties of move-jump and jump-move transitions where applicable.
  return path_IMP.connection_time(oph0, oph1, use_links, mp_jump)
  
# PATH-CONTACT GRAPH

# Each {Path} object {ph} may have two sets of {Contact} objects
# assigned to it, accessded with {get_contacts} below.
# 
# Specifically, for each {ct} in {get_contacts(ph,isd)}, the move which is
# side {isd} of {ct} occurs in {ph}, in any orientation. In that case,
# {ph} and/or its reverse should be in {contact.get_side_paths(ct,isd)}.
#
# This information must be consitent with that provided by
# {contact.get_side_paths}. Specifically, for every contact {ct} in
# {get_contacts(ph,isd)}, the set {contact.get_side_paths(ct,isd)} must
# include {ph} and/or its reverse.

def clear_contacts(oph, isd):
  # Clears the lists of associated contacts that have their side {isd} used by the {Path} object {ph}
  # of the oriented path {oph}.  The orientation of {oph} is irrelevant.
  return path_IMP.clear_contacts(oph, isd)

def get_contacts(oph, isd):
  # Returns a copy of the set of contacts that have side {isd} used by the {Path} object {ph}
  # of the oriented path {oph}, as assigned by {add_contact}.  The orientation of {oph} is irrelevant.
  return path_IMP.get_contacts(oph, isd)

def add_contact(oph, isd, ct):
  # Appends {ct} to the set of contacts that have their side {isd} used by the {Path}
  # object {ph} underlying the oriented path {oph}, if it is not there yet.
  # The path must contain the move that is side {isd} of the contact. 
  # The orientation of {oph} is irrelevant.
  path_IMP.add_contact(oph, isd, ct)
  
# PATH GROUPS

# Each {Path} object can be assigned a {group} code, that can be a non-negative
# integer index or {None}. This field is used, for instance, to indicate
# how filling elements are to be joined into blocks for input to the
# toolpath contruction algorithms.
# 
# Whenever a {Path} object is created, its group code is ser to {None}.

def set_group(oph, igr):
  # Sets the group index of the {Path} object underlying the oriented path {oph} to {igr}.
  # Ignores the orientation of {oph}.  The new group {igr} must be {None} or a non-negative integer.
  path_IMP.set_group(oph, igr)
  
def get_group(oph):
  # Returns the last group index assiged to the {Path} object underlying {oph}
  # by the last call of {set_group}, or the default group (0) if there was no such call.
  # The index must be a non-negative integer.
  return path_IMP.get_group(oph)

def separate_by_group(OPHS):
  # Given a list {OPHS} of disjoint paths, outputs a list {GRPHS} of
  # lists of paths, consisting of the subsets of {OPHS} that have the same
  # group as specified by {path.get_group}; and the number {ngr} of non-empty
  # distinct groups found.
  #
  # More precisely, {GRPHS[igr]} will be the list of all paths of {OPHS}
  # with {path.get_group(oph) == igr}. Elements with group index {None}
  # will be ignored. The paths in {GRPHS[igr]} will be in the same order
  # as in {OPHS}. The list {GRPHS[igr]} will be empty if there is no
  # element with index {igr}.
  # 
  # The length of {GRPHS} will be the highest group index seen in {OPHS},
  # plus one. Note that the returned number {ngr} may be less than {len(GRPHS)}.
  #
  return path_IMP.separate_by_group(OPHS)

# NESTING OF CONTOUR PATHS

# A /contour path/ is a closed oriented path that is part of the slice's
# boundary. "Closed" means that its initial and final points coincide.
# Distinct contour paths must not cross or touch each other.
# 
# A contour path may be the outer boundary of an "island", in which case
# it must be oriented counterclockwise; or the boundary of a "hole", in
# which case it must be oriented clockwise. There may be holes in
# islands, and islands inside holes, nested to arbitrary depth.
# 
# The nesting of a set of contour paths can be precomputed by assigning
# to each contour path {ocr} in the set the list of other contour paths
# of the set that are contained in it, and the list of contour paths
# that contain it. These lists are built by {compute_contour_nesting}
# below. The contour paths must be properly oriented (CCW or CW) for
# this operation, but the lists are associated to the underlying {Path}
# object of each contour path {ocr}, ignoring its orientation. Therefore,
# these lists are shared with {path.rev(ocr)}, but not with any new
# {Path} object created by operations such as {path.shift_contour} or
# {path.displace}.

def compute_contour_nesting(OCRS):
  # Given a list {OCRS} of the contour paths of a slice, determines their nesting 
  # arrangement and saves in each {Path} object the lists used by {inner_contours}
  # and {outer_contours}.  The contour paths must be properly oriented: CCW for
  # island contours, CW for hole contours.
  path_IMP.compute_contour_nesting(OCRS)

def inner_contours(ocr):
  # Given a contour path {ocr}, returns the list (possibly empty) of the contour paths
  # that are immediately contained in it, as determined by the last call to
  # {compute_contour_nesting}. The orientaton of {ocr} is ignored,
  # and the returned paths are always properly oriented. If
  # {compute_contour_nesting} was never applied to {ocr}, returns {None}
  return path_IMP.inner_contours(ocr)

def outer_contour(ocr):
  # Given a contour path {ocr}, returns the innermost contour path
  # that encloses it, as determined by {compute_contour_nesting};
  # or {None} if there is no such contour. The
  # orientaton of {ocr} is ignored, and the returned path is always
  # properly oriented. If {compute_contour_nesting} was never applied to
  # {ocr}, returns {None}
  return path_IMP.outer_contour(ocr)

def contour_nesting(ocr0, ocr1):
  # Given two contour paths {ocr0,ocr1}, returns {+1} if {ocr0} contains
  # {ocr1}, {-1} if {ocr0} is contained in {ocr1} and 0 otherwise
  # (including if {ocr1} and {ocr2} are the same {Path} object). Uses
  # nesting information attached to the paths by
  # {compute_contour_nesting}. The orientatons of {ocr0} and {ocr1} are
  # ignored.
  return path_IMP.contour_nesting(ocr0, ocr1)
  
def shift_contour(ocr, kmv):
  # Given a contour path {ocr}, returns a contour path {ocrk} that has the same moves 
  # in the same sequence, except that it starts with the move {elem(ocr,kmv)}.
  # 
  # More precisely, one will have {elem(ocrk,imv) = elem(ocr,(imv+kmv)%n)} for
  # all {imv} in {0...nmv-1}, where {n=nelems(ocr)}. Note that the shift
  # amount {kmv} is taken modulo {n}. In particular, if {kmv=0} (or {kmv=n})
  # the result {ocrk} will be a copy of {ocr}. In any case, the result
  # will be a new {Path} object that shares the same {Move} objects with
  # {ocr}.
  return path_IMP.shift_contour(ocr, kmv)

# PLOTTING

def plot_to_files(fname, OPHS, CLRS, rwd, wd_axes, grid, deco):
  # The argument {OPHS} must be a list or tuple of oriented paths (or
  # plain {Path} objects). The {CLRS} argument must be a list or tuple
  # of {pyx.color.rgb} objects.
  #
  # Plots each path in {OPHS} using using some default sytle and the
  # correspondng color of {CLRS} for the trace sausages, over a
  # millimeter grid.  Trace axes, dots, and arrowheads will be plotted
  # iff {deco} is true.  Trace axes and jump lines will be drawn with line
  # width {wd_axes}. The fat sausage of a trace {mv} will have width
  # {rwd*move.width(mv)}.
  #
  # Plots a background grid iff {grid} is true.
  #
  # Then writes the plot to files "{fname}.{ext}" where {ext} is "eps",
  # "png", and "jpg".
  #
  # If {CLRS} has a single element, uses that color for all trace
  # sausages. If {CLRS} is {None}, uses some default color,  Otherwise {CLRS}
  # must have the same length as {OPHS}.
  path_IMP.plot_to_files(fname, OPHS, CLRS, rwd, wd_axes, grid, deco)

def plot_standard(c, OPHS, dp, layer, CLRS, rwd, wd_axes, axes, dots, arrows, matter, arrow_sz = 1):
  # The argument {OPHS} must be a list or tuple of oriented paths (or
  # plain {Path} objects). The {CLRS} argument must be a list or tuple
  # of {pyx.color.rgb} objects.
  #
  # Plots each path in {OPHS}, displaced by the vector {dp}, using the
  # correspondng color of {CLRS} for the trace sausages.
  #
  # If {axes} is true, draws the axes of traces, not just of jumps.
  #
  # If {dots} is true, prints dots at the ends of traces, not just of jumps.
  #
  # If {arrows} is true, draws arrowheads on traces, not just on jumps.
  #
  # If {matter} is {true}, shows the estimated area covered by the
  # material extruded during traces.
  #
  # The fat sausage of a trace will have width {rwd*wd} where {wd} is
  # its natural width. The axes of moves (traces or jumps) will be drawn
  # with width {wd_axes}: solid for traces, dashed for jumps. The dots
  # at the endpoints and arrows will be drawn with size proportional to
  # {wd_axes}. The color of these lines and decorations will black for
  # jumps, and darkened version of the sausage color for trace axes.
  #
  # If {CLRS} has a single element, uses that color for all trace
  # sausages. If {CLRS} is {None}, uses some default color,  Otherwise {CLRS}
  # must have the same length as {OPHS}. If {dp} is {None},
  # assumes {(0,0)} (no displacement).
  #
  # The plot is done in 4 passes or /layers/: (0) the material overflow
  # sausages of all traces, if requested; (1) the main sausage of all
  # traces; (2) the axes, dots, and arrows of all traces, as requested;
  # and (3) the axes, dots, and arrows of all jumps. 
  #
  # If {layer} is not {None}, it must be an integer in {0..3}, in which
  # case only that layer is plotted.
  path_IMP.plot_standard(c, OPHS, dp, layer, CLRS, rwd, wd_axes, axes, dots, arrows, matter, arrow_sz)

def plot_layer(c, oph, dp, jmp, clr, rwd, wd, dashed, wd_dots, sz_arrows):
  # Plots selected elements of the oriented path {oph} on the {pyx}
  # context {c}.
  #
  # Plots only the jumps if {jmp} is true, and only the traces if {jmp}
  # is false.
  #
  # If {rwd} and/or {wd} are positive, the axis of each
  # selected move {omv} is drawn as a fat line segment: a rectangle with
  # round caps at the endpoints. The line will be dashed if {dashed} is
  # true. The width of that line will be {rwd*width(omv) + wd}.
  # So, if the move is a jump, it will be just {wd}.
  #
  # If {wd_dots} is true, dots of that diameter will be plotted at both
  # ends of each selected move, even if the axis itself is not drawn.
  #
  # If {sz_arrows} is nonzero, an arrowhead of that size will be drawn
  # halfway along the axis of each selected move that is long enough,
  # even if the axis itself is not drawn.
  #
  # These items are drawn with color {clr}. If {clr} is {None}, the
  # procedure does nothing.
  #
  # The plot will be displaced by the vector {dp} (a 2-tuple
  # of floats). If {dp} is {None}, assumes {(0,0)} (no
  # displacement).
  #
  # If {None} is given as {wd_axis}, {wd_dots}, or {sz_arrow}, the value 0 is
  # assumed.
  path_IMP.plot_layer(c, oph, dp, jmp, clr, rwd, wd, dashed, wd_dots, sz_arrows)

def plot_single(c, oph, dp, split, clr):
  # Plots the path {oph}, displaced by {dp}, on the Pyx canvas {c} with a default style.
  # 
  # If the boolean {split} is false, plots the traces of the whole path with the same color {clr},
  # or  with a default color if {clr} is {None}.
  #
  # If the boolean {split} is true, then {clr} must be {None}; the procedure
  # breaks the path into sub-paths at its jumps, and plots each path 
  # with a different color, from a list of colors generated internally.
  return path_IMP.plot_single(c, oph, dp, split, clr)

# DEBUGGING
 
def validate(oph):
  # Runs some consistency checks on the oriented path {oph}.  
  # Aborts with assert failure if it detects any errors.
  path_IMP.validate(oph)

def compare(oph0, oph1, tol, die):
  # Checks whether the oriented paths {oph0} and {oph1} coincide -- that
  # is, whether they have the same number of moves of the same type
  # (traces or jumps) and all the corresponding move endpoints coincide,
  # apart from differences of {tol} in their coordinates.
  #
  # The two paths may have different orientation bits; the
  # orientation of each path is taken into account when obtaining the
  # order and orientation of its moves.
  #
  # Does NOT compare the {MoveParms} of the moves, the path names, the
  # fabtime-related attributes (like {tcov}), the group indices (as per
  # {get_group}) and the sets of link paths and contacts attached to
  # them with {get_links} and {get_contacts}. 
  #
  # If the paths coincide, {die} is false
  # returns {True} or {False} else fails if they are not;
  return path_IMP.compare(oph0, oph1, tol, die)
  
def check_links(OPHS,OLKS):
  # Given a list of raster paths {OPHS}, checks whether the link information 
  # provided by {get_links} is consistent. In particular, checks whether the
  # set of all links returned by that function is {OLKS}.
  return path_IMP.check_links(OPHS,OLKS)

def has_name(oph):
  # Returns {True} if and only if the name atribute of {oph}'s underlying {Path} object is not {None}.
  return path_IMP.has_name(oph)

def get_name(oph):
  # Given an oriented path {oph}, returns the name attribute of the
  # underlying {Path} object {ph}, if not {None}. If the name of {ph} is
  # {None}, returns instead "P?". In any case, that name is prefixed "~" if {oph} is
  # {rev(ph)}.
  return path_IMP.get_name(oph)

def set_name(oph, name, moves):
  # Given an oriented path {oph}, saves the string {name} as the name
  # attrbute of the underlying {Path} object {ph}. The orientation of
  # {oph} is ignored.
  # 
  # If {moves} is true, also sets the name of every move {mv} in {ph}
  # to "{name}.{KK}" where {KK} is the position of the move in the path.
  # The move name will be the same no matter whether the path uses {mv}
  # or {rev(mv)}.
  path_IMP.set_name(oph, name, moves)
  
def tag_names(OPHS, tag):
  # Prepends the string {tag} (if not {None}) to the names of all paths
  # in {OPHS}. The {Path} objects had better be all distinct.
  path_IMP.tag_names(OPHS, tag)

def show(wr, pref, oph, suff, moves, wna,wnm,wgr):
  # Writes to file {wr} a readable description of the path {oph}.
  #
  # The description has:
  #
  #   * the direction bit of the path, " " or "~"
  #   * the name of the underlying {Path} object, as returned by {get_name}
  #   * the number {nmv} of moves in the path
  #   * the group index as per {get_group}
  #   * the initial and final points
  #   * the names of atached incoming and outgoing link paths
  #   * the names of attached contacts
  #
  # If {moves} is true, the description also includes the list of the
  # names of the moves that comprise the path {oph}; namely, the
  # sequence {move.get_name(path.elem(oph,imv))} for {imv} in
  # {range(path.nelems(oph))}.
  # 
  # The incoming links are obtained with {get_links(oph)}, and the 
  # ougoing ones by reversing the paths returned by {get_links(rev(oph))}.
  # The name of each link is preceded by "<" if incoming, with ">" if
  # outgoing.
  #
  # The names of contacts are obtained with {get_contacts(oph,isd)} for
  # {isd=0} and {isd=1}. The name of each contact that is covered on
  # only one side {isd} by {oph} is followed by ":{isd}". The name of a
  # contact that is closed (covered on both sides) by {oph} is followed
  # by "*".
  # 
  # The name of the path {oph} is padded to {wna} columns and
  # left-aligned. The number of moves is padded to {wnm} columns. The
  # group index is padded to {wgr} columns. The whole description is
  # printed in a single line, is preceded by {pref} (if not {None}, and
  # followed by {suff} (if not {None}). It is NOT terminated by a
  # newline, unless provided by {suff}.
  path_IMP.show(wr, pref, oph, suff, moves, wna,wnm,wgr)

def show_list(wr, pref, OPHS, suff, moves, links):
  # Writes to file {wr} a readable description of the oriented paths in
  # the list {OPHS}. 
  # 
  # The output consists of a header and then the description of each
  # path, one per line, as produced with {show} with the given parameter
  # {moves}. Each line is prefixed by {pref} and the block index in
  # {BCS}, and followed by {suff} and a newline.
  #
  # If {links} is true, the above output is followed by a second
  # table that lists all the links in all those paths, in the same format.
  path_IMP.show_list(wr, pref, OPHS, suff, moves, links)
  
