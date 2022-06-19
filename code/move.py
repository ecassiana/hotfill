# Types and tools to handle moves (traces or jumps).

import move_IMP; from move_IMP import Move_IMP
import move_parms
import pyx

class Move(Move_IMP):
  # An object of the {Move} class represents a /move/, a
  # straight-line movement of the nozzle: either a /trace/ if the nozzle
  # is down and extruding material, or a /jump/ if the nozzle is raised
  # and not extruding.
  #
  # The /axis/ of a {Move} object {mv} is the line segment traversed by
  # the center of the nozzle. The /endpoints/ of the move are those of
  # that segment. The order of the two points is not necessarily the
  # direction of motion of the nozzle (which may be unknown). Each
  # endpoint is a list of 2 float coordinates, in millimeters.
  #
  # A {Move} object {mv} has a read-only set of /move parameters/,
  # encoded as a {Move_Parms} object. See {move_parms.py}. These include
  # the move's /nominal width/ in mm, {width(mv)}, that is mostly used
  # when generating G-code or plotting. For a trace of length {L} (mm),
  # the volume of material deposited will be {V = L*T*width(mv)} (mm^3),
  # where {T} is the slice's thickness (mm). (This formula is valid in
  # the limit of large {L}, since it does not account for the roundish
  # "caps" of material deposited at the ends of the move.) The nominal
  # width is zero if and only if the move is a jump.
  #
  # Another read-only attribute of a {Move} is its /execution time/
  # {fabtime(mv)}, in seconds. It depends on the length of the move and
  # on the dynamic parameters specified by the {Move_Parms} object. The
  # execution time of a jump does NOT include the trace/jump transition
  # penalties; these are accounted for when computing the execution time
  # of multi-move tool-paths.
  #
  # A {Move} object also has a mutable /name/ attribute, that can be useful 
  # in debugging or documentation.  It must be either a string or {None}.
  #
  # An /oriented move/ is a pair {(mv,dr)} where {mv} is a {Move} object
  # and {dr} is 0 or 1 to indicate the direction of motion of the nozzle
  # along the axis. A {Move} object {mv} by itself is assumed to be the
  # oriented move {(mv,0)}.
  pass
  
def make(p0, p1, mp):
  # Creates and returns a {Move} object with endpoints {(p0,p1)}, with
  # the nominal width and timing parameters specified by the
  # {Move_Parms} object {mp}. The latter specifies whether the move is a
  # trace or a jump.  The {name} is set to {None}.
  return move_IMP.make(p0, p1, mp)
  
def parameters(omv):
  # Returns the {Move_Parms} object of the move {omv}.
  return move_IMP.parameters(omv)

def is_jump(omv):
  # Returns {True} if {omv} is an oriented or unoriented jump -- that is,
  # if its nominal width is zero.  Returns {False} if it is {None} or a trace.
  return move_IMP.is_jump(omv)

def is_trace(omv):
  # Returns {True} if {omv} is an oriented or unoriented trace -- that is,
  # if its nominal width is positive.  Returns {False} if it is {None} or a jump.
  return move_IMP.is_trace(omv)
  
def rev(omv):
  # Returns the reversal of the oriented move {omv}, namely the same
  # {Move} object with the orientation bit {dr} complemented.
  return move_IMP.rev(omv)
  
def spin(omv, dr):
  # Applies {rev} to the oriented move {omv}, {dr} times.
  # Thus the result is {omv} if {dr} is even, and {rev(omv)}
  # if {dr} is odd.
  #
  # If the {omv} argument is an unoriented {Move} object {mv}, the
  # result is {(mv,dr)}.
  return move_IMP.spin(omv, dr)

def unpack(omv):
  # If {omv} is an oriented move, returns the underlying 
  # {Move} object and the direction bit as two separate results.
  # If {omv} is already a move object, returns {omv} and {0}.
  # Also checks whether the object is indeed a {Move}.
  return move_IMP.unpack(omv)

# GEOMETRY

def width(omv):
  # Returns the nominal width of the oriented or unoriented move {omv}.
  return move_IMP.width(omv)

def pini(omv):
  # Retuns the initial endpoint of an oriented move {omv}, taking the
  # orientation bit into account.
  return move_IMP.pini(omv)
  
def pfin(omv):
  # Retuns the final endpoint of an oriented move {omv}, taking the
  # orientation bit into account.
  return move_IMP.pfin(omv)
  
def endpoints(omv):
  # Returns the endpoints of the axis of the oriented move {omv}, namely the pair
  # {(pini(omv), pfin(omv))}. Note that the order depends on the
  # orientation bit of {omv}.
  return move_IMP.endpoints(omv)

def length(omv):
  # Returns the Euclidean length of the axis of {omv}, that is, the
  # distance between {pini(omv)} and {pfin(omv)}. Note that it does not
  # consider the round caps at the ends of a trace.
  return move_IMP.length(omv)

def bbox(OMVS):
  # The parameter {OMVS} must be a list of oriented moves. Returns a
  # bounding box for those moves, as a pair of points {(plo, phi)} that
  # its lower left and upper right corners. If the list is empty,
  # returns {None}.
  #
  # The bounding box is the smallest axis-aligned rectangle that
  # contains endpoints (only) of all those moves. Note that decorations
  # such as dots and arrowheads, as well as the nominal "sausages" of
  # traces, may extend outside the box.
  return move_IMP.bbox(OMVS)

def displace(omv, ang, disp, mp):
  # Returns a copy of the oriented move {omv} rotated by {ang} radians
  # around the origin and translated by the vector {disp},
  # in that order.  
  #
  # The move will be a new {Move} object, even if {ang} and {disp} are
  # zero. The parameters are replaced by the {Move_Parms} object {mp},
  # and the execution time is recomputed.
  return move_IMP.displace(omv, ang, disp, mp)
  
def shared_border(omv0, omv1, mvdir, tol):
  # The procedure checks whether the moves {omv0,omv1} are roughly
  # parallel to the unit direction vector {mvdir} (a pair of floats),
  # their fat segments (nominal extents) touch each other, and their
  # projections along {mvdir} overlap by a significant amount. If that is
  # the case, the procedure returns the endpoints {p0,p1} of the line
  # segment that is the intersection of their boundaries, in increasing
  # order along {mvdir}. The orientation of the moves is ignored.
  #
  # If the fat segments of the two moves are separated or overlap
  # in the in the direction perpendicular to {mvdir} by more than {tol}
  # millimeters, returns {None,None}. The same is returned if 
  # either or both moves are jumps.
  # 
  # The procedure uses a small tolerance when deciding whether the moves
  # are parallel, their boundaries touch, and the {mvdir} projections overlap.
  #
  # The {mvdir} can be {None}, in which case the procedure will try to
  # guess a suitable {mvdir} from the mean direction of the two moves.
  # Iin this case, if the guess fails, the procedures returns
  # {None,None}.
  return move_IMP.shared_border(omv0, omv1, mvdir, tol)

def sort_by_midpoint(omv0, om1, ydir):
  # Returns the oriented moves {omv0} and {omv1}, swapped
  # if needed so that the midpoint of {omv0} is lower than that
  # of {omv1} in the direction {ydir}.
  #
  # If {ydir} is {None}, uses {(0,1)}, that is, sort s by the Y coordinate.
  return move_IMP.sort_by_midpoint(omv0, om1, ydir)

# TIMING
  
def fabtime(omv):
  # Returns the execution time of the oriented move {omv}.
  # Does NOT include the nozzle up/down or filament suck/spit 
  # at the start/end of jumps.
  return move_IMP.fabtime(omv)

def ud_penalty(omv):
  # If {omv} is a jump, returns the additional time penalty 
  # for raising/lowerng the nozzle and/or retracting/refeeding the filament if  
  # the jump were to be preceded or followed by a trace in a path.
  # If  {omv} is a trace, returns zero.
  return move_IMP.ud_penalty(omv)

def transition_penalty(omv0, omv1):
  # If either or both of {omv0} or {omv1} is {None}, returns 0.
  # Otherwise they must be oriented moves, and {pfin(omv0)} must be the
  # same as {pini(omv1)}. If one is a trace and the other is a jump,
  # returns the trace/jump transition time penalty that would apply were
  # consecutive elements of a tool-path. Otherwise returns 0.
  return move_IMP.transition_penalty(omv0, omv1)

def cover_time(omv, m):
  # Returns the time when the nozzle passes by the point {m} while
  # executing the oriented move {omv} in the specified direction,
  # counted from the beginning of the execution.
  #
  # Specificaly, returns the time that the nozzle takes to travel from
  # {pini(omv)} to the point on the axis of the move that is closest to
  # point {m}. The result is always between 0 and {fabtime(omv)}.
  #
  # This function is intended for traces.  If used for a jump,
  # it does not include any trace/jump transition penalty time.
  return move_IMP.cover_time(omv, m)

# PLOTTING

def plot_to_files(fname, OMVS, CLRS, rwd, wd_axes):
  # The argument {OMVS} must be a tuple or list of oriented moves (or
  # plain {Move} objects). The {CLRS} argument must be a list or tuple
  # of {pyx.color.rgb} objects.
  #
  # Plots each move in {OMVS} using some default sytle and the
  # correspondng color of {CLRS} for trace sausages, with a millimeter
  # grid in the backgrounds. The trace axes and jump lines will be drawn
  # with line width {wd_axes}.  Each trace is drawn as a sausage of width
  # {rwd*wd} where {wd} is its nominal width.
  #
  # Writes the plot to files "{fname}.{ext}" where {ext}
  # is "png", "eps", and "jpg".
  #
  # If {CLRS} has a single element, uses that color for all trace
  # sausages. If {CLRS} is {None}, uses some default color,  Otherwise {CLRS}
  # must have the same length as {OMVS}.
  return move_IMP.plot_to_files(fname, OMVS, CLRS, rwd, wd_axes)
 
def plot_standard(c, OMVS, dp, layer, CLRS, rwd, wd_axes, axes, dots, arrows, matter):
  # The argument {OMVS} must be a tuple or list of oriented moves (or
  # plain {Move} objects). The {CLRS} argument must be a list or tuple
  # of {pyx.color.rgb} objects.
  #
  # Plots each move {omv} in {OMVS}, displaced by the vector {dp}, using the
  # corresponding color of {CLRS} for the trace sausages.
  #
  # If {omv} is a trace, draws it as a "sausage" -- a thick rectangle
  # with round caps at the end, with width {rwd*width(mv)}. 
  #
  # If {axes} is true, also plots the axis of each trace, not just of jumps. 
  # 
  # If {dots} is true, also plots dots at the endpoints of each trace, not 
  # just of jumps
  #
  # If {arrows} is true, also plots an arrowhead midway along each trace,
  # showing its orientation.
  #
  # If {matter} is {true}, shows the estimated area covered by the
  # material extruded during traces, as a gray sausage slightly wider
  # than the nominal width, painted under the trace.
  #
  # The axes, dots, and arrows of traces are drawn with a darkened
  # version of the color used for the sausage.
  #
  # If {omv} is a jump, the main and overflow sausages is omitted, and
  # the axis always drawn drawn as a dashed black line with dots and
  # arrowhead.
  #
  # The width of the axis of a trace or jump will be {wd_axes}, which
  # should be positive. The sizes of dots and arrows will be fixed
  # proportions of {wd_axes}.
  #
  # If {CLRS} has a single element, uses it for all traces. If {CLRS} is
  # {None}, uses a default color. Otherwise {CLRS} must have the same
  # length as {OMVS}. If {dp} is {None}, assumes {(0,0)} (no
  # displacement).
  #
  # The plot is done in 4 passes or /layers/: (0) the overflow sausage
  # of the estimated extruded material, if {omv} is a trace and {matter}
  # is true, (1) the main sausage, if {omv} is a trace, (2) the axis,
  # dots, and arrowheads, if {omv} is a trace and these items have been
  # requested, and (3) the axis, dots, and arrowhead, if {omv} is a
  # jump. 
  #
  # If the parameter {layer} is not {None}, it must be an integer
  # in {0..3}, in which case only that layer is plotted.
  move_IMP.plot_standard(c, OMVS, dp, layer, CLRS, rwd, wd_axes, axes, dots, arrows, matter)

def plot_layer(c, omv, dp, clr, wd, dashed, wd_dots, sz_arrow):
  # Plots the move (trace or jump) {omv} on the {pyx} context {c}.
  # 
  # Specifically, If {wd} is a positive number, draws the move's axis as
  # a fat line segment: a rectangle of width {wd} with round caps
  # centered at the move's endpoints. The move's nominal width is
  # ignored.
  #
  # If {dashed} is true, the axis line will be dashed.
  #
  # If {wd_dots} is a positive number, plots round dots of that diameter at the 
  # endpoints
  #
  # If {sz_arrow} is a positive number, draws an arrowhead with that size
  # at the midpoint of the axis, to show the move's direction.
  # The arrow is drawn even if the axis is not drawn.
  #
  # All items will be drawn with the color {clr}. If {clr} is {None}, does nothing.
  #
  # The plot will be displaced by the vector {dp} (a 2-tuple of floats).
  # If {dp} is {None}, assumes {(0,0)} (no displacement).
  #
  # If {None} is given as {wd}, {wd_dots}, or {sz_arrow}, the value 0 is
  # assumed.
  move_IMP.plot_layer(c, omv, dp, clr, wd, dashed, wd_dots, sz_arrow)

# DEBUGGING AND TESTING

def has_name(omv):
  # Returns {True} if and only if the name attribute of {omv}'s
  # underlying {Move} object is not {None}.
  return move_IMP.has_name(omv)

def get_name(omv):
  # Given an oriented move {omv}, returns the name attribute of the
  # underlying {Move} object {mv}, if not {None}. If the name of {mv} is
  # {None}, returns instead "T?" or "J?" depending on whether it is a trace
  # or a jump. In any case, that name is prefixed "~" if {omv} is
  # {rev(mv)}.
  return move_IMP.get_name(omv)

def set_name(omv, name):
  # Given an oriented move {omv}, saves the string {name} as the name attrbute of the underlying
  # {Move} object {mv}. The orientation of {omv} is ignored.
  move_IMP.set_name(omv, name)
  
def tag_names(OMVS, tag):
  # Prepends the string {tag} to the names of all moves in {OMVS}.
  # The {Move} objects had better be all  distinct.
  move_IMP.tag_names(OMVS, tag)

def show(wr, pref, omv, suff, wna):
  # Writes to {wr} a readable descrption of the oriented move {omv}.
  #
  # The description has: 
  #
  #   * the orientation bit {dr} of {omv}, as " " or "~"
  #   * the name of the underlying {Move} object, as returned by {get_name}
  #   * the initial and final points
  #   * the nominal width (or "jmp" if it is a jump)
  #   * an optional observation 'vertiacal', 'horizontal', 'trivial', etc.
  #
  # The name is padded to {wna} columns and left-aligned. The whole
  # output fits in one line, is preceded by {pref} (if not {None}, and
  # followed by {suff} (if not {None}). It is NOT terminated by a
  # newline, unless provided by {suff}.
  move_IMP.show(wr, pref, omv, suff, wna)

def show_list(wr, pref, OMVS, suff):
  # Writes to file {wr} a readable summary description of the oriented moves in the
  # list {OMVS}.
  #
  # The output consists of a header and then the description of each
  # move, one per line, produced by {show}. Each line is prefixed by
  # {pref} and its index in {OMVS}, and followed by {suff} and a
  # newline.
  move_IMP.show_list(wr, pref, OMVS, suff)
