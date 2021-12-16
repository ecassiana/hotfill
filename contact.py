# Tools and types for representing contacts between traces.
# Last edited on 2021-11-06 09:32:23 by stolfi

import contact_IMP; from contact_IMP import Contact_IMP
import move
import move_parms
import path
import pyx

class Contact(Contact_IMP):
  # An object of the {Contact} class represents a /contact/, a line
  # segment (possibly a single point) on the boundary of two distinct
  # traces, that must become a physical weld after both are extruded.
  # These traces are the /sides/ of the contact, arbitrarily indexed 0
  # and 1.
  #
  # A {Contact} object {ct} also has two /endpoints/. Their order is
  # arbitrary. Each endpoint is a list of 2 float coordinates, in
  # millimeters.
  #
  # A contact {ct} is /covered/ by a move {mv} if {mv} is one of its 
  # sides.  
  #
  # Acontact {ct} is /covered/ a path {ph} if at least one of its sides
  # is an element of {ph}, and is /closed/ by {ph} if both sides are.
  #
  # Acontact {ct} has a /cooling time limit/ attribute, positive number
  # that is the maximum time (in seconds) that should elapse betweem its
  # first covering and its closing by the tool-path. The default for new contacts is
  # {+inf}, meaning that there is no limit on the cooling time.
  #
  # A {Contact} object also has a mutable /name/ attribute, that is used
  # for debugging and documentation. It must be either a string or {None}.
  #
  pass

# All the procedures below that take an oriented path as parameter will
# also accept an unoriented {Path} object, which is taken in its native
# orientation.

def make(p0, p1, omv0, omv1):
  # Creates and returns a {Contact} object {ct} with endpoints {(p0,p1)}
  # whose sides are the oriented moves {omv0} and {omv1} -- which must be
  # traces,not jumps.  The cooling time limit is set to {+inf} (no limit).
  #
  # The contact's name will be "C({na0}:{na1})" where {na0} and {na1}
  # are the names of the two moves, if they are both defined; otherwise
  # it will be {None}.
  return contact_IMP.make(p0, p1, omv0, omv1)
 
def side_move(ct, isd):
  # The parameter {isd} must be 0 or 1. Returns the (unoriented) {Move}
  # object {mv} that is the trace on side {isd} of the contact.
  return contact_IMP.side_move(ct, isd)

def side_moves(ct):
  # Returns a tuple with the two (unoriented) {Move} objects
  # which are the sides of {ct}, namely {(side(ct,0), side(ct,1))}.
  return contact_IMP.side_moves(ct)

def which_side(omv, ct):
  # The parameter {omv} must be an unoriented move, and {ct} must be a
  # {Contact} object.  If {omv} (ignoring its orientation) is one of the sides of {ct}, returns the
  # index (0 or 1) of that side. Otherwise returns {None}.
  return contact_IMP.which_side(omv, ct)

def endpoints(ct):
  # Returns the endpoints of {ct}, as a pair of points, in no
  # particular order.
  return contact_IMP.endpoints(ct)

def pmid(ct):
  # Retuns the midpoint of the contact {ct}.
  return contact_IMP.pmid(ct)

def side_tcov(ct, isd):
  # The parameter {isd} must be 0 or 1. Returns the (precomputed) time
  # that the nozzle will take to move from the starting point of the move
  # {side(ct,isd)} to the point on the move axis that is closest to
  # {pmid(ct)}.
  return contact_IMP.side_tcov(ct, isd)
  
def bbox(CTS):
  # The parameter {CTS} must be a list of {Contact} objects. Returns a
  # bounding box for those contacts, as a pair of points {(plo, phi)} that
  # its lower left and upper right corners. If the list is empty,
  # returns {None}.
  #
  # The bounding box is the smallest axis-aligned rectangle that
  # contains endpoints (only) of all those contacts. Note that
  # decorations such as tics and arrowheads, as well as bits of the
  # contact line itself, may extend outside the box. Moreover, the
  # traces that are the sides of those contacts are NOT included in the
  # result.
  return contact_IMP.bbox(CTS)

# AUTOMATIC CREATION

def from_moves(omv0, omv1, mvdir, szmin, rszmin, tol):
  # The parameters {omv0} and {omv1} must be oriented moves. If they are
  # both traces with distinct {Move} objects, are approximately parallel
  # to the unit vector {mvdir}, and share a sufficient length of common
  # border, the procedure returns a contact between them spanning that
  # shared segment. Otherwise it returns {None}. The nominal traces
  # ("fat segments") of the two moves may overlap or be separated by up
  # to {tol}. The orientations of the moves are irrelevant.
  #
  # The procedure uses {move.shared_border(omv0,omv1,mvdir,tol)} to check the
  # geometry constraints and compute the shared border.  The
  # {mvdir} parameter may be {None} to mean the approximate 
  # mean direction of the two moves.
  #
  # The contact is created only if its length (in mm) is
  # at least {szmin}. If {rszmin} is positive, the contact
  # is created only if its length is at least {rszmin*L} where {L} is
  # the length of the shortest of the two traces. 
  return contact_IMP.from_moves(omv0, omv1, mvdir, szmin, rszmin, tol)

def from_move_lists(OMVS0, OMVS1, szmin, rszmin, tol, ydir):
  # Returns a list (possibly empty) {CTS} of {Contact}
  # objects between moves in the lists {OMVS0} and {OMVS1}.
  #
  # The parameters {OMVS0} and {OMVS1} must be lists of oriented moves.
  # The procedure calls {from_moves(mv0,mv1,None,szmin,rszmin)} with every
  # pair of distinct {Move} objects {mv0} from {OMVS0} and {mv1} from
  # {OMVS1}.
  #
  # Side 0 of each contact will be the move whose midpoint is lower in
  # the direction {ydir}. If {ydir} is {None}, it defaults to {(0,1)}
  # (the Y coordinate).
  #
  # If two moves {mv0} and {mv1} appear in both lists, only one contact
  # will be created for the pair.
  return contact_IMP.from_move_lists(OMVS0, OMVS1, szmin, rszmin, tol, ydir)

def from_paths(oph0, oph1, szmin, rszmin, tol, ydir):
  # Same as {from_move_lists}, with {OMVS0} and {OMVS1} being the moves
  # of the oriented paths {oph0} and {oph1}, respectively.
  #
  # The procedue also adds any contact created {ct} to the contact sets
  # of {oph0} and {oph1}, with {path.add_contact}; and adds {oph0} and
  # {oph1} to the side paths sets of {ct}, with {add_side_path} --
  # whenever the constact lies on those paths. Note that these
  # connections will not exclude pre-existing ones.
  return contact_IMP.from_paths(oph0, oph1, szmin, rszmin, tol, ydir)

def from_path_lists(OPHS0, OPHS1, szmin, rszmin, tol, ydir):
  # The parameters {OPHS0} and {OPHS1} must be lists of oriented paths.
  # Same as {from_move_lists}, with {OMVS0} and {OMVS1} being the moves in all paths
  # of {OPHS0} and {OPHS1}, respectively.  
  #
  # The proceure also adds all contact-path interactions in with
  # {path.add_contact} and {add_side_path}. Note that these connections
  # will not exclude pre-existing ones.
  return contact_IMP.from_path_lists(OPHS0, OPHS1, szmin, rszmin, tol, ydir)

def from_blocks(bc0, bc1, szmin, rszmin, tol, ydir):
  # Same as {from_path_lists}, with {OPHS0}  and {OPHS1} being
  # all choices of {Block} objects {bc0} and {bc1}, respectively.
  #
  # The proceure also adds all contact-path interactions in with
  # {path.add_contact} and {add_side_path}. Note that these connections
  # will not exclude pre-existing ones.
  return contact_IMP.from_blocks(bc0, bc1, szmin, rszmin, tol, ydir)

# COVERAGE BY PATHS

def path_tcov(oph, imv, ct, isd):
  # The parameter {ct} must be a {Contact} object, {oph} must be an
  # oriented path, and {imv} must be the index such that {side(ct,isd)} is
  # the move {path.elem(oph,imv} or its reverse. 
  #
  # The procedure returns the time {tcs} when the nozzle fabricates that
  # move, in the direction specified in the path, and and passes next to
  # the midpoint {m} of the contact {ct} -- specifically, on the point
  # that is closest to {m} on the axis of that trace.
  #
  # The time is counted from the beginning of the execution of the
  # oriented path {oph}, assuming that it is fabricated in the direction
  # specified.
  # 
  # Note that the resultis usually different for {oph} and for
  # {path.rev(oph)}. In fact, the two times are complementary relative
  # to {path.fabtime(oph)}.
  #
  # For convenience, the procedure returns {None} if {imv} is {None}.
  #
  # This procedure executes in {O(1)} time.
  return contact_IMP.path_tcov(oph, imv, ct, isd)

def path_ixcovs(oph, ct):
  # Returns a pair {ixs} of indices such that the move
  # {path.elem(oph,ixs[isd])} is {side(ct,isd)}, in either orientation. If the
  # move {side(ct,isd)} does not occur in {oph}, {ixs[isd]} is {None}.
  #
  # Note that the indices are usually different for {oph} and for
  # {path.rev(oph)}. In fact, they are complemented relative to {nmv-1},
  # where {nmv} is the number of moves in the path.
  #
  # This procedure executes in {O(nmv}) time.
  #
  # NOTE: this procedure should not be used if the path-contact links
  # {contact.get_side_paths(ct,isd)} have been set, snce those links
  # provide the indices of the moves in {O(1)} time.
  return contact_IMP.path_ixcovs(oph, ct)

def is_relevant(ct, BCS, ich):
  # Returns true if and only if at least one side of the contact {ct} is
  # a trace of at least one {Block} in {BCS}.
  #
  # If {ich} is {None}, considers all choices of every block. Otherwise
  # considers only choice {ich} of each block (if it exists).
  return contact_IMP.is_relevant(ct, BCS)

def path_tcovs(oph, ct):
  # Returns a pair {tcs} of floats, where {tcs[isd] = path_tcov(oph,imv[isd],ct,isd)}
  # where {imv[0],imv[1]} are the indices returned by {path_ixcovs(oph, ct)}.
  # However, if {imv[isd]} is {None}, {tcs[isd]} is {None} too.
  #
  # This procedure executes in {O(nmv)} time because it calls {path_ixcovs}.
  return contact_IMP.path_tcovs(oph, ct)

# PATH-CONTACT GRAPH

# Each side of a contact has two sets of paths that use the
# moves that are the sides of that contact.
# 
# More precisely, for any contact {ct}, the two sets are {PS[isd]=get_side_paths(ct,isd)} for
# {isd} in {0..1}.  Each element of {PS[isd]} is a triple
# {(ph, dr, imv)} such that {oph=(ph,dr)} is an oriented path
# and {imv} is an index such that {path.elem(oph,imv)} is the move that is
# side {isd} of {ct}, in either orientation. 
#
# If a set contains a triple {(ph,dr,imv)} it may or may not contain
# also the triple that refers to the reverse of that pat, namely {(ph,
# 1-dr, n-1-imv)}. Note that, since all the paths in {PS[isd]} share the
# same {Move} object, the final tool-path may include at most one path
# from each set.
#
# This information should be consistent with that provided by
# {path.get_contacts}. Namely, for every oriented path {oph} and every
# contact {ct} of interest, the path {oph} or its reverse must be in
# {contact.get_side_paths(ct,isd)} if and only if the contact {ct} is in
# {path.get_contacts(oph,isd)}

def clear_side_paths(ct, isd):
  # Resets the set of paths associated with side {isd} of {ct} to empty.
  contact_IMP.clear_side_paths(ct, isd)

def add_side_path(ct, isd, oph, imv):
  # Adds the oriented path {oph} to the set of paths associated with side {isd} of {ct}.  The path {oph}
  # must contain the {Move} object {mv} which is side {isd} (0 or 1) of {ct}, in any orientation.
  #
  # if {imv} is not {None}, it should be the index of the move {mv} in {oph}, and the operation has computing cost {O(1)}.  Otherwise
  # the procedure will compute the index with {path.find_move}, which may have cost proportional to the number of moves in the path.
  contact_IMP.add_side_path(ct, isd, oph, imv)

def get_side_paths(ct, isd):
  # Returns a copy of the set of triples {(ph,dr,imv)} that specify the
  # oriented paths of interest that cover side {isd} (0 or 1) {ct}, as
  # created by {add_side_path}.
  return contact_IMP.get_side_paths(ct, isd)

def check_side_paths(CTS,OPHS,ydir):
  # Given a list of contacts {CTS} and a list of raster paths {OPHS},
  # checks whether the contact information provided by {get_side_paths}
  # and {path.get_contacts} are consistent. 
  # 
  # In particular, checks whether the set of all paths returned by
  # {get_side_paths} is {OPHS}, and all contacts returned
  # by {path.get_contacts} are in {CTS}.
  #
  # If {ydir} is not {None}, its should be a unit 2-vector. The
  # procedure then also checks whether the two sides of each contact are
  # ordered by increasing coordinate of their midpoints along {ydir}.
  return contact_IMP.check_side_paths(CTS,OPHS,ydir)

# SIDE INDICES

# A {Contact} object also has two /side indices/, mutable attrbutes
# supposedly associated with its tow sides. They are recoverd and set by
# {get_side_index} and {set_side_index} below. Their value is either
# {None} or an integer. They are set to {None} when the contact is
# created.
#
# These side indices are completely independent from the path-contact
# pointers accessed through {get_side_paths}, {add_side_paths} and
# {clear_side_paths}.

def set_side_index(ct, isd, ix):
  # Sets the side index for side {isd} of {ct} to {ix},
  # which must be {None} or an integer.
  contact_IMP.set_side_index(ct, isd, ix)
  
def get_side_index(ct, isd):
  # Gets the side index for side {isd} of {ct}.
  return contact_IMP.get_side_index(ct, isd)

# COOLING TIME 

def tcool_limit(ct):
  # Returns the cooling time limit of contact {ct}.
  return contact_IMP.tcool_limit(ct)

def set_tcool_limit(ct, tclim):
  # Sets the cooling time limit of contact {ct} to {tclim}.
  contact_IMP.set_tcool_limit(ct, tclim)

def tcool(oph, ct):
  # The parameter {oph} must be an oriented path and {ct} must be a
  # {Contact} object, which is closed by {oph} (meaning that
  # the move on each side of {ct}, or its reverse, occurs in {oph}).  
  #
  # Returns cooling time of {ct} in that path -- namely the fabtime of the
  # part of {oph} between the the covering of the two contacts.
  # 
  # If the path {oph} does not close {ct}, returns {None}.
  return contact_IMP.tcool(oph, ct)

def coldest(oph, CTS, n):
  # Finds the {n} contacts in {CTS} that are closed by the path {oph}
  # and have the maximum cooling time.
  return contact_IMP.coldest(oph, CTS, n)

def est_rcool(oph, ct, mp_jump, quick):
  # The parameter {oph} must be an oriented path and {ct} must be a
  # {Contact} object.  Returns a lower bound for the cooling time ratio
  # of {ct} for any path that begins with {oph} and closes the contact. 
  # 
  # In any case, if {tcool_limit(ct)} is {+inf}, returns 0.
  # 
  # Otherwise, if {oph} covers both sides of {ct}, returns the cooling
  # time {tcool(oph,ct)} of {ct} for that path, divided by
  # {tcool_limit(ct)}.
  #
  # Othwerwise, if the path {oph} covers exactly one side of {ct}, side
  # {isd}, returns a lower bound for {max_rcool(P,ct)} over any path {P}
  # that begins with {oph} and eventually covers side {1-isd} by including
  # some path {oph1} in {get_side_paths(ct,1-isd)}.
  # 
  # The above estimate considers how much fabtime elapsed between the
  # first coverage of {ct} and the end of {oph}, plus the minimum
  # fabtime that would be neeed to go from the end of {oph} to the start
  # of {oph1}, plus the time needed to fabricate {oph1} up to the
  # midpoint of {ct}; the sum divided by {tcool_limit(ct)}. The gap
  # between {oph} and {op1} will assume the link paths associated with
  # {oph0,oph1} if available, otherwise it will assume a jump with
  # parameters {mp_jump}
  #
  # The above estimate uses the path {oph1} in {get_side_paths(ct,1-isd)}
  # that gives the minumum cooling time. The procedure fails if no paths
  # were associated to that side.
  #
  # If the path {oph} does not cover any side of {ct}, returns 0.
  #
  # In any case, if the result is greater than 1, any path that starts
  # with {oph} (using a choice of {bc}, if needed, to close {ct}) will
  # be invalid, because it will exceed the cooling time of that contact.
  # 
  # The above description applies if {quick} is false. If {quick} is
  # true, the procedure may prematurely return {+inf} as soon as it
  # determines that the result would be greater than 1.
  return contact_IMP.est_rcool(oph, ct, mp_jump, quick)
  
def est_max_rcool(oph, CTS, mp_jump, quick):
  # The parameter {oph} must be an oriented path, {CTS} must be a list
  # or tuple of {Contact} objects. Returns the maximum value of
  # {est_rcool(oph,ct,mp_jump,quick)} over every contact {ct} of {CTS}.
  #
  # If the list is empty, or no contact of {CTS} is partially or totally
  # covered bt {oph}, returns {0}.
  #
  # In any case, if the result is greater than 1, any path that starts
  # with {oph} and closes all contacts of {CTS} with the paths
  # associated to them, will be invalid, because it will exceed the
  # cooling time of at least one closed or partially covered contact.
  # 
  # The above applies if {quick} is false. If {quick} is true, the
  # procedure may prematurely return {+inf} if it determines that the
  # result would be greater than 1.
  return contact_IMP.est_max_rcool(oph, CTS, mp_jump, quick)

# PLOTTING

def plot_to_files(fname, CTS, clr, dashpat, ext, OPHS, CLRS, rwd, wd_axes, tics, arrows):
  # The arguments {CTS}, {OPHS}, and {CLRS} must be lists or tuples of
  # {Contact} objects, oriented paths (or plain {path.Path} objects),
  # and {pyx.color.rgb} objects, respectively.
  #
  # Plots each path in {OPHS} using using some default sytle and the
  # correspondng color of {CLRS} for the trace sausages, over a
  # millimeter grid. Trace axes and jump lines will be drawn with line
  # width {wd_axes}.  Each trace is drawn as a sausage of width {rwd*wd}
  # where {wd} is the nominal width.
  #
  # Then draws each contact from {CTS} using the color {clr}, dash
  # patterm {dashpat}, and extension amnount {ext}. Will draw tics at
  # the contacts' midpoints if {tics} is true, or arrows if {arrows} is
  # true. (See {plot_single} for details.)
  #
  # Then writes the plot to files "{fname}.{ext}" where {ext} is "eps",
  # "png", and "jpg".
  #
  # Typically, both sides of every contact in {CTS} should be traces of
  # paths in {OPHS}.
  #
  # If {CLRS} has a single element, uses that color for all trace
  # sausages. If {CLRS} is {None}, uses some default color, Otherwise
  # {CLRS} must have the same length as {OPHS}.
  contact_IMP.plot_to_files(fname, CTS, clr, dashpat, ext, OPHS, CLRS, rwd, wd_axes, tics, arrows)

def plot_single(c, ct, dp, clr, dashpat, ext, wd, sz_tic, arrow):
  # Plots the contact {ct} on the {pyx} context {c}, as a solid line
  # segment of width {wd} and color {clr}, with the nominal endpoints, with round caps. The plot
  # will be displaced by the vector {dp} (a 2-tuple of floats).
  # 
  # If {dashpat} is {None} the line will be solid, else {dashpat} should
  # be a dashing pattern spec as for {hacks.plot_line}. If {ext} is
  # nonzero, the contact endpoints will be displaced out by {ext}.
  #
  # If the float {sz_tic} is positive, draws a short tic perpendicular
  # to the contact line at its midpoint.  However, if the boolean {arrow} is true, 
  # plots a short arrowhead pointing from side 0 to side 1 instead of 
  # the tic.
  contact_IMP.plot_single(c, ct, dp, clr, dashpat, ext, wd, sz_tic, arrow)

# PRINTING AND DEBUGGING

def has_name(ct):
  # Returns {True} if and only {ct}'s name atribute is not {None}.
  return contact_IMP.has_name(ct)

def get_name(ct):
  # Given a {Contact} object {ct}, returns its name attribute, if not {None}. 
  # If that attribute is {None}, returns "C?" instead.
  return contact_IMP.get_name(ct)

def set_name(ct, name):
  # Saves the string {name} as the name attrbute of the {Contact} object {ct}.
  contact_IMP.set_name(ct, name)

def tag_names(CTS, tag):
  # Prepends the string {tag} (if not {None}) to the names of all
  # {Contact} objects in {CTS}. The objects had better be all distinct.
  contact_IMP.tag_names(CTS, tag)

def compare(ct0, ct1, tol, die):
  # Checks if the contacts {ct0} and {ct1} have the same endpoints, and
  # if the moves that are their sides are coincident traces, apart from
  # differences of at most {tol} in each coordinate.
  #
  # Does NOT compare the names and cooling time limits, or the {MoveParms} records
  # of the two moves. Also does not compare the paths associated by
  # {get_side_paths}.
  # 
  # If the contacts coincide, returns {True}, and {die} is irrelevant.
  # If they don't coincide, retuns {False} if {die} is false, and
  # bombs out if {die} is true.
  return contact_IMP.compare(ct0, ct1, tol, die)
  
def show(wr, pref, ct, suff, wna):
  # Writes {ct} on {wr} in a human-readable format.
  #
  # The output has the name and endpoints of the contact, the cooling
  # time limit, and the moves and relevant paths that are the sides of 
  # the contact.
  #
  # Namely, for {isd} in {0,1}, the output will show name of the {Move}
  # object that is side {isd} of the contact, and (between "@{" and "}")
  # the names of the paths attached to the contact that are supposed to
  # contain that move. The index {imv} of the {Move} in each path is
  # also shown as "[{imv}]" or "[~{imv}]" depending on whether it is
  # used in its native or reversed orientation.
  #
  # The name of the contact is padded to {wna} columns, left-aligned. .
  # The whole description uses one line, is preceded by {pref} (if not
  # {None}, and followed by {suff} (if not {None}). It is NOT terminated
  # by a newline, unless provided by {suff}.
  contact_IMP.show(wr, pref, ct, suff, wna)

def show_list(wr, pref, CTS, suff):
  # Writes to file {wr} a readable summary description of the contacts
  # in the list {CTS}. 
  #
  # The output consists of a header and then the description of each
  # contact, one per line, produced by {show}. Each line is prefixed by
  # {pref} and the index of the contact in {CTS}, and followed by {suff}
  # and a newline.
  contact_IMP.show_list(wr, pref, CTS, suff)

def show_times(wr, OPHS, CTS):
  # Prints to {wr}  the fabtimes (total, trace, jump) of all the paths in {OPHS}.
  # If {ophs} has a single path, prints the cooling times of the contacts in {CTS}
  # determined by that path.
  contact_IMP.show_times(wr, OPHS, CTS)

def write_times (wr, OPHS, CTS, execution_time, Nrast, Trast, Nlink, Tlink, Njump, Tjump):
  #   
  contact_IMP.write_times(wr, OPHS, CTS, execution_time, Nrast, Trast, Nlink, Tlink, Njump, Tjump)
