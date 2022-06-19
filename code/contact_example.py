# Procedures for creating examples of contacts between paths.

import contact_example_IMP
import move
import move_parms
import path
import pyx

# CONTACTS

# SAMPLE CONTACTS FOR THE UNIT TEST PROGRAM

# These test data assume that the width of the filling trace is 1.0 mm, since
# then the nominal sausages of the sides actually touch each other.

def misc_A(mp_trace):
  # Returns a single contact {ct} between two new traces.
  #
  # The moves {side(ct,i)} are
  #
  #   i pini  pfin
  #   - ----- -----
  #   0 (1,1) (4,1)
  #   1 (2,2) (5,2)
  #
  # The endpoints of the contact are {(2,1.5)} and {(4,1.5)}.
  return contact_example_IMP.misc_A(mp_trace)

def misc_B(mp_trace, mp_jump):
  # Returns a list {CTS} of eight contacts between the ten traces of the
  # five oriented paths created by {path_example.misc_E} with
  # aruments {mp_trace, mp_jump}. Also returns the list {OPHS} of those
  # oriented paths and the list {TRS} of those ten traces.
  # 
  # The sides {side(CTS[k],i)} of the contacts, and their midpoints, are
  #
  #   k side0  side1  pmid       obs
  #   - ------ ------ ---------- -----
  #   0 TRS[0] TRS[1] (2.5, 3.5)
  #   1 TRS[1] TRS[2] (4.0, 3.5) single point
  #   2 TRS[0] TRS[2] (3.5, 3.0) single point
  #   3 TRS[2] TRS[3] (%%%, 2.0) 
  #   4 TRS[1] TRS[4] (5.0, 4.5) single point
  #   5 TRS[7] TRS[9] (3.5, 7.5)
  #   6 TRS[8] TRS[9] (5.5, 7.5)
  #   7 TRS[7] TRS[5] (3.5, 6.5)
  #  
  # where "%%%" means a bit more than 5.0.
  #
  return contact_example_IMP.misc_B(mp_trace, mp_jump)

def misc_C(mp_trace):
  # Creates two serpentine paths side by side, one with horizontal rasters, one with 
  # vertical rasters, and returns a list {CTS} of the contacts between them.
  # Also returns a list {PHS} of the two paths.
  return contact_example_IMP.misc_C(mp_trace)

def misc_F(alt, mp_trace, mp_jump):
  # Returns two results, a {Contact} object {ct} and a {Path} object {ph}.
  #
  # The contact {ct} lies between two adjacent traces of {ph},
  # that is expected to have the largest cooling time mong such 
  # contacts.
  #
  # The path {ph} is generated by {path.example.misc_F} with the 
  # {alt} parameter.
  #
  return contact_example_IMP.misc_F(alt, mp_trace, mp_jump)
 
def misc_K(mp_trace, mp_jump):
  # Returns a list {CTS} of three contacts between three paths, and 
  # the list {OPHS} of those paths.
  # Special for {contact_TST.py}.
  return contact_example_IMP.misc_K(mp_trace, mp_jump)

# TEST PROBLEMS FOR HOTFILL/HOTPATH

def road_of_blocks(nmv, nmg, org, ix0,ix1, iy0, mp_fill):
  # Returns two results: a list {BCS} of zero or more {Block} objects,
  # and a list {CTS} of zero or more {Contact} objects.
  #
  # The list {BCS} will have a bunch of blocks stacked vertically.
  # each block will have 2 or 4 choices which are snake paths,
  # comprising a total of {nmv} horizontal raster lines.
  #
  # If {nmg} is 1, each raster line will be a separate block, with two
  # possible choices, corresponding to the two orientations.
  #
  # If {nmg} is greater than 1, only a certain number {nms} of the
  # raster lines, in the middle of the road, will be blocks by
  # themseves, as above. The other {nmv-nms} raster lines will be fused
  # together into serpentine path blocks of {nmg} rasters each. Each of
  # these blocks will have 4 choices: the two possible orders for the
  # rasters (bottom to top, top to bottom) and the two possible
  # orientations for the bottom raster. The number {nms} will be at
  # least 1, and will be such that {nmv-nms} is divisible by {2*nmg}. In
  # particular, if {nmv} is less than {4*nmg+1}, each raster will be a
  # block by itself. All choices of each block will share the same
  # {Move} objects.
  #
  # All rasters will have width {wdf = width(mp_fill)}, and will be
  # spaced vertically by that amount. The X coordinates of the raster
  # endpoints will be {org[0] + ix0*wdf} and {org[0] + ix1*wdf}. The Y
  # coordinate of the bottom raster will be {org[1] + iy0*wdf}.
  #
  # The {Move_Parms} objects {mp_jump} and {mp_fill} will be used for
  # jumps and raster traces, respectively.
  #
  # There will be a contact between every pair of consecutive blocks in
  # each road.
  return contact_example_IMP.road_of_blocks(ix0,ix1, iy0, nmv, nmg, nbt, org, mp_fill)

def two_roads_and_islands(nmv, nmg, nis, mp_cont, mp_fill, mp_jump):
  # Returns two results: a list {BCS} of zero or more {Block} objects,
  # and a list {CTS} of zero or more {Contact} objects.
  #
  # The list {BCS} will have some "raster" blocks and some "island"
  # blocks.
  #
  # The raster blocks will form two parallel vertical roads, each
  # created by {road_of_blocks} with the given parameters
  # {nmv,nmg,mp_fill}. There will be {nis} islands arranged in a
  # vertical group at the left of the left road, between the two roads,
  # and at the right of the right road, near the middle of the Y span of
  # the roads.
  #
  # An island block will have a single choice that is a small circular
  # filled contour. It starts and ends at nearly the same point, so
  # there is no point in adding the reversal. There will be {nis}
  # islands arranged in a vertical group at the left of the left road,
  # between the roads, and at the right of the right road.
  # If the island has any filling, there will be a jump before 
  # each filling element, internal to the island path.
  #
  # The {Move_Parms} objects {mp_jump}, {mp_cont}, and {mp_fill} will be
  # used for the internal island jumps, island contours, and fillings, respectively.
  #
  # There will be a contact between every pair of consecutive blocks in
  # each road.
  return contact_example_IMP.two_roads_and_islands(nmv, nmg, nis, mp_cont, mp_fill, mp_jump)

# FROM BLOCKS

def raster_raster_contact(oph0, oph1):
  # Creates and returns a new {Contact} object {ct}
  # between the oriented paths {oph0} and {oph1}.
  #
  # Assumes that {oph0} and {oph1} are "snake" paths made of horizontal
  # raster lines connected by jumps or links, such that the topmost
  # trace {mv0} of {oph0} is adjacent to the bottommost trace {mv1} of
  # {oph1}, and they have overlapping {X} ranges. The contact will be
  # created with sides 0 and 1 set to {mv0} and {mv1}, respectively.
  # 
  # The contact will be undirected. Also adds {oph[i]}
  # to the set {contact.get_side_paths(ct,i)} and adds {ct} to the sets
  # {path.get_contacts(oph[i],i)}, for {i=0} and {i=1}; where
  # {oph=(oph0,oph1)};
  #
  return contact_example_IMP.raster_raster_contact(oph0, oph1)

