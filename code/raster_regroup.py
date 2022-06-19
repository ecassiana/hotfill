# Procedures to group ratser filling elements in various ways.

import raster_regroup_IMP

# The procedures in this module assume that the parameter {OPHS} is a
# list of oriented paths that are elements of a solid raster filling, each
# possibly assigned to a /group/ whose index is accessed through
# {path.get_group}. The procedures split or redefine those groups
# into new groups according to various criteria, and attach the new
# group indices to the paths of {OPHS} with {path.set_group}.
# 
# Each path in {OPHS} should consist of a single raster trace. The paths
# must be pairwise disjoint and all parallel to some unit direction
# vector {xdir}. They must be sorted in /scanline order/, which is
# defined by increasing projection on an axis {ydir} perpendicular to
# {xdir}, ignoring roundoff errors; with ties broken by increasing
# projection along {xdir}. The orientations of the paths and traces are
# ignored by these procedures.
#
# Some of these procedures use the contact-path information provided bt
# {path.get_contacts} and {contact.get_side_paths}, which must be
# consistent with each other and with the list {OPHS}. Namely, for every path {oph} in {OPHS},
# and every contact {ct} in {path.get_contacts(oph,isd)},
# {contact.get_side_path(ct,isd)} is a subset of {OPHS} that has {oph}
# and/or its reverse.  Moreover, for this module, every contact {ct}
# must have side 0 below and side 1 above.
#
# Every procedure that is supposed to change the current group
# assignments, such as {split_at_forks} or {split_by_size}, will
# renumber the new groups with consecutive integers starting at 0, even
# if the partition of {OPHS} into groups remains the same. They return
# the number of new groups created.
#

def merge_all(OPHS):
  # Ignores the current group assignments of the paths in {OPHS}, even
  # those that have group index {None}, and reassigns all to group 0.
  return raster_regroup_IMP.merge_all(OPHS)

def by_contours(OPHS, OCRS):
  # Ignores the current group assignment and instead separates the filler
  # elements from {OPHS} into new groups, one for each contour (closed
  # oriented path) in the list {OCRS}. Each filler will be assigned to
  # the group with index {icr} such that {ocr=OCRS[icr]} is the
  # innermost contour that contains it.
  #
  # It is assumed that all fillers are in the interior of the slice.
  # Therefore, there will not be any filler assigned to a group index
  # {icr} if {icr >= len(OCRS)}, if {OCRS[icr]} is an inner (hole)
  # contour, or if contour {OCRS[icr]} is an outer contour but there is
  # no filler inside it except inside nested contours.
  return raster_regroup_IMP.by_contours(OPHS, OCRS)
  
def split_at_forks(OPHS): 
  # Splits the current groups of fillers from {OPHS} into maximal
  # sub-groups connected by single contacts.
  #
  # Specifically, each of the original groups is split into maximal
  # sub-groups such that, in every sub-group whose fillers are
  # {GPH[0..m-1]}, and for every {iph} in {1..m-1}, there is a contact
  # between {GPH[iph-1]} and {GPH[iph]}, which is the only contact on the
  # upper side of {GPH[iph-1]} and the only contact on the lower side of
  # {GPH[iph]}.
  #
  # Paths of {OPHS} which have group index {None} are not reassigned.
  #
  # This procedure uses the contact-path information provided by
  # {path.get_contacts} and {contact.get_side_paths}, which must be
  # consistent as spcified above.
  return raster_regroup_IMP.split_at_forks(OPHS)
  
def split_by_size(OPHS, ydir, max_lines):
  # Splits each current group, if needed, into sug-groups with at most 
  # {max_lines} whole scan-lines each.  
  #
  # Paths of {OPHS} which have group index {None} are not reassigned.
  #
  # In particular, if {max_lines} is 1, each new group will be contained
  # in a single scan-line. If that operation is applied after
  # {split_at_forks} each group will be a single filler element.
  raster_regroup_IMP.split_by_size(OPHS, ydir, max_lines)
