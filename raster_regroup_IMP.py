# Implementation of module {raster_regroup}
# Last edited on 2021-11-09 02:27:31 by stolfi

import block
import path
import move
import move_parms
import contact
import path_example
import raster
import job_parms

import txt_read

import hacks
import rn

from shapely.geometry import Point, Polygon
import sys
from math import sqrt, sin, cos, log, exp, floor, ceil, inf, nan, pi

def merge_all(OPHS):
  for oph in OPHS: path.set_group(oph, 0)
  return 0 if len(OPHS) == 0 else 1
  # ----------------------------------------------------------------------

def by_contours(OPHS, OCRS):
    
  ncr = len(OCRS)
  POLS = [None]*ncr # Contours converted to {shapely} polygons.
  
  # Convert all contours to {shapely} polygons {POLS[0..ncr-1]}:
  for icr in range(ncr):
    ocr = OCRS[icr]
    ocr_ph, ocr_dr = path.unpack(ocr)
    assert path.pini(ocr) == path.pfin(ocr)
    PTS = [ move.pini(path.elem(ocr,kmv)) for kmv in range(path.nelems(ocr)) ]
    POLS[icr] = Polygon(PTS)
    
  # Assign all fillers to group {ncr} to mean "not assigned"
  for oph in OPHS: path.set_group(oph, ncr)
  
  # Assign each filler from {OPHS} to a group whose index is the 
  # innermost contour that contains it:
  for icr in range(ncr):
    ocr = OCRS[icr]
    ocr_ph, ocr_dr = path.unpack(ocr)
    pg = POLS[icr]
    for oph in OPHS:
      cini = pg.contains(Point(path.pini(oph)))
      cfin = pg.contains(Point(path.pfin(oph)))
      assert cini == cfin
      if cini:
        # Found a contour {ocr} that contains {oph}.
        # Check whether it is the innermost one:
        icr1 = path.get_group(oph)
        if icr1 >= ncr:
          # Not assigned yet:
          path.set_group(oph, icr)
        else:
          pg1 = POLS[icr1]
          if pg1.contains(pg):
            # Current group contains {icr}
            path.set_group(oph, icr)
          else:
            # Current group does not contain {icr}, must be contained in it:
            assert pg.contains(pg1)
  return
  # ----------------------------------------------------------------------

def split_at_forks(OPHS):

  # Separate fillers by group, preserving their order, ignoring {None} group:
  GRPHS_old, ngr_old = path.separate_by_group(OPHS)
  
  ngr_new = 0  # Number of new groups created so far.
  for OPHS_gr in GRPHS_old:
    if len(OPHS_gr) > 0:
      ngr_new = split_rasters_at_forks(OPHS_gr, ngr_new)
  return
  # ----------------------------------------------------------------------
      
def split_by_size(OPHS, ydir, max_lines):
  debug = False
  assert type(max_lines) is int and max_lines >= 1
  
  xdir = ( +ydir[1], -ydir[0] )
  
  # Separate fillers by group, with each group still sorted by scanline order, ignoring {None} group:
  GRPHS_old, ngr_old = path.separate_by_group(OPHS)
  
  # Scan orig groups, splitting as requested:
  ngr_new = 0
  nph_new = 0
  for OPHS_gr in GRPHS_old:
    if len(OPHS_gr) > 0:
      # Split group into scanlines, preserving order:
      ystep, yphase = raster.get_spacing_and_phase(OPHS_gr, xdir, ydir)
      ngr_new = split_group_by_size(OPHS_gr, ngr_new, max_lines, xdir, ydir, ystep, yphase)
  return ngr_new
  # ----------------------------------------------------------------------
      
# INTERNAL PROCEDURES

def split_rasters_at_forks(OPHS, ngr):
  # Used by {split_at_forks}. Splits the rasters in {OPHS} into groups
  # at forks. Assumes that {ngr} groups have been created, and renumbers
  # the new groups starting at {ngr}. Returns the updated group count
  # {ngr}.
  #
  # This procedure uses the contact-path information provided by
  # {path.get_conatcts} and {contact.get_paths}.  See {split_at_forks}
  # in {path.py} for details.
  #
  # Each subgroup will be a maximal subset of {OPHS}, in
  # consecutive scan-lines, such that every two consecutive elements
  # {oph0,oph1} have a contact {ct} between them, and that
  # contact is the only contact of {oph0} on its top side, and the only
  # contact of {oph1} on its bottom side.
  # 

  OPHS_top = set()  # Top elements of still-incomplete new groups.

  # Each element of the set {OPHS_top} is the highest element (in the
  # {ydir} direction), among the elements processed so far, in some new
  # group that has not been completed yet.
  # 
  # Every element in that new group has exactly one contact on its upper
  # side (but not necessarily to an element of the input list {OPHS}).  
  # Every element in that new group, except the first one, has
  # exactly one contact on its lower side, to the element below it in the
  # same new group.
  
  ngr_new = ngr

  # Process the filler elements in of {OPHS} in order, assumed to
  # be lowest to highest in the {ydir} direction:
  for oph in OPHS:
    oph_added = try_adding_to_open_group(oph, OPHS_top, )
    if oph_added == None:
      # Could not add {oph} to any existing new group.  
      # Create another new group just for it:
      path.set_group(oph, ngr_new)
      ngr_new += 1
      OPHS_top.add(oph)

    # Check if the new sub-group is complete:
    CTS_hi = path.get_contacts(oph, 0) # Set of contacts on upper edge of {oph}.
    if len(CTS_hi) != 1:
      # The new sub-group is complete:
      OPHS_top.remove(oph)

  # At this point new groups that are still open must have a single upper contact 
  # with elements in other old groups. Ignore them.

  return ngr_new
  # ----------------------------------------------------------------------

def try_adding_to_open_group(oph, OPHS_top):
  # Used by {split_rasters_at_forks}. 
  # Tries to add the path {oph} to one of the open groups whose
  # top elements are in {OPHS_top}. If it succeeds, replaces that top
  # element by {oph} and returns it. Otherwise, returns {None}.
  
  debug = False
  name_oph = path.get_name(oph)

  CTS_lo = path.get_contacts(oph, 1) # Set of contacts on lower edge of {oph}.

  nct_lo = len(CTS_lo)
  if nct_lo != 1:
    # Element {oph} has no low-side contacts, or more than one such contact  
    # Cannot add to any existing group, no matter what those contacts are:
    if debug: sys.stderr.write("  raster %s has %d down-contacts, break\n" % (name_oph,nct_lo))
    return None
  else:
    # Element {oph} has only one low-side contact.  See if it matches 
    # any of the high-side contacts of the still open new groups:
    [ct_lo] = CTS_lo # Puts the lone element in {ct_lo}.
    name_ct_lo = contact.get_name(ct_lo)
    if debug: sys.stderr.write("  raster %s down contact is %s ..." % (name_oph,name_ct_lo))
    for oph_top in OPHS_top:
      # Consider adding {oph} to the open group of {oph_top}:
      CTS_top_hi = path.get_contacts(oph_top, 0) # Contacts above path {oph_top}
      nct_top_hi = len(CTS_top_hi)
      assert nct_top_hi == 1
      [ct_top_hi] = CTS_top_hi # Puts the lone element in {ct_top_hi}.
      name_ct_top_hi = contact.get_name(ct_top_hi)
      if debug: sys.stderr.write(" %s" % name_ct_top_hi)
      if ct_lo == ct_top_hi:
        # Contacts match. Add to new group:
        if debug: sys.stderr.write(" matched!\n")
        path.set_group(oph, path.get_group(oph_top))
        OPHS_top.remove(oph_top)
        OPHS_top.add(oph)
        return oph
    
    # Not continuation of any open group. Assume that the contact {ct_lo}
    # is with a path of some other old group:
    if debug: sys.stderr.write(" no match.\n")
    return None
  # ......................................................................

def split_group_by_size(OPHS, ngr, max_lines, xdir, ydir, ystep, yphase):
  # The input {OPHS} must be a list of raster flll elements.
  # Splits it into sub-groups consisting of at most {max_lines}
  # scanlines each.  Assumes that group indices {0..ngr-1} have been
  # assigned, and renumbers the new groups starting at {ngr}
  # Returns the updated number of groups {ngr}.
  
  debug = False
  nph_old = len(OPHS)

  SCS = raster.separate_by_scanline(OPHS, xdir, ydir, ystep, yphase)
  nsc_old = len(SCS) # Number of scanlines in group
  if debug: sys.stderr.write("  splitting group with %d scanlines (%d rasters): " % (nsc_old,nph_old))

  ngr_old = ngr
  ngr_new = ngr
  if nsc_old > 0:
    nph_new = 0      # Raster elements processed.
    nsc_new = 0      # Number of scanlines processed.
    isc = 0          # Index of next scanline to process.
    while isc < nsc_old > 0:
      nsc_rm = nsc_old - isc # Number of scanlines left to split
      # Decide number of subgroups {ngr_rm} in remaining scanlines:
      ngr_rm = (nsc_rm + max_lines - 1)//max_lines
      # Decide number of scanlines {nsc_gr} in next group:
      nsc_gr = int(floor(nsc_rm/ngr_rm + 0.5))
      assert nsc_gr <= max_lines
      if debug: sys.stderr.write(" %d (" % nsc_gr)

      # Unite scanlines in next chunk:
      isc_lim = min(nsc_old, isc + nsc_gr)
      while isc < isc_lim:
        # Assign fillers in scanline to group {ngr_new}:
        if debug: sys.stderr.write(" [%d:%d: " % (isc,len(SCS[isc])))
        for iph in SCS[isc]:
          oph = OPHS[iph]
          path.set_group(oph, ngr_new)
          nph_new += 1
          if debug: sys.stderr.write("%s," % path.get_name(oph))
        if debug: sys.stderr.write("]")
        isc += 1
        nsc_new += 1
      if debug: sys.stderr.write(")")
      # One more new group defined:
      ngr_new += 1
      # Discount the scanlines processed:
      nsc_rm = nsc_rm - nsc_gr
  if debug: sys.stderr.write("\n")
  ngr_plus = ngr_new - ngr_old
  if debug: 
    sys.stderr.write("group was split into %d groups %d..%d" % (ngr_plus,ngr_old,ngr_new-1))
    sys.stderr.write(" (%d paths in %d scanlines)" % (nph_new, nsc_new))
  assert nsc_new == nsc_old
  assert nph_new == nph_old
  return ngr_new
  # ----------------------------------------------------------------------

