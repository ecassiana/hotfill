# Bandpath construction for the HotFill algorithm.

import bandpath_IMP
  
def build(OPHS, SCS, z, mp_jump, quick):
  # This procedure tries to build a {z}-type bandpath for some
  # {(i,j)}-band of a solid raster fill, that honors the cooling time
  # limits of internal contacts.
  #
  # The result of this procedure will be a quintuple
  # {bph,CTS_lo,TCVS_lo,CTS_hi,TCVS_hi} where {bph} is the requested
  # bandpath, {CTS_lo,CTS_hi} are the contacts on the bottom and top
  # edges of the band, and {TCVS_lo,TCVS_hi} are the cover times by {bph}
  # of those contacts. If it cannot or will not build such a path, the
  # procedure returns {None,None,...,None} instead.
  # 
  # The paranter {OPHS} should be a list of raster elements for the
  # filling. Each element should be an oriented path consisting of a
  # single horizontal trace. The path should be oriented left to right.
  # These raster elements are assumed to lie on a number of uniformly
  # spaced horizontal scan-lines. All raster traces are assumed to have
  # the same {MoveParms} record, and the nominal width is assumed to be
  # equal to the scanline spacing.
  # 
  # The parameter {SCS} should be a list of {j-i} lists of integers,
  # such that the raster with index {irs} on scanline {isc} (counted
  # from 0 at bottom of the band) is {OPHS[SCS[isc][irs]]}.
  #
  # Each raster element {oph} in {OPHS} should have a set of zero or
  # more link paths, accessed through {path.get_links(oph)}, that all
  # end at {pini(oph)}; and another set, accessed by
  # {path.get_links(rev(oph))}, of link paths that end at {pfin(oph)}.
  # The reverse of each of those link paths should be a link path of
  # some other raster trace, usually on a scan-line adjacent to that o
  # {oph}.
  #
  # Each raster element {oph} in {OPHS} should also have a set of zero
  # or more relevant {Contact} objects, accessed through
  # {path.get_contacts}. In each contact, side 0 should be the lower
  # trace, and side 1 should be the upper one. Thus
  # {path.get_contacts(oph,0)} should be all relevant contacts on the
  # upper side of the raster {oph}, and {path.get_contacts(oph,1)}
  # should be those on the lower side of {oph}.
  #
  # The bandpath {bph} returned by this procedure will includes every
  # trace indentified in {SCS}, possibly reversed. The gap between any
  # two successive elements will be filled by a link path with the
  # proper endpoints, or by a jumps with {Move_Parms} parameters
  # {mp_jump} if there is no such link
  #
  # The path can be of two types, selected by the parameter {z}. If {z}
  # is 0, the path {bph} will begin with the leftmost raster on the bottom
  # scan-line of the band and end with the rightmost raster on the top
  # scan-line, both oriented left-to-right.
  # 
  # If {z} is 1, the path {bph} will begin with the rightmost raster on
  # the bottom scan-line of the band and end with the leftmost raster on
  # the top scan-line, both oriented right-to-left.
  #
  # The result {CTS_lo} will be a list of the contacts on the bottom
  # edge of the band (cut-line {i}), in left-to-right order. The result
  # {TCVS_lo} will be a list with the cover times
  # {contact.path_tcov(bph,imv,ct,1)} of those contacts, in the same
  # order.
  #
  # The result {CTS_hi} will be a list of the contacts on the top edge
  # of the band (cut-line {j}), also in left-to-right order. The result
  # {TCVS_hi} will be a list with the cover times
  # {contact.path_tcov(path.rev(bph),imv,ct,1)} of those contacts, in the
  # same order.
  #
  # All the cover times in {TCVS_lo} and {TCVS_hi} as well as the
  # cooling times of all internal contacts in the path {bph} will be less
  # than or equal to the cooling time limits of those contacts. The
  # procedure will return {None,None,...,None} if it can't build a path
  # {bph} that satisfies these constraints.
  # 
  #
  # If {quick} is false, the procedure uses the greedy
  # algorithm described in the paper, with the meet-in-the-middle
  # optimization.  The running time will be is quadratic on the number
  # of rasters.
  #
  # If {quick} is true, the procedure will use another algorithm that
  # runs in time proposrtinal to the number of rasters.  
  # However, the resulting bandpath may have bigger fabtime,
  # and the procedure may even return {None} when the other 
  # version would return a possible bandpath.  Or vice-versa,
  # it is not cerain.
  return bandpath_IMP.build(OPHS, SCS, z, mp_jump, quick)
