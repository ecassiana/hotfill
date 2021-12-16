# Implementation of the module {bandpath}.
# Last edited on 2021-11-25 15:11:14 by stolfi

import bandpath
import path
import contact
import move
import move_parms
import hacks
import rn
import sys
from math import sqrt, sin, cos, log, exp, floor, ceil, inf, nan, pi

def build(OPHS, SCS, z, mp_jump, quick):
  
  # Count the rasters in the band:
  nrs = 0
  for Si in SCS: nrs += len(Si)
  
  # Get the lists of moves {BPHS0,BPHS1} that will becomethe two halves of the bandpath:
  if quick:
    BPHS0, BPHS1 = get_bandpath_elements_quick(OPHS, SCS, z, mp_jump)
  else:
    BPHS0, BPHS1 = get_bandpath_elements_greedy(OPHS, SCS, z, mp_jump)
  n0 = len(BPHS0); n1 = len(BPHS1)
  assert n0 + n1 == nrs
  
  use_links = True
  assert n0 > 0
  bph0 = path.concat(BPHS0, use_links, mp_jump)
  if n1 > 0:
    bph1 = path.concat(BPHS1, use_links,mp_jump)
    bph = path.concat([bph0, path.rev(bph1)], use_links, mp_jump)
  else:
    bph = bph0
  
  if check_internal_contacts(bph, BPHS0, BPHS1): 
    CTS_lo = get_contacts([OPHS[iph] for iph in SCS[0]], 1)  # Contacts on cut-line {i}
    CTS_hi = get_contacts([OPHS[iph] for iph in SCS[-1]], 0) # Contacts on cut_line {j}

    TCVS_lo = get_cover_times(bph, CTS_lo, 1, z)
    TCVS_hi = get_cover_times(path.rev(bph), CTS_hi, 0, z)
    if TCVS_lo == None or TCVS_hi == None: 
      # Edge cover times too large, will violate cooling constraints:
      bph = None
  else:
    # Internal cooling constraints violated:
    bph = None
      
  if bph == None:
    CTS_lo = None; TCVS_lo = None; CTS_hi = None; TCVS_hi = None
  
  return bph, CTS_lo, TCVS_lo, CTS_hi, TCVS_hi
  # ----------------------------------------------------------------------

def get_contacts(OPHS, isd):
  # Given a list {OPHS} of raster fill elements, returns a list with every contact
  # {ct} that has one of those elements as side {isd}. Uses {path.get_contacts}.
  # The contacts will be sorted by the horizontal coordinates of their midpoints.
  CTS = []
  for oph in OPHS:
    for ct in path.get_contacts(oph, isd):
      CTS.append(ct)
  CTS.sort(key = lambda ct: contact.pmid(ct)[0])
  return CTS
  # ----------------------------------------------------------------------
  
def get_cover_times(bph, CTS, isd, dir):
  # Given a list {CTS} of contacts, returns a list with {contact.path_tcov(bph,imv,ct,isd)}
  # for each contact {ct} in that list. The path must cover side {isd} of every contact.
  # 
  # If {dir=0}, the procedure expects that the path will cover the contacts in roughly the same order
  # as they occur in {CTS}. If {dir=1}, it expects them to be covered in the opposite order.
  # These conditions are not critical, but the procedure will be a bit faster if
  # they are satisfied.
  assert dir == 0 or dir == 1 
  nmv = path.nelems(bph)
  assert nmv  >= 1
  
  TCVS = []
  imv_guess = 0 # Initial guess of move that covers the side of the next contact
  for ct in CTS:
    mv_ct = contact.side_move(ct, isd)
    # Find index {imv} of that move in {bph}:
    imv = None
    for kmv in range(nmv):
      mv_ph, dr_ph = move.unpack(path.elem(bph, imv_guess))
      if mv_ph == mv_ct: imv = imv_guess; break
      imv_guess += (+1 if dir == 0 else -1)
      imv_guess = imv_guess % nmv
      assert imv_guess >= 0
    assert imv != None, "{bph} does not cover {ct}!?"
    # Get the cover time:
    tcv = contact.path_tcov(bph, imv, ct, isd)
    if tcv > contact.tcool_limit(ct): return None
    TCVS.append(tcv)
    
  return TCVS
  # ----------------------------------------------------------------------

def check_internal_contacts(bph, BPHS0, BPHS1):
  # Checks the cooling time limits of contacts that are closed by a bandpath {bph}.
  #
  # The arguments {BPHS0} and {BPHS1} should be lists of paths, which are
  # original raster filling elements or their reverses.  The bandpath
  # {bph} is supposed to be {concat(concat(BPHS0), rev(concat(BPHS1)))}.
  # The relevant contacts covered by an element {oph} of either list
  # should be provided by {path.contacts(oph)}.

  debug = False
  nmv = path.nelems(bph)
  n0 = len(BPHS0); n1 = len(BPHS1)
  assert nmv >= 2*(n0 + n1) - 1 # There must be at leasr 1 link trace or jump between two rasters.
  
  # Initialize to {-1} the index hints of all contacts (intra-band or inter-band) of the raster elemets:
  for bphi in BPHS0 + BPHS1:
    for isd_rs in 0,1:
      for ct in path.get_contacts(bphi, isd_rs):
        for isd_ct in 0,1:
          contact.set_side_index(ct, isd_ct, -1)

  # Collect a list {CTS} of all contacts on the upper sides of all rasters of {bph} 
  # and saves the *guessed* indices of the sides of each contact {ct} on the path {bph}
  # as {contact.get_side_index(ct,isd)}. 
  #
  # That is, side {isd} of {ct} should be {path.elem(bph,imv)} where
  # {imv} should be {contact.get_side_index(ct,isd)} or somewhat ahead
  # of it, if the latter is not {-1}.
  CTS = []
  for irs in range(n0+n1):
    # Get the *raster element* of index {irs} in {bph}. Note that it is 
    # not the *move* of index {irs}, because of links and jumps.
    bphi = BPHS0[irs] if irs < n0 else BPHS1[n0+n1-1-irs]
    for isd in 0,1:
      for ct in path.get_contacts(bphi,isd):
        if isd == 0: CTS.append(ct)
        imv_guess = 2*irs # A lower bound to the index of the move.
        contact.set_side_index(ct, isd, imv_guess)

  # Now find the corrext positions of the sides of each contact in 
  # the path {bph}, and check the cooling time limits of all contacts closed by the path:
  for ict in range(len(CTS)):
    ct = CTS[ict]
    if debug: contact.show(sys.stderr, ("  CTS[%03d] = " % ict), ct, "\n", 0)
    # Get the cover times of {ct} by {bph}
    tcv = [ None, None ] # Cover times of each side by {bph}, or {None}
    for isd in 0,1:
      # Get the index {imv} of side {isd} of {ct} in {bph}.  Beware of garbage:
      imv_guess = contact.get_side_index(ct, isd)
      if debug: sys.stderr.write("    isd = %s index = %s\n" % (isd,imv_guess))
      if imv_guess != -1:
        assert imv_guess >= 0 and imv_guess < nmv
        imv = imv_guess
        while imv < nmv:
          omv = path.elem(bph, imv);
          mv, dr = move.unpack(omv)
          if mv == contact.side_move(ct, isd): break
          imv += 1
        assert imv < nmv
        tcv[isd] = contact.path_tcov(bph, imv, ct, isd)
      else:
        # This side of {ct} is not covered by {bph}:
        tcv[isd] = None
    # Path must cover at least one side:
    assert tcv[0] != None or tcv[1] != None

    # Check cooling time limit:
    if tcv[0] != None and tcv[1] != None:
      tcool = abs(tcv[0] - tcv[1])
      if tcool > contact.tcool_limit(ct): return False
        
  # Passed all checks:
  return True
  # ----------------------------------------------------------------------

# QUICK VERSION:

def get_bandpath_elements_quick(OPHS, SCS, z, mp_jump):
  # Given the list of all raster elements {OPHS} and a list of lists {SCS} 
  # of indices into it, returns a list with the raster elements specified in {SCS}
  # sorted and re-oriented as in a {z}-type bandpath.
  #
  # The elements of {OPHS} and {SCS} should be the same of {build}.
  
  # The result is two lists {BPHS0} and {BPHS1} of raster elements
  # such that the final bandpath will be the concatenation of all
  # elements of {BPHS0} and the reverse of the concatenation of those
  # of {BPHS0}.
  
  nsc = len(SCS) # Number of scan-lines in band.
  
  BPHS0 = []
  BPHS1 = []
  
  # The arrays {rlo} and {rhi} keep track of the still-unused rasters.
  # See {get_band_stack_pair}.
  rlo = [ 0 ] * nsc
  rhi = [ len(S)-1 for S in SCS ]
  
  failed = False  # No more rasters to take?
  v = 0           # Vertical order of the next stack of rasters (0 = bottom to top, 1 = vv).
  while not failed:
    failed = shave_band_stack_pair(OPHS, SCS, rlo, rhi, BPHS0, BPHS1, z, v, mp_jump)
    v = 1-v
    
  return BPHS0, BPHS1
  # ----------------------------------------------------------------------
  
def shave_band_stack_pair(OPHS, SCS, rlo, rhi, BPHS0, BPHS1, z, v, mp_jump):
  # Peels off from the unused part of {SCS} two partial stacks and appends
  # them to {BPHS0,BPHS1}. Also updates {rlo,rhi}.
  #
  # Assumes that the first half of the bandpath will be the rasters of
  # {BPHS0}, in that orider and orientation; while the second half will
  # be the rasters of {BPHS1}, in reverse order and all reversed.
  #
  # The list {OPHS} must have all the raster paths that were given to
  # {HotFill}, all oriented left-to-right.
  #
  # The list of lists {SCS} identifies the raster elements in the band.
  # Each element {SCS[k]} is a list of indices into {OPHS} of the
  # rasters on scan-line {k} of the band, with {k=0} being the lowest
  # scanline in the band. The raster element on that scan-line with
  # index {r}, from left to right, is {OPHS[SCS[k][r]]}.
  #
  # The arrays {rlo} and {rhi} keep track of the rasters in each
  # scan-line of {SCS} that have not been used yet in {BPHS0} or {BPHS1}.
  # Namely, the unused elements of {SCS[k]} are {SCS[k][u]} for {u} in
  # {rlo[k]..rhi[k]}. ]
  #
  # The parameter {z} tells the order of the stacks in the bandpath, and
  # also the orientation of the first raster in each stack:
  #
  #   If {z=0}, the leftmost raster in each scan-line is appended to
  #   {BPHS0}, and the righmost one is appended to {BPHS1}. The first
  #   raster added to {BPHS0} will be oriented left-to-right, and
  #   the first one added to {BPHS1} will be oriented right-to-left.
  #
  #   if {z=1}, the right raster in each scan-line will be added to {BPHS0},
  #   and the leftmost one to {BPHS1}. The first addition to {BPHS0}
  #   will be oriented right-to-left, and the first one to {BPHS1} will be
  #   oriented left-to-right
  #
  # The parameter {v} tells the order of the rasters in each stack.
  # If {v=0} the first raster to be added to {BPHS0} is the lowest one of the
  # stack, otherwise it is the higest one; and the opposite holds
  # for the first raster to be added to {BPHS1}.
  #
  # The procedure returns {True} if it failed to add a single raster
  # to {BPHS0} or {BPHS1}, because all rasteds listed in {SCS} were used already.
  
  def reorient(oph, BPHSX, n):
    # Returns the raster element {oph} re-oriented as appropriate before
    # adding it to the list {BPHSX} (either {BPHS0} or {BPHS1}). Assumed
    # that {n} rasters from the current stack have been added to
    # {BPHSX}. The raster {oph} should be oriented left--to-right for
    # {BPHS0} and right-to-left for {BPHS1}.
    
    if n == 0:
      # Modify the orientation as specified by {z}:
      oph = path.spin(oph, z)
    else:
      # Choose the orientation with smaller connection fabtime from previous raster:
      use_links = True
      t0 = path.connection_time(BPHSX[-1], oph, use_links, mp_jump)
      t1 = path.connection_time(BPHSX[-1], path.rev(oph), use_links, mp_jump)
      if t1 < t0: oph = path.rev(oph)
    return oph
    # ....................................................................

  nsc = len(SCS)
  n0 = 0 # Rasters added to {BPHS0}.
  n1 = 0 # rasters added to {BPHS1}.
  
  for ksc in range(nsc):
    # Decide the scanlines {isc0,isc1} where the next raster of {BPHS0} and {BPHS1} come from:
    if v == 0:
      isc0 = ksc; isc1 = nsc-1-ksc
    else:
      isc0 = nsc-1-ksc; isc1 = ksc
    
    # Grab the next raster for {BPHS0}:
    if rlo[isc0] <= rhi[isc0]:
      if z == 0:
        # Take from the left end of scanline {isc0}:
        oph0 = OPHS[SCS[isc0][rlo[isc0]]]; rlo[isc0] += 1
      else:
        # Take from the right end of scanline {isc0}:
        oph0 = OPHS[SCS[isc0][rhi[isc0]]]; rhi[isc0] -= 1
      oph0 = reorient(oph0, BPHS0, n0)
      BPHS0.append(oph0)
      n0 += 1
    else:
      oph0 = None
    
    # Grab the next raster for {BPHS1}:
    if rlo[isc1] <= rhi[isc1]:
      if z == 0:
        # Take from the right end of scanline {isc1}:
        oph1 = OPHS[SCS[isc1][rhi[isc1]]]; rhi[isc1] -= 1
      else:
        # Take from the left end of scanline {isc1}:
        oph1 = OPHS[SCS[isc1][rlo[isc1]]]; rlo[isc1] += 1
      oph1 = path.rev(oph1)
      oph1 = reorient(oph1, BPHS1, n1)
      BPHS1.append(oph1)
      n1 += 1
    else:
      oph1 = None
        
  failed = (n0+n1 == 0)
  return failed
  # ----------------------------------------------------------------------

# GREEDY VERSION

def get_bandpath_elements_greedy(OPHS, SCS, z, mp_jump):
  # Same as {get_bandpath_elements_quick}, but uses the
  # meet-in-the-middle greedy method.
  #
  # The next raster at each step is selected by the Eucliean distance
  # criterion, ignoring link paths.

  nsc = len(SCS) # Number of scan-lines in band.
  assert nsc >= 1

  BPHS0 = []
  BPHS1 = []

  SCS_used = [ [False]*len(SCS[i]) for i in range(nsc) ]

  # Next candiate for {BPHS0} is {OPHS[SCS[i0][j0]}}, with orientation {z0}:
  i0 = 0 
  j0 = 0 if z == 0 else (len(SCS[i0])-1)
  z0 = z 
  
  # Next candiate for {BPHS'} is {OPHS[SCS[i1][j1]}}, with orientation {z1}:
  i1 = nsc-1
  j1 = (len(SCS[i1])-1) if z == 0 else 0
  z1 = 1 - z

  # Endpoints of the half paths:
  q_end_0 = None; q_end_1 = None

  if i0 == i1 and j0 == j1:
    i1 = None; j1 = None; z1 = None

  while i0  != None or i1 != None:
    assert i0 != None # Since we grab it first.
    if i1 == None:
      # Only one raster left.
      # Add to {BPHS0}, but the best orientation may not be {z0}:
      if q_end_0 != None:
        assert q_end_1 != None # Since both paths must have the same length at this point.
        oph = OPHS[SCS[i0][j0]]
        # Elis's criterion hard to understand. 
        # Roll my own:
        oph_end_0 = BPHS0[-1]
        oph_end_1 = path.rev(BPHS1[-1])
        
        use_links = True
        # Jump/link time if using in native orientation:
        tc0_pre = path.connection_time(oph_end_0, oph, use_links, mp_jump)
        tc0_pos = path.connection_time(oph, oph_end_1, use_links, mp_jump)
        tc0 = tc0_pre + tc0_pos
        # Jump/link time if using in reversed orientation:
        tc1_pre = path.connection_time(oph_end_0, path.rev(oph), use_links, mp_jump)
        tc1_pos = path.connection_time(path.rev(oph), oph_end_1, use_links, mp_jump)
        tc1 = tc1_pre + tc1_pos
        # Choose the orientation with smaller total connection time:
        z0 = 0 if tc0 <= tc1 else 1
    
    if i0 != None:
      SCS_used[i0][j0] = True
      oph0 = path.spin(OPHS[SCS[i0][j0]], z0)
      BPHS0.append(oph0)
      q_end_0 = path.pfin(oph0)
    
    if i1 != None:
      SCS_used[i1][j1] = True
      oph1 = path.spin(OPHS[SCS[i1][j1]], z1)
      BPHS1.append(oph1)
      q_end_1 = path.pfin(oph1)

    i0, j0, z0 = find_nearest_raster(OPHS, SCS, SCS_used, nsc, q_end_0)
    i1, j1, z1 = find_nearest_raster(OPHS, SCS, SCS_used, nsc, q_end_1)
  
  return BPHS0, BPHS1
  # ----------------------------------------------------------------------
  
def find_nearest_raster(OPHS, SCS, SCS_used, nsc, q):
  # Among all the endpoints of all the rasters {OPHS[SCS[i][j]]} that
  # are still unused, finds the one that is nearest to the point
  # {q}. Returns the indices {i_min,j_min} of that raster {oph}, and
  # a flag {z_min} that is 0 if the endpoint is is {pini(oph)} and 1 if
  # it is {pfin(oph).
  #
  # If there are no unused rasters, returns {None,None,None}.
  
  dp2_min = None; i_min = None; j_min = None; z_min = None
  for i in range(nsc):
    for j in range(len(SCS[i])):
      if not SCS_used[i][j]:
        for z in 0, 1:
          oph = path.spin(OPHS[SCS[i][j]], z)
          p = path.pini(oph)
          dp2 = rn.dist_sqr(p, q)
          if (dp2_min == None or dp2 < dp2_min): 
            dp2_min = dp2; i_min = i; j_min = j; z_min = z

  if i_min != None: SCS_used[i_min][j_min] = True
  return i_min, j_min, z_min
  # ----------------------------------------------------------------------
  
