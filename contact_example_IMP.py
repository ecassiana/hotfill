# Implementation of module {contact_example}
# Last edited on 2021-10-31 05:30:41 by stolfi

import contact_example
import move
import move_parms
import path
import path_example
import block
import block_example
import contact
import job_parms
import hacks
import rn
import pyx 
from math import sqrt, sin, cos, log, exp, nan, inf, pi
import sys

parms = job_parms.typical_js()
parms['solid_raster_width'] = 1.00
parms['contour_trace_width'] = 0.50

mp_jump = move_parms.make_for_jumps(parms)
mp_cont = move_parms.make_for_contours(parms)
mp_fill = move_parms.make_for_fillings(parms)

wd_fill = move_parms.width(mp_fill)
wd_cont = move_parms.width(mp_cont)

# CONTACTS

def misc_A(mp_fill):
  
  # Create a contact {ct}:
  p00 = (1,1)
  p01 = (4,1)
  mv0 = move.make(p00, p01, mp_fill)
  
  p10 = (2,2)
  p11 = (5,2)
  mv1 = move.make(p10, p11, mp_fill)
  
  q0 = (2, 1.5)
  q1 = (4, 1.5)
  ct = contact.make(q0, q1, mv0, mv1)
  contact.set_name(ct, "CA")

  return ct
  # ----------------------------------------------------------------------

def misc_B(mp_trace, mp_jump):

  # Create five paths with ten traces and five jumps:
  OPHS, TRS, JMS = path_example.misc_E(mp_trace, mp_jump)
  assert TRS != None
  
  # Make some contacts between the traces:
  q00 = (2.0, 3.5)
  q01 = (3.0, 3.5)
  ct0 = contact.make(q00, q01, TRS[0], TRS[1])

  q10 = (4.0, 3.5)
  q11 = q10
  ct1 = contact.make(q10, q11, TRS[1], TRS[2])

  q20 = (3.5, 3.0)
  q21 = (3.5, 3.0)
  ct2 = contact.make(q20, q21, TRS[0], TRS[2])

  dq3 = (1.0/sqrt(5), 0.5/sqrt(5))
  q30 = tuple(rn.sub(move.pini(TRS[3]), dq3))
  q31 = tuple(rn.add(move.pfin(TRS[2]), dq3))
  ct3 = contact.make(q30, q31, TRS[2], TRS[3])

  q40 = (5.0, 4.5)
  q41 = (5.0, 4.5)
  ct4 = contact.make(q40, q41, TRS[1], TRS[4])
  
  q50 = (3.0, 7.5)
  q51 = (4.0, 7.5)
  ct5 = contact.make(q50, q51, TRS[7], TRS[9])
  
  q60 = (5.0, 7.5)
  q61 = (6.0, 7.5)
  ct6 = contact.make(q60, q61, TRS[8], TRS[9])
  
  q70 = (3.0, 6.5)
  q71 = (4.0, 6.5)
  ct7 = contact.make(q70, q71, TRS[7], TRS[5])
  
  CTS = [ ct0, ct1, ct2, ct3, ct4, ct5, ct6, ct7, ]
  for ict in range(len(CTS)): contact.set_name(CTS[ict], "CB%d" % ict)
   
  return CTS, OPHS, TRS
  # ----------------------------------------------------------------------
  
def misc_C(mp_trace):
  mp_trace = mp_fill
  wd = move_parms.width(mp_trace)
  nrs = [7, 3]        # Number of rasters in each path.
  plo = [None,None] # Lower left corner of bbox each path.
  sz = [ 3*wd, 5*wd]
  plo[0] = ( 2, 1 )
  plo[1] = ( plo[0][0] + sz[0] + wd, 2)

  PHS = [None, None]  # The two paths.
  for iph in range(2):
    axis = iph
    alt = True
    PHSi, TRSi, LJS0i, LJS1i = path_example.raster_rectangle(plo[iph], axis, nrs[iph], alt, sz[iph], wd, mp_trace,None)
    move.tag_names(TRSi, "s%d." % iph)
    move.tag_names(LJS0i, "s%d." % iph)
    move.tag_names(LJS1i, "s%d." % iph)
    path.tag_names(PHSi, "s%d." % iph)
    PHS[iph] = PHSi[0]
  szmin = 0.25
  rszmin = 0.50
  tol = 0.20*wd  # Tolerance for overlaps, tilts, etc.
  CTS = contact.from_paths(PHS[0], PHS[1], szmin, rszmin, tol, None)
  for ict in range(len(CTS)): contact.set_name(CTS[ict], "CC%d" % ict)

  return CTS, PHS
  # ----------------------------------------------------------------------

def misc_F(alt, mp_trace, mp_jump):

  ph = path_example.misc_F(alt, mp_trace, mp_jump)
  
  # Get the adjacent traces, which should be horizontal:
  if alt:
    ixmv = (5, 15)
  else:
    ixmv = (3, 9)
  mv = [None, None] # The two moves.
  plo = [None, None] # The left endpoints of the moves.
  phi = [None, None] # The right endpoints of the moves.
  for kmv in range(2):
    mvk, drk = move.unpack(path.elem(ph, ixmv[kmv]))
    plok = move.pini(mvk)
    phik = move.pfin(mvk)
    assert plok[1] == phik[1]
    if plok[0] > phik[0]: plok, phik = phik, plok
    mv[kmv] = mvk; plo[kmv] = plok; phi[kmv] = phik
  
  # Get the endpoints of the contact:
  y = (plo[0][1] + plo[1][1])/2
  xlo = max(plo[0][0], plo[1][0])
  xhi = min(phi[0][0], phi[1][0])

  ct = contact.make((xlo,y), (xhi,y), mv[0], mv[1])
  contact.set_name(ct, "CF")

  return ct, ph
  # ----------------------------------------------------------------------

def misc_K(mp_trace, mp_jump):

  wd = move_parms.width(mp_trace)
  
  # Create a path with an internal contact:
  pas = ( (0,0), (2,0), (5,0), (9,0), (9,9), (8,9), (8,7), (6,7), (1,7), (0,7), )
  Pa = path.from_points(pas, mp_trace, mp_jump)
  path.set_name(Pa, "Pa", True)
  assert path.nelems(Pa) == len(pas)-1

  # Get the two moves of {Pa} that make the internal contact:
  omva3 = path.elem(Pa, 3); assert move.pini(omva3) == (9,0)
  omva5 = path.elem(Pa, 5); assert move.pini(omva5) == (8,9)
  # Get two other moves of {Pa} that make contacts with the other two paths:
  omva1 = path.elem(Pa, 1); assert move.pini(omva1) == (2,0)
  omva7 = path.elem(Pa, 7); assert move.pini(omva7) == (6,7)
  
  # Create another path that makes two contacts with {Pa}:
  pbs = ( (0,2), (2,1), (5,1), (1,6), (3,6), )
  Pb = path.from_points(pbs, mp_trace, mp_jump)
  path.set_name(Pb, "Pb", True)
  assert path.nelems(Pb) == len(pbs)-1
  
  # Get the two moves that make contact with {Pa}
  omvb1 = path.elem(Pb, 1); assert move.pini(omvb1) == (2,1)
  omvb3 = path.elem(Pb, 3); assert move.pini(omvb3) == (1,6)

  # Make another path {Pc} that shares one of those moves:
  pcs = ( (2,2), (2,1), (5,1), (5,4), )
  omvc0 = move.make(pcs[0], pcs[1], mp_trace)
  omvc2 = move.make(pcs[2], pcs[3], mp_trace)
  Pc = path.from_moves((omvc0, omvb1, omvc2,))
  path.set_name(Pc, "Pc", True)
  
  OPHS = [ Pa, Pb, Pc ]
  
  # Create contacts:
  
  def mkcontact(omv0, omv1, name):
    szmin = 0.25
    rszmin = 0.01
    e0 = move.pini(omv0)
    e1 = move.pfin(omv0)
    mvdir, mlen = rn.dir(rn.sub(e1, e0))
    tol = 0.20*wd  # Tolerance for overlaps, tilts, etc.
    ct = contact.from_moves(omv0, omv1, mvdir, szmin, rszmin, tol)
    contact.set_name(ct, name)
    return ct
    # ......................................................................

  ctA = mkcontact(omva3, omva5, "ctA")
  ctB = mkcontact(omva1, omvb1, "ctB")
  ctC = mkcontact(omva7, omvb3, "ctC")

  # Add the path-contact information:
  def set_path_contact_info(P, imv, ct, isd):
    path.add_contact(P, isd, ct)
    contact.add_side_path(ct, isd, P, imv)
    contact.add_side_path(ct, isd, path.rev(P), path.nelems(P)-1-imv)
    
  set_path_contact_info(Pa, 3, ctA, 0)
  set_path_contact_info(Pa, 5, ctA, 1)
                              
  set_path_contact_info(Pa, 1, ctB, 0)
  set_path_contact_info(Pb, 1, ctB, 1)
  set_path_contact_info(Pc, 1, ctB, 1)
                              
  set_path_contact_info(Pa, 7, ctC, 0)
  set_path_contact_info(Pb, 3, ctC, 1)
  
  CTS = [ ctA, ctB, ctC, ]
  
  return CTS, OPHS
  # ----------------------------------------------------------------------

# FROM PATHS AND BLOCKS

def raster_raster_contact(oph0, oph1):
  # The {X} range of the contact will be {[xlo _ xhi]}
  xlo = -inf 
  xhi = +inf
  oph = (oph0, oph1)
  mv = [None,None]    # Raster traces at top of {oph[0]} and bottom of {oph[1]}.
  imv = [None,None]   # {imv[kph]} is the index of {mv[kph]} in {oph[k]}.
  
  for kph in range(2):
    dir = 1 - 2*kph
    mv[kph], imv[kph] = get_check_end_move(oph[kph], None, dir)
    
    # Find the {X} span and intersect with current {[xlo _ xhi]}
    pk, qk = move.endpoints(mv[kph])
    assert pk[1] == qk[1] # Trace must be horizontal.
    xlok = min(pk[0], qk[0])
    xhik = max(pk[0], qk[0])
    xlo = max(xlo, xlok)
    xhi = min(xhi, xhik)

  # Check coordinates:
  y0 = move.pini(mv[0])[1]; wd0 = move.width(mv[0])
  y1 = move.pini(mv[1])[1]; wd1 = move.width(mv[1])
  assert abs((y0 + wd0/2) - (y1 - wd1/2)) < 0.01*(wd0+wd1), "traces are not adjacent"
  assert xlo < xhi, "{X} randges do not overlap"
  
  # Create the contact:
  ymd = (y0 + y1 + (wd0 - wd1)/2)/2
  ct = contact.make((xlo,ymd), (xhi, ymd), mv[0], mv[1])
  
  # Set the path-contact info:
  contact.add_side_path(ct, 0, oph0, imv[0])
  contact.add_side_path(ct, 1, oph1, imv[1])
  
  path.add_contact(oph0, 0, ct)
  path.add_contact(oph1, 1, ct)
  
  return ct
  # ----------------------------------------------------------------------

def two_roads_and_islands(nmv, nmg, nis, mp_cont, mp_fill, mp_jump):

  wdc = move_parms.width(mp_cont)
  wdf = move_parms.width(mp_fill)
  
  org = (2,2)   # Lower left corner of whole thing.

  # Key dimensions and coords, in {wdf}  units, rel to {org}:
  nx_rd = 5                 # Length of rasters on the roads.
  nxy_is = 4                # Width and height of an island plus margin.
  
  nx_gp = nxy_is + 2        # Width of gap between roads.

  ixstep_rd = nx_rd + nx_gp # Scan-column increment between roads.
  ixstep_is = ixstep_rd     # Scan-column increment between island centers.
  
  assert nx_gp % 2 == 0     # To avoid half-scanlines.
  ixctr_is0 = nx_gp//2      # Scan-column of center of left island group.
  ix0_rd0 = nx_gp           # Left scan-column of left road.

  nytot_is = nis*nxy_is                       # Total height of islands plus margins.
  assert nytot_is % 2 == 0                    # To avoid half-scanlines.
  iyctr_is = max(nmv//2, nytot_is//2)         # Index of scanline at center of islands.
  iyctr_is0 = iyctr_is - (nytot_is-nxy_is)//2 # Index of scanline at center of bottom island.
  iy0_rd = iyctr_is - nmv//2                  # Index of scanline at bottom of roads. 
  iystep_is = nxy_is                          # Scanline step between centers of islands.

  # Absolute road dimensions and coords in {mm} (do not include the width of traces):
  Ric = wdf + wdc               # Radius of endpoints in island contour.
  assert Ric + wdc/2 < nxy_is*wdf/2    
  Rif = 0                       # Min inradius of island.
    
  def make_island_block(kx, ky):
    # Returns block that is a single island, specifically
    # island {ky} from bottom to top in group {kx} from left to right.
    # See {block_example.onion}.
    dx0 = (ixctr_is0 + kx*ixstep_is)*wdf
    dy0 = (iyctr_is0 + ky*iystep_is)*wdf
    ctr = rn.add(org, (dx0, dy0))
    nch = 4  # Number of choices of starting point of contour.
    bc = block_example.onion(nch, ctr, Ric, mp_cont, Rif, mp_fill, mp_jump)
    return bc

  BCS = []
  CTS = []
  
  # Compute the number {nbc} of multiraster blocks at each end of a road:
  nbc = ((nmv - 1)//(2*nmg))
  assert nbc >= 0
  for ird in range(2):
    # Decide range of scan-columns {ix0..ix1} of this road:
    ix0 = ix0_rd0 + ird*ixstep_rd
    ix1 = ix0 + nx_rd
    BCSi, CTSi = road_of_blocks(nmv, nmg, org, ix0,ix1, iy0_rd, mp_fill)
    BCS += BCSi
    CTS += CTSi
    
  # Add island blocks:
  for kx in range(3):
    for ky in range(nis):
      bc_this = make_island_block(kx, ky)
      BCS.append(bc_this)
  return BCS, CTS
  # ----------------------------------------------------------------------

def get_check_end_move(oph, mv_check, dir):
  # Given a snake path {oph} of horizontal rasters, gets the {Move}
  # object {mv} of the raster of {oph} that has the highest Y if {dir}
  # is {+1}, or the lowest Y if {dir} is -1, independently of the
  # orientation of {oph}. Also returns the index {imv} of that move in
  # {oph} Assumes that it is either the first or the last move of {oph}.
  # 
  # If {mv_check} is not None, bombs out if {mv} is not equal to
  # {mv_check}.
  
  nmv = path.nelems(oph)
  mv0, dr0 = move.unpack(path.elem(oph,0))
  mv1, dr1 = move.unpack(path.elem(oph,nmv-1))
  dy = move.pini(mv1)[1] - move.pini(mv0)[1]
  mv, imv = (mv1, nmv-1) if dy*dir >= 0 else (mv0, 0)
  if mv_check != None: assert mv == mv_check
  return mv, imv
  # ----------------------------------------------------------------------
  
def road_of_blocks(nmv, nmg, org, ix0,ix1, iy0, mp_fill): 

  wdf = move_parms.width(mp_fill)
  szmin = 0.9*wdf
  rszmin = 0.0

  nbc = ((nmv - 1)//(2*nmg))  # Number of snake blocks at each end of road.
  nbt = 2*nbc + (nmv - 2*nbc*nmg) # Total number of blocks in road.

  sys.stderr.write("iy0 = %d nmv = %d nmg = %d nbc = %d\n" % (iy0,nmv,nmg,nbc))
  
  BCS = []
  CTS = []
  bc_prev = None # Previous block in same road.
  iy1_prev = iy0 - 1       # Scanline index of top of previous block.
  for ib in range(nbt):
    sys.stderr.write("iy1_prev = %d\n" % iy1_prev)
  
    # Decide number of raster traces in this block:
    nmvi = nmg if ib < nbc or ib >= nbt - nbc else 1

    # Decide scan line span {iy0_this..iy1_this} for block {ib}:
    iy0_this = iy1_prev + 1
    iy1_this = iy0_this + nmvi - 1
    assert iy0_this <= iy1_this and iy1_this < iy0 + nmv

    bc_this = make_road_block(org, ix0,ix1, iy0_this,iy1_this, mp_fill)
    BCS.append(bc_this)
    if bc_prev != None:

      # Create a contact between the two blocks {bc_prev,bc_this}:
      mv_prev, imv_prev = get_check_end_move(block.choice(bc_prev,0), None, +1);  # Top {Move} of {bc_prev}
      mv_this, imv_this = get_check_end_move(block.choice(bc_this,0), None, -1);  # Bottom {Move} of {bc_this}
      mvdir = (1,0)
      tol = 0.20*wdf  # Tolerance for overlaps, tilts, etc.
      ct = contact.from_moves(mv_prev, mv_this, mvdir, szmin, rszmin, tol) # Contact between the two blocks.
      assert ct != None
      CTS.append(ct)

      # Set the path-contact links:
      for isd in range(2):
        dir = 1 - 2*isd # {+1} for side 0, {-1} for side 1.
        bc = (bc_prev,bc_this)[isd] # The block on side {isd} of contact
        mv_exp = (mv_prev,mv_this)[isd]
        for ich in range(block.nchoices(bc)):
          ophi = block.choice(bc,ich)
          mvi, imvi = get_check_end_move(ophi, mv_exp, dir);  # Should be the move of {ophi} on side {isd} of {ct}.
          contact.add_side_path(ct, isd, ophi, imvi)
          path.add_contact(ophi, isd, ct)

    bc_prev = bc_this
    iy1_prev = iy1_this
  sys.stderr.write("iy1_prev = %d\n" % iy1_prev)
  assert iy1_prev == iy0 + nmv - 1
  return BCS, CTS
  # ----------------------------------------------------------------------

def make_road_block(org, ix0,ix1, iy0,iy1, mp_fill):
  # Returns a block that is a serpentine raster path
  # spanning scan-lines with numbers {iy0} to {iy1} and scan-columns {ix0} to {ix1},
  # assuming scan-lines and scan-columns are spaced {wdf=with(mp_fill)} apart
  # and the indices starts at the point {org}.
  # The block has {iy0=iy1}, the block is a single raster line. See 
  # {block_example.raster_rectangle}.
  sys.stderr.write("make_road_block(%s,%s, %s,%s)\n" % (ix0,ix1,iy0,iy1))
  wdf = move_parms.width(mp_fill)
  plo = rn.add(org, (ix0*wdf, iy0*wdf))
  nx_rd = ix1 - ix0 + 1
  ny_rd = iy1 - iy0 + 1
  hor = True
  ver = False
  alt = True
  bc = block_example.raster_rectangle(plo, nx_rd, ny_rd, hor, ver, alt, mp_fill, None)
  return bc
  # ----------------------------------------------------------------------

