# Implementation of module {block_example}.
# Last edited on 2021-11-09 16:52:15 by stolfi

import block_example
import path_example
import block
import path
import move
import move_parms
import contact
import hacks
import rn
import pyx
import sys
from math import sqrt, sin, cos, acos, floor, ceil, inf, nan, pi

def single_raster(xlo, xhi, y, mp_trace):
  p = (xlo, y)
  q = (xhi, y)
  mv = move.make(p, q, mp_trace)
  ph = path.from_moves((mv,)) # Must be a trace
  path.set_name(ph, "Ph", False)
  bc = block.from_paths((ph, path.rev(ph),))
  block.validate(bc)
  block.set_name(bc, "Bsr")
  return bc

def raster_rectangle(plo, nx, ny, hor, ver, alt, mp_trace, mp_jump):
  assert nx >= 1 and ny >=1 
  assert hor or ver, "at least one of {hor} and {ver} must be true"
  
  wd = move_parms.width(mp_trace)
  
  szx = (nx-1)*wd
  szy = (ny-1)*wd
  OPHS = []
  if hor:
    assert nx >= 2
    PHSh,TRSh,LJS0h,LJS1h = path_example.raster_rectangle(plo, 0, ny, alt, szx, wd, mp_trace,mp_jump)
    move.tag_names(TRSh, "H.")
    move.tag_names(LJS0h, "H.")
    move.tag_names(LJS1h, "H.")
    path.tag_names(PHSh, "H.")
    if ny == 1:
      OPHS += [ PHSh[0], path.rev(PHSh[0]), ] 
    else:
      OPHS += [ PHSh[0], path.rev(PHSh[0]), PHSh[1], path.rev(PHSh[1]), ] 
  if ver:
    assert ny >= 2
    PHSv,TRSv,LJS0v,LJS1v = path_example.raster_rectangle(plo, 1, nx, alt, szy, wd, mp_trace,mp_jump)
    move.tag_names(TRSv, "V.")
    move.tag_names(LJS0v, "V.")
    move.tag_names(LJS1v, "V.")
    path.tag_names(PHSv, "V.")
    if nx == 1:
      OPHS += [ PHSv[0], path.rev(PHSv[0]), ] 
    else:
      OPHS += [ PHSv[0], path.rev(PHSv[0]), PHSv[1], path.rev(PHSv[1]), ] 
  
  bc = block.from_paths(OPHS)
  block.set_name(bc, "Brr")
  return bc

def onion(nch, ctr, Rc, mp_cont, Rf, mp_fill, mp_jump):
  wdc = move_parms.width(mp_cont)
  wdf = move_parms.width(mp_fill)
  dimp_max = 0.1*wdc  # Max deviation allowed from circle to polygon.
  dang_max = 2*acos(1 - dimp_max/Rc) # Max angle between polygon vertices.
  nt = max(3, 4*ceil(pi/2/dang_max)) # Number of traces in contour.
  OPHS = []
  for ich in range(nch):
    phase = 2*pi*ich/nch
    ph, Rin, Rex = path_example.onion(ctr, Rc, mp_cont, Rf, mp_fill, phase, nt, mp_jump)
    path.set_name(ph, "Pon%02d" % ich, False)
    OPHS.append(ph)
  bc = block.from_paths(OPHS)
  block.set_name(bc, "Bon")
  return bc

def spiral_rectangle(plo, szx, szy, mp_trace):
  PHS = []
  sz = [szx, szy ] # Signed sizes of next path's rectangle
  ax = 0           # Starting axis of nex path.
  pini = [ plo[0], plo[1] ]
  for kph in range(4):
    ph = path_example.spiral_rectangle(pini, sz[0], sz[1], ax, mp_trace)
    path.set_name(ph, "Psp%d" % kph, False)
    PHS.append(ph)
    # Move {pini} to end of first move:
    pini[ax] += sz[ax]
    # Invert direction for next start in same axis:
    sz[ax] = -sz[ax]
    # Swich axis:
    ax = 1 - ax
  bc = block.from_paths(PHS)
  block.set_name(bc, "Bsp")
  return bc
  # ----------------------------------------------------------------------

def misc_A(mp_trace, mp_jump):
  
  OPHS, TRS, JMS = path_example.misc_E(mp_trace, mp_jump)
  
  pha = OPHS[0]
  phb = OPHS[1]
  phc = OPHS[2]
  phd = OPHS[3]
  phe = OPHS[4]
  
  bc0 = block.from_paths((pha, path.rev(pha), phb, path.rev(phb),))
  block.validate(bc0)

  bc1 = block.from_paths([phd, phe])
  bc2 = block.from_paths([ phc,])
  
  BCS = [ bc0, bc1, bc2, ]
  for ibc in range(len(BCS)): block.set_name(BCS[ibc], "BA%d" % ibc)
  return BCS, OPHS, TRS, JMS
  # ----------------------------------------------------------------------

def misc_B(mp_trace, mp_jump):
  
  OPHS, TRS, JMS = path_example.misc_E(mp_trace, mp_jump)

  bc0 = block.from_paths([ OPHS[0], path.rev(OPHS[0]), OPHS[4], path.rev(OPHS[4]), OPHS[2], ])
  bc1 = block.from_paths([ OPHS[1], path.rev(OPHS[1]), ])
  bc2 = block.from_paths([ OPHS[3], ])
  
  BCS = [ bc0, bc1, bc2, ] 
  for ibc in range(len(BCS)): block.set_name(BCS[ibc], "BB%d" % ibc)
  
  return BCS, OPHS, TRS, JMS
  # ----------------------------------------------------------------------

def misc_C(mp_trace):
  
  wd = move_parms.width(mp_trace)
  n_sma = 4; sz_sma = (n_sma-1)*wd
  n_big = 7; sz_big = (n_big-1)*wd

  alt = True     # Serpentine paths 
  mp_jump = None # No jumps needed.
  
  plo0 = ( 1, 1 )
  bc0 = block_example.raster_rectangle(plo0, n_sma, n_big, True, False, alt, mp_trace, mp_jump)
  
  plo1 = rn.add(plo0, (n_sma*wd, 0))
  bc1 = block_example.raster_rectangle(plo1, n_big, n_big, False, True, alt, mp_trace, mp_jump)
  
  plo2 = rn.add(plo1, (n_big*wd, 0))
  bc2 = block_example.raster_rectangle(plo2, n_sma, n_big, True, False, alt, mp_trace, mp_jump)
   
  plo3 = rn.add(plo2, (n_sma*wd, 0))
  bc3 = block_example.raster_rectangle(plo3, n_big, n_big, True, False, alt, mp_trace, mp_jump)
  
  plo4 = rn.add(plo3, (n_big*wd, 0))
  bc4 = block_example.raster_rectangle(plo4, n_sma, n_big, False, True, alt, mp_trace, mp_jump)
  
  plo5 = rn.add(plo1, (0, n_big*wd))
  bc5 = block_example.raster_rectangle(plo5, n_big, n_sma, False, True, alt, mp_trace, mp_jump)
  
  plo6 = rn.add(plo1, (0, -n_sma*wd))
  bc6 = block_example.raster_rectangle(plo6, n_big, n_sma, True, False, alt, mp_trace, mp_jump)
 
  BCS = [ bc0, bc1, bc2, bc3, bc4, bc5, bc6, ] 
  for ibc in range(len(BCS)): block.set_name(BCS[ibc], "BC%d" % ibc)
  
  # Gather all choices in {OPHS}:
  OPHS = []
  for bc in BCS:
    for ich in range(block.nchoices(bc)):
      oph = block.choice(bc, ich)
      OPHS.append(oph)
  
  # # Add links:
  # xdir = (1,0)
  # ydir = (0,1)
  # mp_link = mp_trace
  # OLKS = raster.create_all_raster_raster_links(OPHS, xdir, ydir, mp_link)
  # 
  # # Add contacts:
  # CTS = raster.create_all_raster_raster_contacts(OPHS, xdir, ydir)
  
  return BCS
  # ----------------------------------------------------------------------

def misc_D(mp_trace, contacts):
  
  wd = move_parms.width(mp_trace)

  # Dimensions and counts. The width and height of blocks do not include
  # the overflow of nominal trace sausages.

  nx_toe =  1; szx_toe = nx_toe*wd                 # Gap count and width of serif of base and top blocks.
  nx_sma =  3; szx_sma = (nx_sma-1)*wd             # Trace count and width of roads.
  nx_gap =  2; szx_gap = nx_gap*wd                 # Gap count and width of space between roads.
  nx_big =  2*(nx_toe + (nx_sma-1)) + nx_gap + 1;  # Trace count of base and top blocks.
  szx_big = (nx_big-1)*wd                          # Width of base and top blocks.

  ny_big =  5; szy_big = (ny_big-1)*wd             # Trace count and height of roads 
  ny_sma =  3; szy_sma = (ny_sma-1)*wd             # Trace count and height of base and top blocks.

  alt = True     # Serpentine paths, no jumps.
  mp_jump = None # No jumps needed.
    
  # Bottom wide block:
  plo0 = ( 1, 1 )
  bc0 = block_example.raster_rectangle(plo0, nx_big, ny_sma, True, False, alt, mp_trace, mp_jump)
  
  # Left road block:
  plo1 = rn.add(plo0, ( szx_toe, szy_sma + wd ))
  bc1 = block_example.raster_rectangle(plo1, nx_sma, ny_big, True, False, alt, mp_trace, mp_jump)
  
  # Right road block:
  plo2 = rn.add(plo1, ( szx_sma + szx_gap, 0 ))
  bc2 = block_example.raster_rectangle(plo2, nx_sma, ny_big, True, False, alt, mp_trace, mp_jump)
   
  # Top block:
  plo3 = rn.add(plo0, ( 0, szy_sma + wd + szy_big + wd ))
  bc3 = block_example.raster_rectangle(plo3, nx_big, ny_sma, True, True, alt, mp_trace, mp_jump)

  BCS = [ bc0, bc1, bc2, bc3, ] 
  for ibc in range(len(BCS)): block.set_name(BCS[ibc], "BD%d" % ibc)
  
  CTS = []
  if contacts:
    # Create and attach contacts between the blocks:
    szmin = 0.95*wd
    rszmin = 0.001
    nbc = len(BCS)
    for ibc0 in range(nbc):
      bc0 = BCS[ibc0]
      for ibc1 in range(ibc0):
        bc1 = BCS[ibc1]
        assert bc0 != bc1
        tol = 0.20*wd  # Tolerance for overlaps, tilts, etc.
        CSS01 = contact.from_blocks(bc0, bc1, szmin, rszmin, tol, None)
        CTS += CSS01
  return BCS, CTS
  # ----------------------------------------------------------------------

def misc_E(mp_trace, mp_jump):
  wd = move_parms.width(mp_trace)
  n = 4; sz = (n-1)*wd
  
  plo0 = ( 1, 1 )
  alt0 = True
  bc0 = block_example.raster_rectangle(plo0, n, n, True, False, alt0, mp_trace, mp_jump)
  
  plo1 = rn.add(plo0, (n*wd, 0))
  alt1 = False
  bc1 = block_example.raster_rectangle(plo1, n, n, False, True, alt1, mp_trace, mp_jump)
  
  BCS = [ bc0, bc1, ] 
  for ibc in range(len(BCS)): block.set_name(BCS[ibc], "BE%d" % ibc)
  
  return BCS
  # ----------------------------------------------------------------------
 

def misc_G(mp_cont,mp_fill,mp_jump):

  PHS, TRS02, TRS1 = path_example.misc_G(mp_cont,mp_fill,mp_jump)
  
  bcA = block.from_paths( [ PHS[0], path.rev(PHS[0]), PHS[2], path.rev(PHS[2]), ])
  bcB = block.from_paths( [ PHS[1], path.rev(PHS[1]), ])

  BCS = [ bcA, bcB, ]
  for ibc in range(len(BCS)): block.set_name(BCS[ibc], "BG%d" % ibc)
  
  return BCS, PHS, TRS02, TRS1
  # ----------------------------------------------------------------------

def raster_raster_block_contact(bc0, bc1):
  # The {X} range of the contact will be {[xlo _ xhi]}
  xlo = -inf 
  xhi = +inf
  mv = [None,None]    # Top raster traces at top of {bc0} and bottom of {bc1}.
  for ksd in range(2):
    bck = (bc0, bc1)[ksd][0]
    # Pick some choice of the block, {bc0}  or {bc1}
    ophk = block.choice(bck, 0)
    omvk = path.elem(ophk, (bc0, bc1)[ksd][1])
    
    # Check if they are both traces of the same width
    wdk = move.width(omvk)
    assert wdk > 0 # Must be trace not jump.
    # Save the unoriented move into {mv[ksd]}:
    mv[ksd], drk = move.unpack(omvk)
    # Find the {X} span and intersect with current {[xlo _ xhi]}
    pk = move.pini(mv[ksd])
    qk = move.pfin(mv[ksd])
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
  return ct
  # ----------------------------------------------------------------------

