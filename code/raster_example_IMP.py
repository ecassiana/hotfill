# Implementation of module {raster_example}

import raster_example
import move
import move_example
import move_parms
import contact
import path
import path_example
import raster
import hacks
import rn
import sys
from math import sqrt, hypot, sin, cos, floor, ceil, inf, nan, pi

def patch_array(cols, rows, nx,ny, islands, mp_cont, mp_fill, mp_link):

  colper = 4 if islands else 3  # Periodicity of columns.
  
  wdf = move_parms.width(mp_fill)
  
  assert type(rows) is int and rows >= 1
  assert type(cols) is int and cols % colper == 1
  assert type(nx) is int and nx >= 1
  assert type(ny) is int and ny >= 5
  
  org = (2*wdf, 2*wdf)  # Coordinate of bottom corner of filling.
  xdir = (1, 0); ydir = (0, 1)          # Axis aligned rasters.
  
  # Dimensions in multiples of {wdf} assuming traces have width 0:
  nx_gp = 3 if islands else nx + 2 # Width of gap between roads.
  ixstep = nx + nx_gp              # Scan-column step between array columns.
  iystep = ny + 1                  # Scan-line step between array rows.
  nx_tt = colper*ixstep            # Period of columns.
  nxtot = (cols-1)*ixstep + nx     # Total width of everything.
  nytot = rows*iystep              # Total height of everything.
  
  def make_raster(iy, row, col, irp, igr):
    # If {col} and {irp} are not {None}, creates a raster that is part
    # of the patch on row {row} and column {col} of the array;
    # specifically, the raster on patch scanline {irp} (in {0..ny-1}).
    # If {col} and {irp} are none, creates instead a beam spanning the
    # whole array, assumed to be at the bottom of row {row}. 
    #
    # In either case, the raster will be on global scanlne {iy} (in
    # {0..nytot}) and will be assigned to group {igr}. The contact lists
    # and link lists will be cleared.
    assert (col == None) == (irp == None)
    assert 0 <= row and row <= rows
    if col != None:
      ix0 = col*ixstep
      ix1 = ix0 + nx  
      name = "%d.%d.%d" % (row,col,irp)
    else:
      ix0 = 0
      ix1 = nxtot
      name = "%d" % row
    p0 = rn.mix(1, org, wdf, (ix0, iy))
    p1 = rn.mix(1, org, wdf, (ix1, iy))
    rmv = move.make(p0, p1, mp_fill)
    move.set_name(rmv, "T" + name)
    rph = path.from_moves((rmv,))
    path.set_name(rph, "P" + name, False)
    path.set_group(rph, igr)
    return rph
    # ......................................................................
  
  OPHS = []
  CTS = [] 
  rph_prev = [None]*cols # Previous raster in each column, possibly a beam.
  
  def make_bottom_beam(row):
    # Creates the beam at the bottom of array row {row} (in {0..cols}).
    # Adds the beam to {OPHS} and saves it into all entries
    # of {rph_prev}.  Adds simple link paths between the beam and the
    # patches below it.
    iy = row*iystep
    bph = make_raster(iy, row, None, None, 0)  # The beam at the bottom of the current row.
    
    # Create links between {bph} to rasters below it:
    for iend in 0, 1:
      col_below = (0, cols-1)[iend] # Coumn index of patch below end {iend}
      rph_below = rph_prev[col_below]
      if rph_below != None:
        name = "%d.%d.%d" % (iend,col_below,iy)
    # Add to list:
    OPHS.append(bph)

    # Reset {rph_prev}
    for col in range(cols):
      rph_prev[col] = bph
    return
    # ......................................................................

  def make_patches(row):
    # Creates the rasters of all patches of row {row} of the array.
    # Assumes that all entries {rph_prev[0..cols-1]]} are set 
    # to the beam raster just below that row. 
    # Adds the patch rasters to {OPHS} and updates {rph_prev}
    # to be the rasters at the top of each patch, or {None}.
    
    # Enumberate patches of this row:
    for col in range(cols):
      # Enumerate the scan-lines of this array slot:
      for irp in range(ny):
        iy = row*iystep + irp + 1   # Scanline index 
        # Decide if the patch on column {col} extends into this scanline:
        rc = col % colper
        island_patch = islands and rc == 2
        skip_bot = (island_patch or rc == colper-1) and (irp == 0 or irp == 1)
        skip_top = (rc == 1 or island_patch) and (irp == ny-2 or irp == ny-1)
        skip = skip_bot or skip_top
        if not skip:
          # Determine the group index {igr} of this patch:
          igr = len(OPHS)
          # Create the raster path with index {irp} in this patch:
          rph_this = make_raster(iy, row, col, irp, igr)
          OPHS.append(rph_this)
          # Remember this patch raster:
          rph_prev[col] = rph_this
        else:
          # Skip raster:
          rph_prev[col] = None
    return
    # ----------------------------------------------------------------------
   
  # Create the rasters:
  for row in range(rows):
    # Beam at bottom of row:
    make_bottom_beam(row)
    make_patches(row)
  make_bottom_beam(rows)
  
  # Add links and contacts:
  dmax = 3.0*wdf
  dxmid = None # Use straight single-trace links.
  LKS, CTS = raster.create_all_raster_raster_links_and_contacts(OPHS, xdir, ydir, dmax, None, mp_link)
  
  # Create the contours:
  OCRS = []
  
  def contour_point(ix, iy, xoff, yoff):
    # Creates a point at scan-grid point {(ix,xy)} wth offset {(xoff,yoff)}
    return rn.mix3(1, org, wdf, (ix,iy), 1, (xoff,yoff))
    # ......................................................................
  
  def hole_contour(ix0,iy0,off):
    # Returns a contour path for the hole in a patch of the array.
    iy1 = iy0 + 3
    iy3 = iy0 + ny + 1
    iy2 = iy3 - 3

    ix1 = ix0 + nx_gp
    ix2 = ix0 + ixstep
    ix6 = ix0 + nx_tt - ixstep
    ix7 = ix1 + nx_tt - ixstep
    ix5 = ix7 - ixstep

    # Hole contour:
    h00 = contour_point(ix0, iy0, +off, +off)
    h01 = contour_point(ix1, iy0, -off, +off)
    h02 = contour_point(ix1, iy2, -off, +off)
    h03 = contour_point(ix2, iy2, +off, +off)
    h04 = contour_point(ix2, iy0, +off, +off)
    h05 = contour_point(ix7, iy0, -off, +off)
    h06 = contour_point(ix7, iy3, -off, -off)
    h07 = contour_point(ix6, iy3, +off, -off)
    h08 = contour_point(ix6, iy1, +off, -off)
    h09 = contour_point(ix5, iy1, -off, -off)
    h10 = contour_point(ix5, iy3, -off, -off)
    h11 = contour_point(ix0, iy3, +off, -off)
    ocrh = path.from_points((h00,h01,h02,h03,h04,h05,h06,h07,h08,h09,h10,h11,h00,), mp_cont, None)
    return ocrh
    # ....................................................................


  def island_contour(ix0,iy0,off):
    ix1 = ix0 + nx_gp
    ix2 = ix0 + ixstep
    ix3 = ix1 + ixstep
    ix4 = ix2 + ixstep

    iy1 = iy0 + 3
    iy3 = iy0 + ny + 1
    iy2 = iy3 - 3

    d0 = contour_point(ix3, iy1, -off, -off)
    d1 = contour_point(ix4, iy1, +off, -off)
    d2 = contour_point(ix4, iy2, +off, +off)
    d3 = contour_point(ix3, iy2, -off, +off)

    ocrq = path.from_points((d0,d1,d2,d3,d0,), mp_cont, None)
    path.set_name(ocrq, "CRd.%d.%d" % (row,4*cg+2), False)
    return ocrq
    # ....................................................................

  # Outer contour:
  wdc = move_parms.width(mp_cont)
  off = (wdc + wdf)/2  # Offset of contour from filling. 
  p0 = contour_point(    0,     0, -off, -off)
  p1 = contour_point(nxtot,     0, +off, -off)
  p2 = contour_point(nxtot, nytot, +off, +off)
  p3 = contour_point(    0, nytot, -off, +off)
  ocro = path.from_points((p0, p1, p2, p3, p0,), mp_cont, None)
  path.set_name(ocro, "CRo", False)
  OCRS.append(ocro)
  
  # Contours of holes and patches:
  for row in range(rows):
    for cg in range(cols//colper):
      # Contours of holes and islands in row {row}, between patches in columns {4*cg} and {4*cg+4}:
      ix0 = colper*cg*ixstep + nx
      iy0 = row*iystep
      ocrh = hole_contour(ix0,iy0,off)
      
      path.set_name(ocrh, "CRh.%d.%d" % (row,4*cg), False)
      OCRS.append(ocrh)
      
      if islands:
        # Island contour:
        ocrq = island_contour(ix0,iy0,off)
        OCRS.append(ocrq)
  
  path.compute_contour_nesting(OCRS)

  return OCRS, OPHS, LKS, CTS
  # ----------------------------------------------------------------------

# MISCELLANEOUS RASTER FILLS FOR TESTING

def rasters_A(mp_fill, xdir, ydir, yphase, eps):

  # Endpoints in the standard coord system, unit step, zero phase, no perturbation:

  p00 = (0, 0); q00 = (3, 0)
  p20 = (5, 0); q20 = (9, 0)
  
  p01 = (0, 1); q01 = (9, 1)
  
  p02 = (1, 2); q02 = (2, 2)
  p22 = (4, 2); q22 = (6, 2)
  p42 = (8, 2); q42 = (9, 2)

  p03 = (1, 3); q03 = (2, 3)
  p23 = (4, 3); q23 = (6, 3)
  
  p = [ 
    [ p00, p20, p01, p02, p22, p42, p03, p23, ],
    [ q00, q20, q01, q02, q22, q42, q03, q23, ],
  ]
  
  # Scale, add phase, and rotate the points to the {xdir,ydir} system:
  wd = move_parms.width(mp_fill)
  OPHS = []
  for k in range(len(p[0])):
    prk = [None, None] # Scaled, rotated, perturbed endpoints
    for iend in range(2):
      pke = rn.scale(wd, p[iend][k])
      ke = 2*k + iend
      dpke = ( eps*sin(17*ke + 3), yphase + eps*cos(31*ke +5) ) # Perturbation to apply.
      ppke = rn.add(pke, dpke)
      prk[iend] = rn.mix(ppke[0], xdir, ppke[1], ydir)
    mvk = move.make(prk[0], prk[1], mp_fill);
    phk = path.from_moves((mvk,))
    OPHS.append(phk)

  return OPHS
  # ----------------------------------------------------------------------

def rasters_B(nph, xdir, mp_trace, mp_link):

  wdt = move_parms.width(mp_trace)
  ydir = (-xdir[1], +xdir[0])

  sz_tr = 3*wdt # Length of traces.

  szx = sz_tr        # Filled region width in {xdir,ydir} system.
  szy = (nph-2)*wdt  # Filled region height in {xdir,ydir} system.

  shift = (1 + min(0, szx*xdir[0], szy*ydir[0]), 1 + min(0, szx*xdir[1], szy*ydir[1]))
  
  dmax = 3.0*wdt
  dxmid = 0.5*wdt

  # Creates {nph+1} traces {TRS[0..nph+1]}:
  TRS = [] # The traces
  for ksc in range(nph+2):
    pk = rn.add(rn.mix(0,    xdir, ksc*wdt,ydir), shift);
    qk = rn.add(rn.mix(sz_tr,xdir, ksc*wdt,ydir), shift);
    TRSk = move.make(pk, qk, mp_trace)
    TRS.append(TRSk)

  # Creates {nph} raster paths with all the traces:
  PHS = []
  for ksc in range(nph+2):
    ph = path.from_moves([TRS[ksc],])
    path.set_name(ph, ("P%d" % ksc), True)
    path.validate(ph)
    PHS.append(ph)
      
  # Create links and contacts for all the traces:
  LKS, CTS = raster.create_all_raster_raster_links_and_contacts(PHS, xdir,ydir, dmax, dxmid, mp_link)
  assert len(LKS) == 2*(nph+1)
  assert len(CTS) == nph+1

  # Remove the first and last paths from {contact.side_paths}:
  contact.clear_side_paths(CTS[0],0)
  contact.clear_side_paths(CTS[nph],1)
  
  # Return only the middle {nph} paths:
  PHS = PHS[1:1+nph]

  return TRS, PHS, LKS, CTS
  # ----------------------------------------------------------------------

