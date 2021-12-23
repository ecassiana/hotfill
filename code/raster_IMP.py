# Implementation of {raster}.
# Last edited on 2021-11-25 11:33:14 by stolfi

import raster
import path
import contact
import move
import move_parms
import block
import rootray
import rootray_shape
import rootray_cart
import hacks

import rn
import rmxn

import sys
from math import sqrt, sin, cos, atan2, log, exp, floor, ceil, inf, nan, pi

def get_spacing_and_phase(OPHS, xdir, ydir):

  nph = len(OPHS)
  assert nph > 0, "no rasters given"
  
  # Guess the scanline spacing {ystep}:
  totwd = 0
  for oph in OPHS:
    assert path.nelems(oph) == 1
    mv = path.elem(oph, 0)
    assert not move.is_jump(mv)
    totwd += move.width(mv)
  ystep = totwd/nph
  assert ystep > 1.0e-6
  
  # Find the lowest raster in {ydir} direction and its projection {ymin}:
  ymin = +inf
  for oph in OPHS:
    xm, ym = path.mean_projections(oph, xdir, ydir)
    ymin = min(ymin, ym)
    
  if nph == 1:
    # We can use {ymin} as the phase:
    yphase = ymin
  else:
    # Take {ymin} as the guesed phase and refine it so that the mean dev is zero:
    sum_yr = 0
    isc_min = +inf
    for oph in OPHS:
      xm, ym = path.mean_projections(oph, xdir, ydir)
      yd = ym - ymin
      isc = floor(yd/ystep + 0.5)
      yr = yd - ystep*isc
      assert abs(yr) < 0.4*ystep, "raster deviation from scanline too large"
      isc_min = min (isc, isc_min)
      sum_yr += yr
    yphase = ymin + sum_yr/nph + isc_min*ystep

    # Check for quantized spacing:
    for oph in OPHS:
      xm, ym = path.mean_projections(oph, xdir, ydir)
      yd = ym - yphase
      yr = yd - ystep*floor(yd/ystep + 0.5)  # Deviation from ideal position.
      assert abs(yr) < 0.2*ystep

  return ystep, yphase
  # ----------------------------------------------------------------------

def point_scanline_index(p, xdir, ydir, ystep, yphase):
  ym = rn.dot(p, ydir)
  yd = ym - yphase
  isc = int(floor(yd/ystep + 0.5))
  return isc
  # ----------------------------------------------------------------------

def path_scanline_index(oph, xdir, ydir, ystep, yphase):
  xm, ym = path.mean_projections(oph, xdir, ydir)
  yd = ym - yphase
  isc = int(floor(yd/ystep + 0.5))
  return isc
  # ----------------------------------------------------------------------

def sort_by_scanline(OPHS, xdir, ydir, ystep, yphase):

  nph = len(OPHS)
  if nph <= 1: return OPHS.copy()
  
  def scanpos(oph):
    # Returns a pair {(isc, xm)} where {isc} is the index
    # of the scanline nearest to the raster fill element {oph},
    # and {xm} is its mean coordinate along the direction {xdir}.
    xm, ym = path.mean_projections(oph, xdir, ydir)
    isc = int(floor((ym - yphase)/ystep + 0.5))
    return (isc, xm)
    
  OPHS_srt = sorted(OPHS, key = scanpos)
  return OPHS_srt
  # ----------------------------------------------------------------------

def separate_by_scanline(OPHS, xdir, ydir, ystep, yphase):
  
  SCS = []       # List of lists of indices of fillers in each scan-line.
    
  isc_prev = None # Scanline index of {oph_prev}.
  for iph in range(len(OPHS)):
    oph_this = OPHS[iph]
    # Get scanline index:
    isc_this = path_scanline_index(oph_this, xdir, ydir, ystep, yphase)
    assert isc_this >= 0, "negative scanline index -- check {yphase}"
    if isc_this == isc_prev :
      # Looks like they are on the same scanline:
      assert len(SCS) > isc_this
      assert len(SCS[isc_this]) > 0
    elif isc_prev == None or isc_this > isc_prev:
      # Looks like we got to a new scanline. Extend {SCS}:
      while len(SCS) <= isc_this: SCS.append([].copy())
    else:
      assert False, "rasters are not in scanline order"
    SCS[isc_this].append(iph)
    isc_prev = isc_this
  return SCS
  # ----------------------------------------------------------------------

def make_raster_raster_link(oph0, oph1, xdir, dmax, dxmid, mp_link):
  
  omv0 = path.elem(oph0,0)
  omv1 = path.elem(oph1,0)
  wdl = move_parms.width(mp_link)
  p0 = move.pini(omv0); v0 = rn.sub(move.pfin(omv0), p0); dx0 = rn.dot(v0, xdir) 
  p1 = move.pini(omv1); v1 = rn.sub(move.pfin(omv1), p1); dx1 = rn.dot(v1, xdir)

  if dx0*dx1 > 1.0e-8 and rn.dist(p0,p1) <= dmax:
    if dxmid == None:
      lk = path.from_points([p0, p1,], mp_link, None)
    else:
      dxm = -dxmid if dx0 > 0 else +dxmid
      pm = rn.mix3(0.5,p0, 0.5,p1, dxm,xdir)
      lk = path.from_points([p0, pm, p1,], mp_link, None)

    path.add_link(oph0, path.rev(lk))
    path.add_link(oph1, lk)
  else:
    lk = None
  return lk
  # ----------------------------------------------------------------------

def make_raster_raster_contact(oph0,oph1, xdir, szmin, rszmin, tol):
  assert path.nelems(oph0) == 1
  assert path.nelems(oph1) == 1
  mv0, dr0 = move.unpack(path.elem(oph0,0))
  mv1, dr1 = move.unpack(path.elem(oph1,0))
  mvdir = xdir
  ct = contact.from_moves(mv0, mv1, mvdir, szmin, rszmin, tol)
  if ct != None:
    for isd in range(2):
      oph = (oph0,oph1)[isd]
      path.add_contact(oph, isd, ct)
      contact.add_side_path(ct, isd, oph, 0)
  return ct
  # ----------------------------------------------------------------------

def analyze_fabtime(OPHS, ydir, ytol, minlen):
  debug = False
  
  Trast = 0
  Tlink = 0
  Tjump = 0
  Nrast = 0
  Nlink = 0
  Njump = 0
  for oph in OPHS:
    # Accumulate times from the moves of {oph}:
    Trast_ph = 0
    Tlink_ph = 0
    Tjump_ph = 0
    nmv = path.nelems(oph)
    for imv in range(nmv):
      omv = path.elem(oph, imv)
      p = move.pini(omv)
      q = move.pfin(omv)
      Tfab_omv = move.fabtime(omv)
      if move.is_jump(omv):
        udp = move.ud_penalty(omv)
        if imv > 0 and not move.is_jump(path.elem(oph,imv-1)): 
          # Add the trace-jump penalty:
          Tjump_ph += udp
          if debug: sys.stderr.write("up   %12.8f  %12.8f %12.8f  %12.8f\n" % (udp, Trast_ph,Tlink_ph,Tjump_ph))
        # Count jump time:
        Tjump_ph += Tfab_omv
        if debug: sys.stderr.write("jump %12.8f  %12.8f %12.8f  %12.8f\n" % (Tfab_omv, Trast_ph,Tlink_ph,Tjump_ph))
        if imv < nmv-1 and not move.is_jump(path.elem(oph,imv+1)): 
          # Add the jump-trace penalty:
          Tjump_ph += udp
          if debug: sys.stderr.write("dn   %12.8f  %12.8f %12.8f  %12.8f\n" % (udp, Trast_ph,Tlink_ph,Tjump_ph))
        Njump += 1  
      elif abs(rn.dot(p,ydir) - rn.dot(q,ydir)) <= ytol:
        # Assum move is a raster:
        Trast_ph += Tfab_omv
        if debug: sys.stderr.write("rast %12.8f  %12.8f %12.8f  %12.8f\n" % (Tfab_omv, Trast_ph,Tlink_ph,Tjump_ph))
        # Check for short rasters:
        dpq = rn.dist(p,q)
        if dpq < minlen:
          move.show(sys.stderr, "  short move: ", omv, "", 0)
          sys.stderr.write(" length = %.8f\n" % dpq)
        Nrast += 1  
      else:
        # Assume move is a link:
        Tlink_ph += Tfab_omv
        if debug: sys.stderr.write("link %12.8f  %12.8f %12.8f  %12.8f\n" % (Tfab_omv, Trast_ph,Tlink_ph,Tjump_ph))
        Nlink += 1  
    
    # Paranoia:
    if debug:
      Tpath_exp = Trast_ph + Tlink_ph + Tjump_ph
      Tpath_cmp = path.fabtime(oph)
      sys.stderr.write("path %12.8f %12.8f  %12.8f %12.8f  %12.8f\n" % (Tpath_cmp,Tpath_exp,Trast_ph,Tlink_ph,Tjump_ph))
    assert abs(path.fabtime(oph) - (Trast_ph + Tlink_ph + Tjump_ph)) < 1.0e-8
    
    Trast += Trast_ph
    Tlink += Tlink_ph
    Tjump += Tjump_ph

    # Check whether any links associated with {oph} have 
    # horizontal moves that could be confused with rasters
    # when {oph} is incorporated in a full tool-path:
    for olk in path.get_links(oph) + path.get_links(path.rev(oph)):
      nmv = path.nelems(olk)
      for imv_lk in range(nmv):
        omv_lk = path.elem(olk, imv_lk)
        if abs(move.pini(omv_lk)[1] - move.pfin(omv_lk)[1]) <= ytol:
          move.show(sys.stderr, "  horizontal move in link: ", omv_lk, "\n", 0)
  
  #return Trast, Tlink, Tjump
  return Nrast, Trast, Nlink, Tlink, Njump, Tjump
  # ----------------------------------------------------------------------

def create_all_raster_raster_links_and_contacts(OPHS, xdir, ydir, dmax, dxmid, mp_link): 

  debug = False
  LKS = []
  CTS = []
  
  if len(OPHS) > 0:

    # Paranoia:
    mp_fill = move.parameters(path.elem(OPHS[0],0))
    wd_fill = move_parms.width(mp_fill)
    for oph in OPHS:
      assert path.nelems(oph) == 1
      omv = path.elem(oph, 0)
      assert not move.is_jump(omv)
      assert move.parameters(omv) == mp_fill
    wdf = move_parms.width(mp_fill)

    ystep, yphase = get_spacing_and_phase(OPHS, xdir, ydir)
    OPHS = sort_by_scanline(OPHS, xdir, ydir, ystep, yphase)
    SCS = separate_by_scanline(OPHS, xdir, ydir, ystep, yphase)
    nsc = len(SCS) # Number of scanlines.

    for isc in range(nsc-1):
      # Get the raster lists of two conseutive scanlines:
      SC0 = SCS[isc];   
      SC1 = SCS[isc+1]; 

      if len(SC0) > 0 and len(SC1) > 0:
        for irs0 in range(len(SC0)):
          oph0 = OPHS[SC0[irs0]]
          nars0 = ("%d" % isc) if len(SC0) == 1 else ("%d.%d" % (isc,irs0))
          if debug: path.show(sys.stderr, "    oph0 = ", oph0, "\n", True, 0,0,0);
          for irs1 in range(len(SC1)):
            oph1 = OPHS[SC1[irs1]]
            nars1 = ("%d" % (isc+1)) if len(SC0) == 1 else ("%d.%d" % (isc+1,irs1))
            if debug: path.show(sys.stderr, "    oph1 = ", oph1, "\n", True, 0,0,0)
          
            # Create link paths:
            for kdr0 in range(2):
              oph0r = path.spin(oph0,kdr0)
              nadr0 = "az"[kdr0]
              for kdr1 in range(2):
                oph1r = path.spin(oph1,kdr1)
                nadr1 = "az"[kdr1]
                lk = make_raster_raster_link(oph0r, oph1r, xdir, dmax, dxmid, mp_link)
                if lk != None:
                  LKS.append(lk)
                  name = "L[%s%s:%s%s]" % (nars0,nadr0,nars1,nadr1)
                  path.set_name(lk, name, True)
                  if debug: path.show(sys.stderr, "  created link ", lk, "\n", True, 0,0,0)
            
            # Create the contacts:
            szmin = 1.1 * wd_fill
            rszmin = 0.01
            tol = 0.20*wd_fill
            ct = make_raster_raster_contact(oph0,oph1, xdir, szmin, rszmin, tol)
            if ct != None:
              name = "C[%s:%s]" % (nars0,nars1)
              contact.set_name(ct, name)
              if debug: contact.show(sys.stderr, "  created contact", ct, "\n", 0)
              CTS.append(ct)

  return LKS, CTS
  # ----------------------------------------------------------------------

def endpoints_from_polys(PTCS, xdir,ydir, ystep,yphase):
  # Get the bounding box of the points, expanded:
  
  B = None
  for PTS in PTCS:
    for p in PTS:
      B = rn.box_include_point(B, p)
  mrg = (2,2) # Extra box margin.
  B = rn.box_expand(B, mrg, mrg) 
  
  eps = 1.0e-5 # Fudge factor to compensate roundoff errors.

  # Create a CSG tree for the region {D}, assuming that the outermost 
  # polygon is CCW (an island):
  SH = rootray_shape.from_polygons(PTCS, +1)
  
  # Obtain the rasters by stabbing it:
  PTRS = endpoints_from_shape(SH, B, xdir,ydir, ystep,yphase)
  return PTRS
  # ----------------------------------------------------------------------

def endpoints_from_shape(S, B, xdir,ydir, ystep,yphase):

  debug = False
  
  # Compute approx range of X and of scan-line indices (possibly too large):
  Adir = ((xdir[0], ydir[0]), (xdir[1], ydir[1])) # Matrix that maps (xdir,ydir) to standard basis.
  Brot = rmxn.box_affine(B, Adir, (0,0))
  xlo = Brot[0][0] - 10;       xhi = Brot[1][0] + 10;
  ylo = Brot[0][1] - 2*ystep;  yhi = Brot[1][1] + 2*ystep;

  isc_lo = int(floor((ylo - yphase)/ystep))
  isc_hi = int(ceil((yhi - yphase)/ystep))
  nsc = isc_hi - isc_lo + 1       # Number of scan-lines in filling.
 
  # Clip the shape {S} to the box {B}:
  K = rootray_shape.box(B[0], B[1])
  S = rootray_shape.intersection([S, K])
  
  PTRS = []
  nsg = 0 # Number of line segments (endpoint pairs).

  sys.stderr.write("filling has %d scan-lines {%d..%d}\n" % (nsc, isc_lo, isc_hi))
  
  for di in range(nsc):
    isc = isc_lo + di
    y = isc*ystep + yphase  # Ordinate of axis of scan-line.
    # Map to {xdir,ydir} system and compute scan-line ref points {p,q}:
    p = rn.mix(xlo, xdir, y, ydir)
    q = rn.mix(xhi, xdir, y, ydir)
    # sys.stderr.write("@@ scanline %d y = %d p = %s q = %s\n" % (isc, y, str(p), str(q)))
    # Intersect scanl-line axis with shape {S}:
    sgn1, TS = rootray_shape.trace(p, q, S)
    assert sgn1 == +1
    nt = len(TS)
    assert nt % 2 == 0, "odd ray-polygon intersection count"
    # Collect the pairs of endpoints:
    PTRSi = [] # List of endpoint pairs on scanline {ìsc}
    for jrs in range(nt//2):
      t0 = TS[2*jrs]
      t1 = TS[2*jrs + 1]
      # sys.stderr.write("  raster %d t = [%12.8f _ %12.8f]\n" % (jrs, t0, t1))
      pr = rn.mix(1-t0, p, t0, q)
      qr = rn.mix(1-t1, p, t1, q)
      PTRSij = (pr, qr)
      PTRSi.append(PTRSij)
      nsg += 1
    PTRS.append(PTRSi)
    
  # Trim empty scan-lines from each end:
  while len(PTRS) > 0 and len(PTRS[0]) == 0: PTRS.pop(0)
  while len(PTRS) > 0 and len(PTRS[-1]) == 0: PTRS.pop()
  if debug: sys.stderr.write("{endpoints_from_shape}: found %d non-empty scanlines and %d segments\n" % (len(PTRS), nsg))
  return PTRS
  # ----------------------------------------------------------------------
  
def from_endpoints(PTRS, mp):

  debug = False
  
  OPHS = []
  for isc in range(len(PTRS)):
    PTRSi = PTRS[isc] # rasters on scan-line {isc}
    for jrs in range(len(PTRSi)):
      PTSRij = PTRSi[jrs] # Endpoints of raster {jrs} on scan-line {isc}.
      assert len(PTSRij) == 2
      p = PTSRij[0]
      q = PTSRij[1]
      mv = move.make(p, q, mp)
      mvname = "T%03d.%02d" % (isc, jrs)
      move.set_name(mv, mvname)
      ph = path.from_moves([mv,])
      phname = "R%03d.%02d" % (isc, jrs)
      path.set_name(ph, phname, False)
      OPHS.append(ph)

  if debug: sys.stderr.write("{from_endpoints}: built %d raster paths\n" % (len(OPHS)))
  return OPHS
  # ----------------------------------------------------------------------

def link_points_from_raster_endpoints(PTCS, PTRS, ydir):
  # Slice the perimeter of {D} by the bands between successive scanlines:
  
  debug = False
  
  PTLS = []
  nsc = len(PTRS) # Number of scan-lines.
  nlk = 0         # Number of link paths
  nbd = 0         # Number of inter-line bands with some links.
  for isc0 in range(nsc-1):
    isc1 = isc0 + 1
    # Create list {LS} of link paths in the strip between scan-lines {isc0} and {isc1}
    LS = []
    RS0 = PTRS[isc0] # List of raster endpoint pairs on scan-line {isc0}
    RS1 = PTRS[isc1] # List of raster endpoint pairs on scan-line {isc1}

    if len(RS0) != 0 and len(RS1) != 0: 
      # Get Y coordinates of the two scan-lines:
      y0 = rn.dot(RS0[0][0], ydir)
      y1 = rn.dot(RS1[0][0], ydir)
      for CS in PTCS: 
        LSj = hacks.cut_poly_with_band(CS, ydir, y0, y1)
        LS += LSj
    nlki = len(LS)
    nlk += nlki
    if nlki > 0: nbd += 1
    PTLS.append(LS)

  if debug: sys.stderr.write("{link_points_from_raster_endpoints}: processed %d scanlines got %d link paths in %d bands\n" % (nsc,nlk,nbd))
  return PTLS
  # ----------------------------------------------------------------------
