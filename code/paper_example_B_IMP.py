#! /usr/bin/python3
# Implementation of module {paper_figures_B}

import paper_example_B
import move
import move_example
import move_parms
import contact
import path
import path_example
import rootray_shape
import rootray
import raster
import block
import hacks
import pyx
import rn
import sys
from math import sqrt, hypot, sin, cos, atan2, floor, ceil, inf, nan, pi

######################################################################
### Simple example of contacts being covered by paths:

def make_simple_cover(mp_fill, mp_jump):

  p0 = (  0,  1)
  p1 = (  5,  3)
  p2 = (  6,  1)
  p3 = ( 11,  3)
  p4 = ( 12,  1)
  
  dp23 = rn.sub(p3, p2)
  vp23,mp23 = rn.dir((-dp23[1], +dp23[0]))
  p6 = rn.mix3(1.00,p2, 1.50,dp23, 1.00,vp23)
  p7 = rn.mix3(1.00,p2, 0.50,dp23, 1.00,vp23)
  p8 = rn.add(p7, (-1, +2))
  p5 = rn.add(p6, (+1, -2))

  P = path.from_points(( ( p0, p1, p2, p3, p4, ), ( p5, p6, p7, p8, ) ), mp_fill, mp_jump); path.set_name(P, "P", True)
  

  dp01 = rn.sub(p1, p0)
  vp01,mp01 = rn.dir((-dp01[1], +dp01[0]))
  q1 = rn.mix3(1.00,p0, 0.40,dp01, 1.00,vp01)
  q2 = rn.mix3(1.00,p0, 1.2,dp01, 1.00,vp01)
  q3 = rn.mix(1.00,q2, 1.6,vp01)
  q0 = rn.mix(1.00,q1, 3.0,vp01)
  
  Q = path.from_points(( q0, q1, q2, q3, ), mp_fill, None); path.set_name(Q, "Q", True)
  
  OPHS = [ P, Q, ]
  
  wdf = move_parms.width(mp_fill)
  tol = 0.20*wdf  # Tolerance for overlaps, tilts, etc.

  mvdir = None # Let {from_moves} determine it.
  ct0 = contact.from_moves(path.elem(P,0), path.elem(Q,1), mvdir, 0.05, 0.00, tol) 
  assert ct0 != None 
  contact.set_name(ct0, "C0")
  path.add_contact(P, 0, ct0); contact.add_side_path(ct0, 0, P, None)
  path.add_contact(Q, 1, ct0); contact.add_side_path(ct0, 1, Q, None)

  ct1 = contact.from_moves(path.elem(P,2), path.elem(P,6), mvdir, 0.05, 0.00, tol)
  assert ct1 != None
  contact.set_name(ct1, "C1")
  path.add_contact(P, 0, ct1); contact.add_side_path(ct1, 0, P, None)
  path.add_contact(P, 1, ct1); contact.add_side_path(ct1, 1, P, None)

  CTS = [ ct0, ct1, ] 
  
  return OPHS, CTS
  # ----------------------------------------------------------------------

######################################################################
### Simple example of moves:

def make_simple_moves(mp_fill, mp_jump):
  p = (1,1)
  q = (6,2)
  tr = move.make(p, q, mp_fill); move.set_name(tr, "T")
  jm = move.make(p, q, mp_jump); move.set_name(jm, "J")
  phtr = path.from_moves([tr,]); path.set_name(phtr, "PT", False)
  phjm = path.from_moves([jm,]); path.set_name(phtr, "PJ", False)
  
  OPHS = [ phtr, phjm ]
  return OPHS
  # ----------------------------------------------------------------------

######################################################################
### Simple example of paths

def make_simple_path(mp_fill, mp_jump):
  pd = (9,2)
  pe = (9,0)
  pf = (5,1)
  pg = (1,0)
  ph = (0,2)
  
  mvde = move.make(pd, pe, mp_fill); move.set_name(mvde, "Tde")
  mvef = move.make(pe, pf, mp_fill); move.set_name(mvef, "Tef")
  mvfg = move.make(pf, pg, mp_jump); move.set_name(mvfg, "Jfg")
  mvgh = move.make(pg, ph, mp_fill); move.set_name(mvgh, "Tgh")

  MVS = [ mvde, mvef, mvfg, mvgh ]
  ph = path.from_moves(MVS)
  path.set_name(ph, "P", False)
  OPHS = [ path.rev(ph), ]
  return OPHS
  # ----------------------------------------------------------------------

######################################################################
### Simple example of contacts

def make_simple_contacts(mp_fill):
  
  pa0 = (3, 1)
  pa1 = (7, 1)
  pha = path.from_points((pa0, pa1,), mp_fill, None)
  path.set_name(pha, "PA", True)
  
  pb0 = (2, 2)
  pb1 = (9, 2)
  phb = path.from_points((pb0, pb1,), mp_fill, None)
  path.set_name(phb, "PB", True)
  
  pc0 = (8, 4)
  pc1 = (7, 3)
  pc2 = (5, 3)
  pc3 = (4, 4)
  pc4 = (3, 3)
  pc5 = (1, 3)
  pc6 = (0, 4)
  phc = path.from_points((pc0, pc1, pc2, pc3, pc4, pc5, pc6, ), mp_fill, None)
  path.set_name(phc, "PC", True)
  
  omva0 = path.elem(pha,0)
  omvb0 = path.elem(phb,0)
  omvc1 = path.elem(phc,1)
  omvc4 = path.elem(phc,4)
  
  wdf = move_parms.width(mp_fill)
  tol = 0.20*wdf  # Tolerance for overlaps, tilts, etc.

  mvdir = (1,0)
  cta0b0 = contact.from_moves(omva0, omvb0, mvdir, 0.1, 0.0, tol); contact.set_name(cta0b0, "Ca0b0")
  contact.add_side_path(cta0b0, 0, pha, 0); path.add_contact(pha, 0, cta0b0)
  contact.add_side_path(cta0b0, 1, phb, 0); path.add_contact(phb, 1, cta0b0)

  ctb0c1 = contact.from_moves(omvb0, omvc1, mvdir, 0.1, 0.0, tol); contact.set_name(ctb0c1, "Cb0c1")
  contact.add_side_path(ctb0c1, 0, phb, 0); path.add_contact(phb, 0, ctb0c1)
  contact.add_side_path(ctb0c1, 1, phc, 1); path.add_contact(phc, 1, ctb0c1)

  ctb0c4 = contact.from_moves(omvb0, omvc4, mvdir, 0.1, 0.0, tol); contact.set_name(ctb0c4, "Cb0c4")
  contact.add_side_path(ctb0c4, 0, phb, 0); path.add_contact(phb, 0, ctb0c4)
  contact.add_side_path(ctb0c4, 1, phc, 4); path.add_contact(phc, 1, ctb0c4)

  OPHS = [pha, phb, phc, ]
  CTS = [cta0b0, ctb0c1, ctb0c4, ]
  
  return OPHS, CTS
  # ----------------------------------------------------------------------

######################################################################
### The "turkey" part

def key_coords(dil):
  # The raster line endpoints lie on the boundary of a region {D}.
  # Returns the key coordinates and radii of the {D} region, dilated by {dil}.
  # Assumes that the raster spacing {wdf} is 1.
  
  # Dimensions and coordinates of undilated {D}:
  dxms =  5.00   # X distance between centers of holes 'm' and 's'.
  dxgi =  3.00   # X distance between centers of slot and right shoulders.
  dyrm =  2.00   # Y distance between centers of slot bottom and hole 'm'.
  wdtb =  8.00   # Width of right tab, slot center to shoulder center.
  ygap =  0.10   # Extra gap between top/bot scanline and top/bot edges of {D}

  rado =  4.500 + ygap  # Radius of hole rims and shoulders.
  radi =  2.499         # Radius of holes and half-width of notch and slot.
  radt =  radi  # Radius of fillet between head and neck.

  # Dimensions and parameters of dilated {D}:
  
  ya =  1.00 - ygap - dil  # Y of bottom edge of dilated {D}.
  yz = 18.00 + ygap + dil  # Y of top edge of dilated {D}
  
  ri = radi - dil           # Inner radius of holes. 
  assert ri > 0, "hole radius is too small"

  ro = rado + dil    # Outer radius of hole rims and smooth shoulders.

  rt = radt - dil    # Radius of neck-head fillets
  
  yr = ya + ro              # Y of centers of slot/notch bottom arcs.
  ym = yr + dyrm            # Y of center of hole 'm'.
  ys = yz - ro              # Y of center of hole 's'.
  
  xm = ceil(radi + 3)       # X of center of hole 'm'.
  xs = xm + dxms            # X of center of hole 's'.
  xc = xs + ro              # X of left edge of notch.
  xn = xc + ri              # X of midline of notch.
  
  # Compute the center {(xt,yt)} of the 
  # fillet of radius {rt} that between the 'head' and the 'neck'
  xt, yt = hacks.filet_center((xm,ym), (xn,yr), ro, rt)
  
  xf = xn + ri              # X of right edge of notch.
  xg = xf + ro              # X of midline of slot.
  xx = xg - ri              # X of left edge of slot.
  xz = xg + ri              # X of right edge of slot.
  xi = xg + wdtb            # X of center of right shoulders.
  xj = xi + ro              # X of right edge of part.

  # Sine and co-sine of angles of left perimeter
  dms = rn.dist((xm,ym), (xs,ys))           # Dist between centers of hole 'm' and hole 's'.
  cms = (xs - xm)/dms; sms = (ys - ym)/dms  # Cos and sin of CW tilt of upper straight line from horz.

  dmt = rn.dist((xm,ym), (xt,yt))           # Dist between centers of hole 'm' and filet 't'.
  cmt = (ym - yt)/dmt; smt = (xt - xm)/dmt  # Cos and sin of CCW incl of 'm'-'t' from vert.

  dnt = rn.dist((xn,yr), (xt,yt))           # Dist between centers of filet 't' and bot 'n' of notch.
  cnt = (xn - xt)/dnt; snt = (yr - yt)/dnt  # Cos and sin of CCW tilt of 'n'-'t' from horz.
  
  return xm,xs,xt,xc,xn,xf,xg,xi,xj, ya,yr,ym,ys,yt,yz, ri,ro,rt, cms,sms,cmt,smt,cnt,snt
  # ----------------------------------------------------------------------
  
def trace_circle(xc,yc,rd, a0,a1,step, ylo,yhi):
  # Returns points on the circle with center {(xc,yc)} and radius {rd},
  # sampled between angles {a0} and {a1} in increments of arc length
  # {step}, clippled to the Y range {[ylo _ yhi]}. The tracing is CCW if
  # {step} is positive, CW if negative.
  PTS = []
  na = ceil(rd*(a1-a0)/step)      # Number of steps.
  if na > 0:
    da = (a1 - a0)/na # Angle increment for each step.
    for ia in range(na+1):
      a = a0 + da*ia
      x = xc + rd*cos(a)
      y = yc + rd*sin(a)
      if y >= ylo and y <= yhi:
        PTS.append((x,y))
  return PTS
  # ----------------------------------------------------------------------

def contour_points(dil, step):
  # The ponts are on the boundary of the filling region {D} expanded by {dil}.
  # Points on curve sections will be spaced by {step}.
  # 
  # The result is a list {PTCS} of lists. Each element of {PTCS} is a
  # list of points that describe a connected component of the contour.
  # Each component has the last point repeating the first point.

  xm,xs,xt,xc,xn,xf,xg,xi,xj, ya,yr,ym,ys,yt,yz, ri,ro,rt, cms,sms,cmt,smt,cnt,snt = key_coords(dil)
  
  ams = atan2(sms,cms) # CW tilt angle of top straight seg on left side from horz.
  amt = atan2(smt,cmt) # CCW incl angle of 'm'-'t' line from vert.
  ant = atan2(snt,cnt) # CCW tilt angle of 'n'-'t' line from horz.
  
  # ....................................................................
  # Outer perimeter:
  PTOS = []
  
  yeps = 0.05 # Min Y space between arc points and scan-lines.
  
  # Bottom left edge:
  PTOS.append((xn, ya))
  PTOS.append((xi, ya))
  
  # Bottom right shoulder:
  PTS = trace_circle(xi,yr,ro, -0.5*pi,00.0*pi,+step, ya+yeps, +inf)
  PTOS += PTS

  # Top right shoulder:
  PTS = trace_circle(xi,ys,ro, 00.0*pi,+0.5*pi,+step, -inf, yz-yeps)
  PTOS += PTS

  # Top right edge:
  PTOS.append((xi, yz))
  PTOS.append((xg, yz))
  
  # Top right shoulder of notch: 
  PTS = trace_circle(xg,ys,ro, +0.5*pi,+1.0*pi,+step, -inf, yz-yeps)
  PTOS += PTS

  # Bottom semicircle of notch:
  PTS = trace_circle(xn,yr,ri, 00.0*pi,-1.0*pi,-step, -inf, +inf)
  PTOS += PTS
  
  # Rim of hole 's': 
  PTS = trace_circle(xs,ys,ro, 00.0*pi,(0.5*pi+ams),+step, -inf, +inf)
  PTOS += PTS
  
  # Rim of hole in hole 'm': 
  PTS = trace_circle(xm,ym,ro, (0.5*pi+ams),(1.5*pi+amt),+step, -inf, +inf)
  PTOS += PTS
  
  # Rim of filet 't': 
  PTS = trace_circle(xt,yt,rt, (0.5*pi+amt),(0.0*pi+ant),-step, -inf, +inf)
  PTOS += PTS
  
  # Rim of bootom left shoulder: 
  PTS = trace_circle(xn,yr,ro, (-1.0*pi+ant),-0.5*pi,+step, ya+yeps, +inf)
  PTOS += PTS
  
  # Close curve:
  PTOS.append(PTOS[0])
  # ....................................................................
  
  # ....................................................................
  # Hole in hole col 1:
  PTH1S = trace_circle(xm,ym,ri, +2.0*pi,00.0*pi,-step, -inf, +inf)
  # Fix possible roundoff errors in last point:
  if PTH1S[0] != PTH1S[-1]: PTH1S[-1] = PTH1S[0]
  # ....................................................................
  
  # ....................................................................
  # Top hole in hole col 2:
  PTH2S = trace_circle(xs,ys,ri, +2.0*pi,00.0*pi,-step, -inf, +inf)
  # Fix possible roundoff errors in last point:
  if PTH2S[0] != PTH2S[-1]: PTH2S[-1] = PTH2S[0]
  # ....................................................................

  # ....................................................................
  # Slot:
  PTSS = []

  # Bottom half-circle of slot: 
  PTS = trace_circle(xg,yr,ri, 00.0*pi,-1.0*pi,-step, -inf, +inf)
  PTSS += PTS
  
  # Top semicircle of slot:
  PTS = trace_circle(xg,ys,ri, +1.0*pi,00.0*pi,-step, -inf, +inf)
  PTSS += PTS
  
  PTSS.append(PTSS[0])
  # ....................................................................

  PTCS = [ PTOS, PTH1S, PTH2S, PTSS ]
  
  # Remove close points and rotate 180:
  xctr = (xm - ro + xj)/2
  yctr = (ya + yz)/2
  dmin = 1.0e-3
  for iph in range(len(PTCS)):
    PTS = PTCS[iph]
    PTS = hacks.poly_cleanup(PTS, dmin)
    if PTS[0] != PTS[-1]: PTS.append(PTS[0])
    PTS = [ (2*xctr - p[0], 2*yctr - p[1]) for p in PTS ]
    PTCS[iph] = PTS

  return PTCS
  # ----------------------------------------------------------------------

def raster_endpoints(PTCS, ystep, yphase):
  # Computes the endpoints of the raster lines.
  #
  # The result is a list {PTRS} such that {PTRS[isc][jrs][k]} is endpoint
  # {k} (0=left, 1=right) of raster {jrs} on scan-line {isc} from bottom.
  #
  # The scan-line spacing will be {wdf}, and the Y coords of scan-lines
  # are integer multiples of {wdf}.
  
  xdir = (1,0)
  ydir = (0,1)
  PTRS = raster.endpoints_from_polys(PTCS, xdir,ydir, ystep,yphase)
  return PTRS
  # ----------------------------------------------------------------------
 
def link_points(PTCS, PTRS):
  # Returns the points that comprise the link paths.]
  #
  # The result is a list {PTLS} such that {PTLS[isc][jlk][k]} is the list of points
  # that define the link path of index {jlk} between scan-lines {isc} and {isc+1}.
  #
  # The ponts lie on the boundary of the filling region {D}.
  
  ydir = (0,1)
  PTLS = raster.link_points_from_raster_endpoints(PTCS, PTRS, ydir)
  return PTLS
  # ----------------------------------------------------------------------
 
def contour(mp_cont, mp_fill, mp_jump):
  # Returns the list {OCRS} of contour paths for {paper_figures_B}.  See the comments in that
  # function for the meaning of the parameters.

  wdc = move_parms.width(mp_cont)
  wdf = move_parms.width(mp_fill)
  
  # Dilate {D} so that the contour touches the tips of rasters:
  dil = (wdf + wdc)/2
  step = 0.10*wdf     # Step to trace circular arcs (mm).
  PTCS = contour_points(dil, step)

  OCRS = [] # List of contour paths.
  for icr in range(len(PTCS)):
    PSi = PTCS[icr]
    ph = path.from_points(PSi, mp_cont, None)
    pname = ("C%02d" % icr)
    path.set_name(ph, pname, False)
    for jmv in range(path.nelems(ph)):
      mname = pname + (".%03d" % jmv)
      move.set_name(path.elem(ph, jmv), mname)
    OCRS.append(ph)
  
  path.compute_contour_nesting(OCRS)
  return OCRS
  # ----------------------------------------------------------------------

def filling(mp_fill, mp_jump):
  # Returns the lists {OPHS} and {OLKS} for {make_figure}.  See the comments in that
  # function for the meaning of the parameters.
  
  wdf = move_parms.width(mp_fill)
  
  # Get the contours of the undilated region {D}:
  dil = 0
  step = 0.10*wdf
  PTCS = contour_points(dil, step)
  
  # Compute the raster endpoint by stabbing {D} with horizontal lines:
  ystep = 1.00*wdf
  yphase = 0.00*wdf
  PTRS = raster_endpoints(PTCS, ystep, yphase)

  # Make them into links:
  OPHS = raster.from_endpoints(PTRS, mp_fill)
  
  # Break the boundary of {D} at scan-lines:
  PTLS = link_points(PTCS, PTRS)
  
  # Turn the link point lists into link paths:
  OLKS = []
  for isc in range(len(PTLS)):
    LSi = PTLS[isc]  # Links between scan-line {isc} and scan-line {isc+1}
    for jlk in range(len(LSi)):
      LSij = LSi[jlk]  # Points on link path {jlk} above scan-line {isc}
      LSij = snap_points_to_raster_ends(LSij, PTRS)
      assert len(LSij) > 0
      pname = "L%02d.%d" % (isc,jlk)
      ph_raw = path.from_points(LSij, mp_fill, mp_jump)
      nmv_raw = path.nelems(ph_raw)

      # Simplify the link:
      stol = 0.030*wdf
      ph = path.simplify(ph_raw, stol, mp_jump)
      nmv = path.nelems(ph)
      if nmv != nmv_raw:
        T_raw = path.fabtime(ph_raw)
        T = path.fabtime(ph)
        sys.stderr.write("link reduced from %d moves (%.6f s) to %d  (%.6f s)\n" % (nmv_raw,T_raw,nmv,T))
      elif nmv != 1:
        T = path.fabtime(ph)
        sys.stderr.write("link has %d moves (%.6 s)\n" % (nmv,T))

      path.set_name(ph, pname, False)
      for kmv in range(path.nelems(ph)):
        mname = pname + (".%d" % kmv)
        move.set_name(path.elem(ph, kmv), mname)
      OLKS.append(ph)
      
  # Attach the links to the raster paths, discardng sub-attached ones:
  PLSnew = []
  for lk in OLKS:
    lk = attach_link(lk, OPHS);
    if lk != None:
      PLSnew.append(lk)
  OLKS = PLSnew
  return OPHS, OLKS
  # ----------------------------------------------------------------------
   
def snap_points_to_raster_ends(PS, PTRS):
  # Given a list {PS} of points, and a list {PTRS} of all raster endpoint pairs
  # (as returned by {raster_endpoints}), returns a copy of {PS} where
  # every point that is close to a raster endpoint is replaced by that endpoint.
  # Repated points are removed from the list.
  tol = 1.0e-2
  PSnew = []
  for p in PS:
    for RSi in PTRS: # Raster endpoint pairs in a scan-line.
      for RSij in RSi: # A raster endpoint pair.
        for q in RSij: # A raster endpoint.
          if rn.dist(p, q) < tol: p = q
    if len(PSnew) == 0 or p != PSnew[-1]:
      PSnew.append(p)
  return PSnew
  # ----------------------------------------------------------------------

def attach_link(lk, OPHS):
  # If both ends of the link path {lk} coincide with endpoints of
  # raster paths in {OPHS}, attaches {lk} to those paths, in the proper orientations,
  # and returns {lk}.  Otherwise (if only 0 or 1 of the ends of {lk} is attachable)
  # does nothing and returns {None}
  
  p = [ path.pini(lk), path.pfin(lk) ]
  PAS = [ paths_starting_at(OPHS, p[elk]) for elk in range(2) ]
 
  if len(PAS[0]) == 0 or len(PAS[1]) == 0:
    sys.stderr.write("  discarded not-biconnected link path\n")
    path.show(sys.stderr, "    ", lk, "\n", False, 0,0,0)
    return None
  else:
    # Connect {lk} to the paths {PAS[0],PAS[1]}:
    for elk in range(2):
      # Orient {lk} so that it ends at {p[elk]}:
      olk = path.rev(lk) if elk == 0 else lk
      for ors in PAS[elk]:
        path.add_link(ors, olk)
    return lk
  # ----------------------------------------------------------------------
 
def paths_starting_at(OPHS, p):
  # Returns the paths in the list {OPHS} that have an endpoint 
  # practically coincident to {p}.  The paths are reveresed as
  # needed so that they starts at {p}.
  
  tol = 1.0e-6
  PAS = []
  for rs in OPHS:
    ors = None
    if rn.dist(p, path.pini(rs)) < tol: ors = rs
    if rn.dist(p, path.pfin(rs)) < tol: ors = path.rev(rs)
    if ors != None:
      PAS.append(ors)
  return PAS
  # ----------------------------------------------------------------------
  
def contacts(OPHS):
  # Creates contacts between the paths in {OPHS}
  wdf = move.width(path.elem(OPHS[0],0)) # Expects all to be the same width.
  tol = 0.20*wdf  # Tolerance for overlaps, tilts, etc.

  CTS = contact.from_path_lists(OPHS, OPHS, 0.25, 0.01, tol, None)
  return CTS
  # ----------------------------------------------------------------------

def contact_graph(OPHS, CTS):
  # Computes the vertices {VGS} and its edges {EGS}
  # of the contact graph for the filling, given the list of filling
  # raster paths {OPHS} and the list of contacts {CTS}.
  # 
  # The elements of {VGS} are points.  The elements of {EGS} are
  # pairs of indices into {VGS}.
  
  # The vertices are the midpoints of the raster elements:
  VGS = []
  for oph in OPHS:
    p = path.pini(oph)
    q = path.pfin(oph)
    m = rn.mix(0.5, p, 0.5, q)
    VGS.append(m)
    
  # The edges are the pairs {(iph,jsc)} such that there is contact
  # between {OPHS[iph]} and {OPHS[jsc]}:
  EGS = []
  for cc in CTS:
    mvs = contact.side_moves(cc)
    e = [None,None]
    for iph in range(len(OPHS)):
      ophi = OPHS[iph]
      assert path.nelems(ophi) == 1
      omvi = path.elem(ophi,0)
      mvi, dri = move.unpack(omvi)
      for kep in range(2):
        if mvs[kep] == mvi: e[kep] = iph
    EGS.append(tuple(e))
  
  return VGS,EGS
  # ----------------------------------------------------------------------

def make_turkey(mp_cont, mp_fill, mp_jump):
    
  OCRS = contour(mp_cont, mp_fill, mp_jump)
  OPHS,OLKS = filling(mp_fill, mp_jump)
  CTS = contacts(OPHS) #
  VGS, EGS = contact_graph(OPHS, CTS)
  
  ncr = len(OCRS)
  nph = len(OPHS)
  nlk = len(OLKS)
  nct = len(CTS)
  sys.stderr.write("{paper_example_B}: %d contours, %d rasters, %d links, %d contacts\n" % (ncr,nph,nlk,nct))
  assert len(VGS) == nph
  assert len(EGS) == nct
  
  return OCRS,OPHS,OLKS,CTS,VGS,EGS
  # ----------------------------------------------------------------------
