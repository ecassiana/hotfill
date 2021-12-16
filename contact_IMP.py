# Implementation of module {contact}
# Last edited on 2021-11-06 09:23:01 by stolfi

import contact
import move
import move_parms
import path
import block
import hacks; from hacks import fbomb
import rn
import pyx 
from math import nan, inf, sqrt
import sys

class Contact_IMP:
  # The endpoints of the contact are {ct.pt[0]} and {ct.pt[1]}, in 
  # arbitrary order.
  #
  # The contact is between traces {ctside[0]} and {ctside[1]}, 
  # two distinct unoriented {Move} objects that are traces (not jumps).
  #
  # The field {ct.tcov} is a pair (2-tuple) of times, such
  # that {ct.tcov[isd]} is the time for the nozzle to go past the midpont of
  # the contact when tracing the move {ctside[isd]} in its native direction.

  def __init__(self, p0, p1, mv0, tcov0, mv1, tcov1):
    # Read-only fields:
    self.pts = (p0, p1)               # Endpoints.
    self.side = (mv0, mv1)            # {Move} objects on each side.
    self.tclim = +inf                 # Max allowed cooling time.
    self.tcov = (tcov0, tcov1)        # {otcov[isd]] is fabtime from start of side {isd} to contact midpoint.
    self.name = None                  # Name of contact, for debugging.
    # Mutable fields:
    self.sideix = [None, None]        # Arbitrary integers associated with the sides.
    self.sdpaths = [set(), set()]     # {side_paths[isd]} is a set of triples {(ph,dr,imv)} for the paths that contain {side[isd]}.
    self.blocks = [None, None]        # Two {Block} objects that contain the sides of this contact.
    
def make(p0, p1, omv0, omv1):
  assert hacks.is_point(p0); p0 = tuple(p0) # Make sure it is immutable.
  assert hacks.is_point(p1); p1 = tuple(p1) # Make sure it is immutable.
  mv0,dr0 = move.unpack(omv0); assert not move.is_jump(mv0)
  mv1,dr1 = move.unpack(omv1); assert not move.is_jump(mv1)
  assert mv0 != mv1
  m = rn.mix(0.5, p0, 0.5, p1)
  tcov0 = move.cover_time(mv0, m)
  tcov1 = move.cover_time(mv1, m)
  ct = contact.Contact(p0, p1, mv0, tcov0, mv1, tcov1)
  
  # Set a default name:
  if move.has_name(mv0) and move.has_name(mv1):
    na0 = move.get_name(mv0)
    na1 = move.get_name(mv1)
    set_name(ct, "C(%s:%s)" % (na0,na1))
  return ct

def endpoints(ct):
  assert isinstance(ct, contact.Contact)
  return ct.pts

def tcool_limit(ct):
  return ct.tclim

def set_tcool_limit(ct, tclim):
  assert type(tclim) is float or type(tclim) is int, "invalid {tclim}"
  tclim = float(tclim)
  assert tclim >= 0.001, "{tclim} too small"
  ct.tclim = tclim

def pmid(ct):
  return rn.mix(0.5, ct.pts[0], 0.5, ct.pts[1])
  
def side_move(ct, isd):
  return ct.side[isd]
  
def side_moves(ct):
  return ct.side

def side_tcov(ct, isd):
  return ct.tcov[isd]

def which_side(omv, ct):
  mv, dr = move.unpack(omv)
  for isd in range(2):
    if ct.side[isd] == mv:
      return isd
  return None 

def from_moves(omv0, omv1, mvdir, szmin, rszmin, tol):
  if move.is_jump(omv0) or move.is_jump(omv1): return None
  mv0,dr0 = move.unpack(omv0)
  mv1,dr1 = move.unpack(omv1)
  if mv0 == mv1: return None
  p0, p1 = move.shared_border(omv0, omv1, mvdir, tol)
  if p0 == None: return None
  assert p1 != None
  dp = rn.dist(p0, p1)
  if dp == 0: return None
  if dp < szmin: return None
  L0 = move.length(omv0)
  L1 = move.length(omv1)
  if rszmin > 0 and dp < rszmin*min(L0,L1): return None
  ct = make(p0, p1, omv0, omv1)
  return ct
  # ----------------------------------------------------------------------

def from_move_lists(OMVS0, OMVS1, szmin, rszmin, tol, ydir):
  assert type(OMVS0) is tuple or type(OMVS0) is list
  assert type(OMVS1) is tuple or type(OMVS1) is list
  sys.stderr.write("{from_move_lists} OMVS0 = %d OMVS1 = %d\n" % (len(OMVS0),len(OMVS1)))
  
  MVPAIRS = set() # Sets of distinct sorted move pairs from the two lists
  CTS = [ ]
  for imv0 in range(len(OMVS0)):
    for imv1 in range(len(OMVS1)):
      omv0 = OMVS0[imv0]
      omv1 = OMVS1[imv1]
      omv0, omv1 = move.sort_by_midpoint(omv0, omv1, ydir)
      mv0, dr0 = move.unpack(omv0)
      mv1, dr1 = move.unpack(omv1)
      mvpair = (mv0,mv1)
      if mv0 != mv1 and not (mvpair in MVPAIRS):
        mvdir = None # Let the procedure chose a {mvdir} aligned with the moves themselves.
        ct = from_moves(mv0, mv1, mvdir, szmin, rszmin, tol) # Already sorted by {ydir}
        if ct != None: 
          CTS.append(ct)
        MVPAIRS.add(mvpair)
        
  # Set their names, if not defined:
  for ict in range(len(CTS)): 
    ct = CTS[ict]
    if not has_name(ct):
      set_name(ct, "C%d" % ict)

  return CTS
  # ----------------------------------------------------------------------

def from_paths(oph0, oph1, szmin, rszmin, tol, ydir):
  CTS = from_path_lists([oph0,], [oph1,], szmin, rszmin, tol, ydir)
  return CTS
  # ----------------------------------------------------------------------

def from_path_lists(OPHS0, OPHS1, szmin, rszmin, tol, ydir):
  # Collect the two sets of {Move} objects:
  OMVS = [set(), set()]
  for iph in range(2):
    for oph in (OPHS0,OPHS1)[iph]:
      for kmv in range(path.nelems(oph)):
        mvk, drk = move.unpack(path.elem(oph,kmv))
        OMVS[iph].add(mvk)
    # Convert to list in a fixed order:
    OMVS[iph] = list(OMVS[iph])
    OMVS[iph].sort(key = move.pini)

  # create the contacts:
  CTS = from_move_lists(OMVS[0], OMVS[1], szmin, rszmin, tol, ydir)

  # Add the contact-path information:
  for ct in CTS:
    for oph in OPHS0 + OPHS1:
      for imv in range(path.nelems(oph)):
        omv = path.elem(oph, imv)
        mv, dr = move.unpack(omv)
        for isd in range(2):
          if side_move(ct, isd) == mv:
             path.add_contact(oph, isd, ct)    # It is a set, so no problem if repeated {ct}.
             add_side_path(ct, isd, oph, imv)  # It is a set, so no problem if repeated {Path} object.
  return CTS
  # ----------------------------------------------------------------------

def from_blocks(bc0, bc1, szmin, rszmin, tol, ydir):
  # Collect the choices:
  OPHS = [None,None]
  for ibs in range(2):
    bci = (bc0, bc1)[ibs]
    nch = block.nchoices(bci)
    OPHS[ibs] = [ block.choice(bci,ich) for ich in range(nch) ]
  # Now get the contacts:
  CTS = from_path_lists(OPHS[0],OPHS[1], szmin,rszmin, tol, ydir)
  return CTS
  # ----------------------------------------------------------------------

def bbox(CTS):
  B = None
  for ct in CTS:
    B = rn.box_include_point(B, ct.pts[0])
    B = rn.box_include_point(B, ct.pts[1])
  return B
  # ----------------------------------------------------------------------
  
def path_ixcovs(oph, ct):
  ixs = [None, None]
  n = path.nelems(oph)
  for imv in range(n):
    omv = path.elem(oph, imv)
    mv, dr = move.unpack(omv)
    for isd in range(2):
      if side_move(ct, isd) == mv:
        assert ixs[isd] == None, "repeated move in path"
        ixs[isd] = imv
  return tuple(ixs)
  # ----------------------------------------------------------------------

def is_relevant(ct, BCS, ich):
  assert ich == None or type(ich) is int
  assert type(BCS) is list or type(BCS) is tuple
  assert isinstance(ct, contact.Contact)
  for bc in BCS:
    nch = nchoices(bc)
    if ich == None or ich < nch:
      choices = (ich,) if ich != None else range(nch)
      for jch in choices:
        oph = choice(bc, jch)
        ixs = path_ixcovs(oph, ct)
        if ixs[0] != None or ixs[1] != None:
          return True
  return False
  # ----------------------------------------------------------------------

def path_tcov(oph, imv, ct, isd):
  assert isinstance(ct, contact.Contact)
  ph, dr_ph = path.unpack(oph) # For the typechecking.
  if imv == None:
    tc = None
  else:
    mv_ct = ct.side[isd]
    omv = path.elem(oph, imv)
    mv_ph, dr_mv_ph = move.unpack(omv)
    assert mv_ph == mv_ct
    if dr_mv_ph == 0:
      tc =  path.tini(oph, imv) + ct.tcov[isd]
    else:
      tc = path.tfin(oph, imv) - ct.tcov[isd]
  return tc

def path_tcovs(oph, ct):
  ixs = path_ixcovs(oph, ct)
  assert len(ixs) == 2
  tcs = [ None, None ]
  for isd in range(2):
    imv = ixs[isd]
    if imv != None:
      omv = path.elem(oph, imv)
      mv, dr = move.unpack(omv)
      if dr == 0:
        tcs[isd] = path.tini(oph, imv) + ct.tcov[isd]
      else:
        tcs[isd] = path.tfin(oph, imv) - ct.tcov[isd]
  return tuple(tcs)
  # ----------------------------------------------------------------------

# PATHS ASSOCIATED TO CONTACT SIDE

def clear_side_paths(ct, isd):
  assert isinstance(ct, contact.Contact)
  ct.sdpaths[isd] = set()
  # ----------------------------------------------------------------------

def add_side_path(ct, isd, oph, imv):
  assert isinstance(ct, contact.Contact)
  mv = ct.side[isd]
  if imv == None:
    imv = path.find_move(oph, mv)
    assert imv != None, "path does not cover side {isd} of {ct}"
  else:
    omv_c = path.elem(oph, imv)
    mv_c, dr_c = move.unpack(omv_c)
    assert mv_c == mv, "contact side is not the specified move"
  ph, dr = path.unpack(oph)
  ct.sdpaths[isd].add((ph,dr,imv))
  # ----------------------------------------------------------------------

def get_side_paths(ct, isd):
  assert isinstance(ct, contact.Contact)
  return frozenset(ct.sdpaths[isd])
  # ----------------------------------------------------------------------

def check_side_paths(CTS,OPHS,ydir):
  
  def showbug(oph, ct, msg):
    path.show(sys.stderr, "    oph = ", oph, "\n", True, 0,0,0);
    sys.stderr.write("ct =\n");
    show(sys.stderr, "    ct =  ", ct, "\n", 0);
    assert False, msg
    # ....................................................................
  
  # Check from {get.side_paths}
  PHS_org = set((path.unpack(oph))[0] for oph in OPHS)
  for ct in CTS:
    for isd in range(2):
      mvi = side_move(ct, isd)
      for phj, drj, imvj in get_side_paths(ct, isd):
        if not (phj in PHS_org): 
          showbug(phj, ct, "contact refer to path not in {OPHS}")
        if not ct in path.get_contacts(phj,isd):
          showbug(phj, ct, "contact not in contacts of its attached path")
        ophj = path.spin(phj, drj)
        mvk, drk = move.unpack(path.elem(ophj, imvj))
        if mvk != mvi:
          showbug(ophj, ct, "atached path/move do not match contact's side")

  # Check from {path.get_contacts}:
  CTS_org = set(CTS)
  for oph in OPHS:
    ph, dr = path.unpack(oph)
    for isd in range(2):
      for ct in path.get_contacts(oph, isd):
        if not ct in CTS_org:
          showbug(oph, ct, "path refers to contact not in {CTS}")
        mvi = side_move(ct, isd)
        imv = path.find_move(oph, mvi)
        if imv == None:
          showbug(oph, ct, "side of attached contact does not occur in path")
        if not (ph, dr, imv) in get_side_paths(ct, isd):
          showbug(oph, ct, "attached path/moves do not include a path of {OPHS}")
  
  if ydir != None:
    # Check Y-order of moves:
    for ct in CTS:
      mv0 = contact.side_move(ct, 0)
      mv1 = contact.side_move(ct, 1)
      # Compute move midpoints:
      pm0 = rn.mix(0.5, move.pini(mv0), 0.5, move.pfin(mv0))
      pm1 = rn.mix(0.5, move.pini(mv1), 0.5, move.pfin(mv1))
      # Check order of projections:
      assert rn.dot(pm0, ydir) <= rn.dot(pm1, ydir)

  return
  # ----------------------------------------------------------------------

# SIDE INDICES

def set_side_index(ct, isd, ix):
  assert ix == None or type(ix) is int
  ct.sideix[isd] = ix
  return
  # ----------------------------------------------------------------------
  
def get_side_index(ct, isd):
  return ct.sideix[isd]
  # ----------------------------------------------------------------------

# COOLING ESTIMATORS

def tcool_closed(ct, tcs):
  # Returns the cooling time ratio of {ct} in some path {oph} that closes {ct},
  # given the corresponding cover times {tcs[0..1]}.  Fails if either of 
  # the times is {None}.
  assert type(tcs) is list or type(tcs) is tuple
  assert len(tcs) == 2
  assert tcs[0] != None and tcs[1] != None
  tc = abs(tcs[0] - tcs[1])
  return tc
  # ----------------------------------------------------------------------

def tcool(oph, ct):
  tcs = path_tcovs(oph, ct)
  if tcs[0] == None or tcs[1] == None:
    tc = None
  else:
    tc = tcool_closed(ct, tcs)
  return tc
  # ----------------------------------------------------------------------

def coldest(oph, CTS, n):
  WORST = [] # List of pairs {(ict, tc)} where {ict} is index into {CTS} and {tc} is its cooling time.
  for ict in range(len(CTS)):
    ct = CTS[ict] 
    tc = contact.tcool(oph, ct)
    if tc != None:
      WORST.append((ict, tc))
  WORST.sort(key = lambda x: -x[1])
  if len(WORST) > n: WORST = WORST[0:n]
  CTScold = []
  if len(WORST) > 0:
    sys.stderr.write("coldest contacts (%d):\n" % len(WORST))
    for ict, tc in WORST:
      ct = CTS[ict]
      sys.stderr.write("  %10.6f s %s\n" % (tc, contact.get_name(ct)))
      CTScold.append(ct)
  else:
    assert False, "no contacts are closed by the path."
  return CTScold
  # ----------------------------------------------------------------------

def est_rcool_closed_by_path(ct, i0, oph0, tc0, oph1, tc1, mp_jump):
  # Returns a lower bound for the cooling time ratio of {ct} in some path that includes {oph0} and {oph1}, 
  # in that order and direction; where {oph0} and {oph1} cover sides {i0} and {1-i0} of {ct}, respectively.
  # 
  # Assumes that {tc0} is the fabtime of {oph0} from the midpoint of
  # {ct} to the end of {oph0}, and {tc1} is the fabtime of {oph1} from
  # the beginning to the midpoint of {ct}. Fails if either is {None}.
  tc_lim = tcool_limit(ct)
  if tc_lim == +inf: return 0
  use_links = True
  tconn_est = path.connection_time(oph0, oph1, use_links, mp_jump) 
  tcool_min = tc0 + tconn_est + tc1
  rcool_min = tcool_min/tc_lim 
  return rcool_min
  # ----------------------------------------------------------------------
  
def est_rcool_closed_by_other_paths(ct, oph, isd, tc, mp_jump, quick):
  # Assumes that the oriented path {oph} covers side {isd} of {ct}.
  # Returns a lower bound for the cooling time ratio of {ct} among all
  # paths {P} that start with {oph} and include some paths {oph1} among
  # the paths that have been associated to side {1-isd} of {ct}. The
  # procedure fails if there are no such paths.
  # 
  # Assumes that {tc} (which must not be {None}) is the fabtime of
  # {oph} from the midpoint of {ct} to the end of {oph}.
  # 
  # The {quick} parameter has the same meaning as in {est_rcool}.

  tc_lim = tcool_limit(ct)
  if tc_lim == +inf: return 0

  if quick and tc > tc_lim: return +inf  # Time to finish {oph} already exceeds limit.

  # Compute the min fabtime {tc1} to close {ct} with some choice of {bc} after finishing {oph}:
  OPHS1 = get_side_paths(ct, 1-isd)
  rcool_min = +inf 
  for oph1 in OPHS1:
    tcs1 = path_tcovs(oph1, ct)
    assert tcs1[isd] == None and tcs1[1-isd] != None
    tc1 = tcs1[1-isd]
    rcool_min = min(rcool_min, est_rcool_closed_by_path(ct, isd, oph, tc, oph1, tc1, mp_jump))
  assert rcool_min < +inf, "no paths in {get_side_paths} closes {ct}"
  return rcool_min 
  # ----------------------------------------------------------------------

def est_rcool(oph, ct, mp_jump, quick):
  tclim = tcool_limit(ct)
  if tclim == +inf: 
    rc = 0
  else:
    tcs = path_tcovs(oph, ct)
    assert type(tcs) is list or type(tcs) is tuple
    assert len(tcs) == 2
    if tcs[0] != None and tcs[1] != None:
      # Contact {ct} is closed by {oph}, the exact cooling time is known:
      tc = tcool_closed(ct, tcs)
      rc = tc/tclim
    elif tcs[0] != None or tcs[1] != None:
      # Only one side of {ct} is covered bt {oph}.
      # Get the fabtime from {ct} to the end of {oph}:
      isd = 0 if tcs[0] != None else 1
      tc0 = path.fabtime(oph) - tcs[isd]
      # Compute the cloosing time ratio for the best choice of {bc}
      rc = est_rcool_closed_by_other_paths(ct, oph, isd, tc0, mp_jump, quick)
    else:
      rc = 0
  return rc
  # ----------------------------------------------------------------------

def est_max_rcool(oph, CTS, mp_jump, quick):
  assert type(CTS) is list or type(CTS) is tuple
  est_max_rc = 0
  for ct in CTS:
    rc_ct = est_rcool(oph, ct, mp_jump, quick)
    est_max_rc = max(est_max_rc, rc_ct)
    if quick and est_max_rc > 1: 
      est_max_rc = +inf; break
  return est_max_rc
  # ----------------------------------------------------------------------

def plot_to_files(fname, CTS, clr, dashpat, ext, OPHS, CLRS, rwd, wd_axes, tics, arrows):

  assert type(CTS) is list or type(CTS) is tuple
  assert type(OPHS) is list or type(OPHS) is tuple
  assert len(OPHS) > 0

  # Compute the plot's bounding box:
  B = path.bbox(OPHS)
  B = rn.box_join(B, bbox(CTS))

  dp = None

  # No frame, because it may confuse with contour:
  c, szx, szy = hacks.make_canvas(hacks.round_box(B,0.5), None, None, False, True, 1, 1)

  # Plot the paths:
  axes = True
  dots = True
  arrows_ph = True
  matter = True

  path.plot_standard(c, OPHS, None, None, CLRS, rwd, wd_axes, axes, dots, arrows_ph, matter)
  
  # Plot the contacts:
  wd_ct = 1.5*wd_axes
  sz_tics = wd_ct if tics else 0
  arrows_ct = arrows
  for ct in CTS:
    plot_single(c, ct, None, clr, dashpat, ext, wd=wd_ct, sz_tic=sz_tics, arrow=arrows_ct)
  hacks.write_plot(c, fname)
  return
  # ----------------------------------------------------------------------

def plot_single(c, ct, dp, clr, dashpat, ext, wd, sz_tic, arrow): 
  p = ct.pts[0]
  q = ct.pts[1]
  u, distpq = rn.dir(rn.sub(q, p))
  if distpq < 1.0e-6:
    if ext != 0: sys.stderr.write("!! warning: cannot extend zero-length contact")
    # Perturb the points slightly:
    peps = 0.01*wd
    p = (p[0]-peps, p[1]-peps)
    q = (q[0]+peps, q[1]+peps)
  elif ext != 0:
    p = rn.mix(1,p, -ext,u)
    q = rn.mix(1,q, +ext,u)
  hacks.plot_line(c, p, q, dp, clr, wd, dashpat)

  # Should we plot a transversal tic or arrowhead?
  if sz_tic == None: sz_tic = 0 # Simplification.
  if sz_tic > 0 or arrow:
    # Plot the transversal tic or arrowhead:
    m = rn.mix(0.5, p, 0.5, q)  # Midpoint.
    u = get_perp_dir(m, ct.side[0], ct.side[1])
    sz_arrow = 3*wd if arrow else 0
    # We need a tic with a certain min size for the arrowhead:
    sz_tic = max(sz_tic, 0.80*sz_arrow)
    a = rn.mix(1.0, m, -0.5*sz_tic, u)
    b = rn.mix(1.0, m, +0.5*sz_tic, u)
    sty_basic = [ pyx.style.linecap.round, clr, ]
    sty_tic = sty_basic
    if sz_arrow > 0:
      # Add the arrowhead to the tic.
      arrowpos = 0.5  # Position of arrow on transversal line.
      wd_arrow = sz_arrow/5 # Linewidth for stroking the arrowhead (guess).
      sty_arrow = sty_basic + [ 
        pyx.deco.stroked([pyx.style.linewidth(wd_arrow), pyx.style.linejoin.round]),
        pyx.deco.filled([])
      ]
      sty_tic = sty_tic + \
        [ pyx.deco.earrow(sty_arrow, size=sz_arrow, constriction=None, pos=arrowpos, angle=35) ]
      sys.stderr.write("sz_arrow = %.3f wd_arrow = %3f sz_tic = %.3f\n" % (sz_arrow, wd_arrow, sz_tic))
    
    sty_tic = sty_tic + [ pyx.style.linewidth(wd), ]
    if dp != None: sty_tic.append(pyx.trafo.translate(dp[0], dp[1]))
    c.stroke(pyx.path.line(a[0], a[1], b[0], b[1]), sty_tic)
  return
  # ----------------------------------------------------------------------

def get_perp_dir(m, omv0, omv1):
  # Returns the direction from trace {mv0} towards trace{mv1}
  # at the point {m}, assumed to be the midpoint of a contact
  # between them.

  sys.stderr.write("m = ( %.3f %.3f )\n" % ( m[0], m[1],))
  assert hacks.is_point(m)
  mv0, dr0 = move.unpack(omv0)
  mv1, dr1 = move.unpack(omv1)
  assert mv0 != mv1, "both sides on same move?"
  a = [None,None]
  for imv in range(2):
    mvi = (mv0,mv1)[imv]
    p0i, p1i = move.endpoints(mvi)
    r = min(1, max(0, rn.pos_on_line(p0i, p1i, m))) # Nearest rel pos in move to {m}
    a[imv] = rn.mix(1-r, p0i, r, p1i)
  assert a[0] != a[1]
  sys.stderr.write("a = ( %.3f %.3f ) ( %.3f %.3f )\n" % ( a[0][0], a[0][1], a[1][0], a[1][1],))
  u, da = rn.dir(rn.sub(a[1],a[0]))
  return u
  # ----------------------------------------------------------------------

def has_name(ct):
  return ct.name != None
  # ----------------------------------------------------------------------

def get_name(ct):
  assert isinstance(ct, contact.Contact)
  name = ct.name
  if name == None: name = "C?"
  return name
  # ----------------------------------------------------------------------

def set_name(ct, name):
  assert type(name) is str
  ct.name = name
  return
  # ----------------------------------------------------------------------

def tag_names(CTS, tag):
  if tag != None and tag != "":
    assert type(tag) is str
    for ct in CTS:
      ct.name = tag + get_name(ct)
  return
  # ----------------------------------------------------------------------

def compare(ct0, ct1, tol, die):

  # Check endpoints:
  ends0 = endpoints(ct0)
  ends1 = endpoints(ct1)
  for ipt in range(2):
    if not hacks.same_point(ends0[ipt], ends1[ipt], tol): return fbomb(die)

  for isd in range(2):
    # Check moves on side {isd}:
    omv0 = side_move(ct0, isd)
    omv1 = side_move(ct1, isd)
    if not hacks.same_point(move.pini(omv0), move.pini(omv1), tol): return fbomb(die)
    if not hacks.same_point(move.pfin(omv0), move.pfin(omv1), tol): return fbomb(die)

    # Check paths on side {isd}
    OPHS0 = get_side_paths(ct0, isd)
    OPHS1 = get_side_paths(ct1, isd)
    assert len(OPHS0) == len(OPHS1)
    # ??? Should compare the paths ???
  return True
  # ----------------------------------------------------------------------

def show(wr, pref, ct, suff, wna):
  if pref != None: wr.write(pref)
  wr.write("%-*s" % (wna,get_name(ct)))

  # Endpoints:
  pts = endpoints(ct)
  for ipt in range(2):
    pti = pts[ipt]
    wr.write(" ( %8.3f, %8.3f )" % (pti[0], pti[1]))

  tclim = tcool_limit(ct)
  xtclim = "" if tclim == +inf else ("%.1f" % tclim)
  wr.write(" %5s" % xtclim)

  # Side moves and side paths:
  for isd in range(2):
    wr.write(" ")
    mvi = ct.side[isd]
    wr.write(move.get_name(mvi))
    SDPSi = contact.get_side_paths(ct, isd)
    if len(SDPSi) > 0:
      wr.write("@{")
      sep = ''
      for phi, dri, imvi in SDPSi:
        ophi = path.spin(phi,dri)
        wr.write("%s%s" % (sep, path.get_name(ophi)))
        mvj, drj = move.unpack(path.elem(ophi,imvi))
        xrev = "~" if drj == 1 else ""
        wr.write("[%s%d]" % (xrev,imvi))
        if mvi != mvj:
          sys.stderr.write("\n")
          move.show(sys.stderr, "mvi =  ", mvi, "\n", 0)
          path.show(sys.stderr, "ophi = ", ophi, "\n", True, 0,0,0)
          move.show(sys.stderr, "mvj =  ", mvj, "\n", 0)
          assert mvj == mvi, "inconsistent move index in side path"
        sep = ','
      wr.write("}")
  
  if suff != None: wr.write(suff)
  return
  # ----------------------------------------------------------------------

def show_list(wr, pref, CTS, suff):
  assert type(CTS) is list or type (CTS) is tuple
  nct = len(CTS)
  if nct == 0: return
  wna = 4 # Width of "name" column; min 4 because of the header.
  for ct in CTS: wna = max(wna, len(get_name(ct)))
  wix = len(str(nct-1)) # Num digits in index.
  
  wr.write("\n")

  # Write header:
  wpt = 22 # Width of a point coordinates column.
  wr.write("%*s%*s %-*s %*s %*s %*s sides\n" % (len(pref),'',wix,"k",wna,"name",wpt,"pini",wpt,"pfin",5,"tclim"))
  wr.write("%*s%s %s %s %s %s --------------\n" % (len(pref),'',"-"*wix,"-"*wna,"-"*wpt,"-"*wpt,"-"*5))
  
  # Write contacts:
  for kct in range(len(CTS)):
    ct = CTS[kct]
    if pref != None: wr.write(pref)
    wr.write("%*d " % (wix,kct))
    show(wr, None, ct, suff, wna)
    wr.write("\n")

  wr.write("\n")
  return 
  # ----------------------------------------------------------------------
  
def show_times(wr, OPHS, CTS):
  # Compute and show the fabrication times:
  tfab_tot = 0 # Total fabtime of paths in {OPHS}
  tfab_trc = 0 # Total fabtime of traces in {OPHS}
  for oph in OPHS:
    tfab_tot += path.fabtime(oph)
    OSQS, OJMS = path.split_at_jumps(oph)
    for osq in OSQS:
      tfab_trc += path.fabtime(osq)
  wr.write("total fabtime =  %7.3f\n" % tfab_tot)
  wr.write("extrusion time = %7.3f\n" % tfab_trc)
  wr.write("air time =       %7.3f\n" % (tfab_tot-tfab_trc))
    
  if len(OPHS) == 1:
    # Compute and show the cooling times:
    oph = OPHS[0]
    for ct in CTS:
      wr.write("contact %12s" % contact.get_name(ct))
      tcool = contact.tcool(oph, ct)
      if tcool == None:
        wr.write(" not closed\n")
      else:
        wr.write(" tcool = %7.3f\n" % tcool)
  return
  # ----------------------------------------------------------------------

def write_times (wr, OPHS, CTS, execution_time, Nrast, Trast, Nlink, Tlink, Njump, Tjump):

  fname = "tests/out/" + wr + ".txt"
  sys.stderr.write("writing %s ...\n" % fname)
  wr = open(fname, "w")

  # Compute and show the fabrication times:
  tfab_tot = 0 # Total fabtime of paths in {OPHS}
  tfab_trc = 0 # Total fabtime of traces in {OPHS}
  for oph in OPHS:
    tfab_tot += path.fabtime(oph)
    OSQS, OJMS = path.split_at_jumps(oph)
    for osq in OSQS:
      tfab_trc += path.fabtime(osq)
  wr.write("FabTime: %7.5f\n" % tfab_tot)
  wr.write("ExtTime: %7.5f\n" % tfab_trc)
  wr.write("AirTime: %7.5f\n" % (tfab_tot-tfab_trc))
  wr.write("CpuTime: %7.5f\n" % (execution_time))
  if len(OPHS) == 1:
    # Compute and show the cooling times:
    oph = OPHS[0]
    for ct in CTS:
      #wr.write("contact %12s" % contact.get_name(ct))
      tcool = contact.tcool(oph, ct)
      if tcool == None:
        wr.write("CoolTime: not_closed\n")
      else:
        wr.write("CoolTime: %7.5f\n" % tcool)
  wr.write("Nrast: %d\n" % (Nrast))
  wr.write("RastTime: %7.5f\n" % (Trast))
  wr.write("Nlink: %d\n" % (Nlink))
  wr.write("LnkTime: %7.5f\n" % (Tlink))
  wr.write("Njump: %d\n" % (Njump))
  wr.write("JumpTime: %7.5f\n" % (Tjump))
  wr.close()
  return
  # ----------------------------------------------------------------------


