# Implementation of module {path}

import path
import move
import move_parms
import contact
import block
import hacks; from hacks import fbomb
import rn 

from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
import pyx 

import sys
from math import nan, inf, sqrt, sin, cos, pi, floor

class Path_IMP:
  # The initial and final points of the {Path} object, in its native
  # orientation, are {ph.endpoints[0]} and {ph.endpoints[1]}, respectively.
  #
  # The list of oriented moves (traces and jumps) that comprise a path
  # {ph}, in the native order, is {ph.OMVS}. Namely, tracing the oriented
  # path {(ph,0)} means executing the oriented moves
  # {ph.OMVS[0],ph.OMVS[1],...,ph.OMVS[nmv-1]}, in this sequence, where
  # {nmv=len(ph.OMVS)}. Tracing the oriented path {(ph,1)} means executing
  # {rev(ph.OMVS[nmv-1]),rev(ph.OMVS[nmv-2]),...,rev(ph.OMVS[0])}, in this
  # sequence.
  #
  # The field {ph.cumtex} is a list such that {ph.cumtex[kmv]} is the
  # cumulative time to fabricate moves {ph.OMVS[0]} through {ph.OMVS[kmv]},
  # for each {kmv} in {0..nmv-1}. These times include the penalty times for
  # each trace/move or move/trace transition. The transition penalty is
  # always assumed to be added  to the execution time of the jump.
  #
  # These times assume that the nozzle is stationary at the beginning
  # and end of each move, so that {ph.cumtex[kmv]} is just the sum of the
  # execution times of the individual moves plus the applicable
  # transition penalties.
  
  def __init__(self, p, q, OMVS, cumtex):
    # Read-only fields:
    self.OMVS = OMVS          # Ordered list of oriented moves.
    self.cumtex = cumtex      # Cumulative execution time up to each move.
    self.endpoints = (p, q)   # Endpoints of the path.
    self.cr_children = None   # List of contour paths that are immediately contained in {self}.
    self.cr_parent = None     # Contour path that imediately conatains {self}.
    self.name = None
    # Mutable fields used for various purposes:
    self.links = [[], []]         # A pair of lists of links that connect to this path.
    self.contacts = [set(),set()] # A pair of sets such that elem {isd} is a set of contacts whose side {isd} is in this {Path}.
    self.group = None             # Group of path in a set of filling elements.
    # Mutable fields used during the {HotPath} heuristic:
    self.block = None             # The owning block of path that is a block choice.
    self.fill_BCS = None          # A list of candidate blocks immediately contained in this contour path.
  # ----------------------------------------------------------------------

# ATTRIBUTES

def nelems(oph):
  ph, dr = unpack(oph)
  return len(ph.OMVS)
  # ----------------------------------------------------------------------
  
def elem(oph, imv): 
  ph, dr = unpack(oph)
  # sys.stderr.write("oph = %s ph = %s dr = %s\n" % (str(oph), str(ph), str(dr)))
  if dr == 0:
    return ph.OMVS[imv]
  else:
    nmv = len(ph.OMVS)
    return move.rev(ph.OMVS[nmv-1-imv])
  # ----------------------------------------------------------------------
    
def find_move(oph, omv):
  nmv = nelems(oph)
  mv, dr = move.unpack(omv)
  for imv in range(nmv):
    omvk = elem(oph,imv)
    mvk, drk = move.unpack(omvk)
    if mvk == mv: return imv
  return None 
  # ----------------------------------------------------------------------

# GEOMETRY

def pini(oph):
  ph, dr = unpack(oph)
  return ph.endpoints[dr]
  # ----------------------------------------------------------------------
  
def pfin(oph):
  ph, dr = unpack(oph)
  return ph.endpoints[1-dr]
  # ----------------------------------------------------------------------

def endpoints(oph):
  ph, dr = unpack(oph)
  if dr == 0:
    return tuple(ph.endpoints)
  else:
    return (ph.endpoints[1], ph.endpoints[0])
  # ----------------------------------------------------------------------
  
def points(oph):
  PTS = [ pini(oph) ]
  for imv in range(nelems(oph)):
    PTS.append(move.pfin(elem(oph, imv)))
  return PTS
  # ----------------------------------------------------------------------

def bbox(OPHS):
  B = None
  for oph in OPHS:
    ph, dr = unpack(oph)
    B = rn.box_include_point(B, pini(oph)) # In case the path is trivial.
    B = rn.box_join(B, move.bbox(ph.OMVS))
  return B
  # ----------------------------------------------------------------------

def find_nearest_move_endpoint(oph, p):
  nmv = nelems(oph)
  if nmv > 0 and tuple(pini(oph)) == tuple(pfin(oph)):
    # Closed non-trivial path. No need to include {pfin(oph)}:
    imv_min = None
    d_min = +inf
  else:
    # Open or trivial path. Consider {pfin(oph)} explicitly:
    imv_min = nmv
    d_min = rn.dist(p, pfin(oph))
  for imv in range(nmv):
    q = move.pini(elem(oph, imv))
    d = rn.dist(p, q)
    if d < d_min:
      imv_min = imv
      d_min = d
  assert imv_min != None and d_min < +inf
  return imv_min, d_min
  # ----------------------------------------------------------------------

def find_nearest_midline_point(oph, p):
  nmv = nelems(oph)
  if nmv > 0 and tuple(pini(oph)) == tuple(pfin(oph)):
    # Closed non-trivial path. No need to include {pini(oph)}:
    imv_min = None
    d_min = +inf
    t_min = None
    q_min = None
  else:
    # Open or trivial path. Consider {pini(oph)} explicitly:
    imv_min = 0
    q_min = pini(oph)
    d_min = rn.dist(p, q_min)
    t_min = 0.0
  t = 0.0 # Distance along midline from {pini(oph)} to {pini(elem(oph,imv))}
  for imv in range(nmv):
    mvi = elem(oph, imv)
    q0 = move.pini(mvi)
    q1 = move.pfin(mvi)
    di, qi = rn.pt_seg_dist(p, q0, q1)
    if di < d_min:
      imv_min = imv
      d_min = di
      q_min = qi
      t_min = t + rn.dist(q0, qi)
    t += rn.dist(q0, q1)
  assert imv_min != None and d_min < +inf and t_min <= t
  return q_min, imv_min, t_min, d_min
  # ----------------------------------------------------------------------

def find_midline_point_by_length(oph, t):

  if t < 0 : return None, None
  if t == 0 : return 0, pini(oph)

  nmv = nelems(oph)

  imv = None
  q = None

  t_prev = 0 # Total length of previous moves.
  for jmv in range(nmv):
    omvj = elem(oph, jmv)
    ej = move.length(omvj)
    tj = t_prev + ej
    if tj > t:
      imv = jmv
      r = max(0, min(1, (t - t_prev)/ej))
      q = rn.mix(1-r,move.pini(omvj), r,move.pfin(omvj))
      break
    t_prev = tj
  if imv != None: assert imv >= 0 and imv < nmv # Paranoia.
  
  return imv, q
  # ----------------------------------------------------------------------

def extract_section(oph, imv0,q0, imv1,q1):
  nmv = nelems(oph)

  # Correct for pot-extremal cases:
  if q0 == pini(oph) and imv0 == -1: imv0 = 0
  if q1 == pini(oph) and imv1 == -1: imv0 = 0
  if q0 == pfin(oph) and imv0 == nmv: imv0 = nmv-1
  if q1 == pfin(oph) and imv1 == nmv: imv1 = nmv-1
  
  pname = path.get_name(oph)
  if pname != None and len(pname) == 0: pname = None

  if imv0 > imv1: 
    # Crossed over. Don't bother to check for zero-length moves.
    return None
  elif imv0 == imv1:
    omv = elem(oph,imv0)
    # ??? should check whether {q0} and {q1} are on that move ???
    if q0 == q1: return make_trivial(q0)
    pa = move.pini(omv); pz = move.pfin(omv)
    if q0 == pa and q1 == pz: return from_moves([omv,])
    # Check whether {q0} is ahead of {q1}:
    if rn.dot(rn.sub(q1, q0), rn.sub(pz, pa)) < 0: return None
    # New move:
    omv_new = move.make(q0, q1, move.parameters(omv))
    move.set_name(omv_new, pname + (".X%d" % imv0))
    return from_moves([omv_new,])
  else:
    # New ends are in distinct moves:
    OMVS = []
    
    # Take first move or part thereof:
    omv0 = elem(oph, imv0); p0a = move.pini(omv0); p0z = move.pfin(omv0)
    if q0 == p0a:
      # Take whole of {omv0}
      OMVS.append(omv0)
    elif q0 != p0z:
      # Take a piece of it:
      mv0_trim = move.make(q0, p0z, move.parameters(omv0))
      move.set_name(mv0_trim, pname + (".X%d" % imv0))
      OMVS.append(mv0_trim)
    else:
      # Take nothing of it:
      pass
      
    # Copy the moves strictly beween {imv0} and {imv1}:
    for imv in range(imv0+1,imv1):
      OMVS.append(elem(oph,imv))
      
    # Take last move or part thereof:
    omv1 = elem(oph, imv1); p1a = move.pini(omv1); p1z = move.pfin(omv1)
    if q1 == p1z:
      # Take whole of {omv1}
      OMVS.append(omv1)
    elif q1 != p1a:
      # Take a piece of it:
      mv1_trim = move.make(p1a, q1, move.parameters(omv1))
      move.set_name(mv1_trim, pname + (".X%d" % imv1))
      OMVS.append(mv1_trim)
    else:
      # Take nothing of it:
      pass

    # Make the new path:
    if len(OMVS) == 0:
      # Must be {q0} at end of one move and {q1} at start of the next one:
      assert q0 == p0z and q1 == p1a
      sec = make_trivial(q0)
    else:
      sec = from_moves(OMVS)
    return sec
  # ----------------------------------------------------------------------

def trim(oph, d0, d1):

  assert d0 >= 0 and d1 >= 0, "trim amounts cannot be negative"
  
  nmv = nelems(oph)
  if nmv == 0:
    if d0 == 0 and d1 == 0:
      return make_trivial(pini(oph))
    else:
      return None
  
  # Find the endpoints {q_end[0..1]} of the new path and the indices {imv_end[0..1]}
  # of the moves that contain them. Makes sure that any move before {imv_end[0]}
  # or after {imv_end[1]} is to be completely excluded.
  q_end = [None,None]  # Endpoints of new path.
  imv_end = [None,None] # Indices of moves that contain those ends.
  for kdr in 0, 1:
    ophk = spin(oph, kdr)
    dk = (d0, d1)[kdr]
    imvk, qk = find_midline_point_by_length(ophk,dk)
    if imvk == None: return None
    assert imvk >= 0 and imvk < nmv
    q_end[kdr] = qk
    while imvk < nmv and qk == move.pfin(elem(ophk,imvk)): imvk += 1
    imv_end[kdr] = imvk if kdr == 0 else nmv - 1 - imvk

  # Build the list {OMVS} of moves of the new path:
  imv0 = imv_end[0]; q0 = q_end[0]
  imv1 = imv_end[1]; q1 = q_end[1]
  
  oph_new = extract_section(oph, imv0,q0, imv1,q1)
  return oph_new
  # ----------------------------------------------------------------------

def mean_projections(oph, xdir, ydir):
  p = pini(oph)
  q = pfin(oph)
  if xdir != None:
    xp = rn.dot(p, xdir) 
    xq = rn.dot(q, xdir)
    xm = (xp + xq)/2
  else:
    xm = None
    
  if ydir != None:
    yp = rn.dot(p, ydir)
    yq = rn.dot(q, ydir)
    ym = (yp + yq)/2
  else:
    ym = None
  return ( xm, ym )
  # ----------------------------------------------------------------------

# CREATION

def make_trivial(p):
  return path.Path(p, p, (), ())
  # ----------------------------------------------------------------------
  
def from_move(p, q, mp):
  mv = move.make(p, q, mp)
  return path.Path(p, q, ( mv, ), ( move.fabtime(mv), ) )
  # ----------------------------------------------------------------------

def from_moves(OMVS):
  assert type(OMVS) is list or type(OMVS) is tuple
  nmv = len(OMVS) # Number of moves.
  assert nmv >= 1, "must have at least one move"
  p = move.pini(OMVS[0])   # Initial point of final path.
  q = move.pfin(OMVS[nmv-1]) # Final point of final path.
  p_prev = p
  t_prev = 0 
  cumtex = []
  for imv in range(nmv):
    omvk = OMVS[imv]
    mvk, drk = move.unpack(omvk) # For type checking.
    pk, qk = move.endpoints(omvk)
    assert pk == p_prev, "moves are not connected: %s %s" % (str(p_prev),str(pk))
    tfink = t_prev + move.fabtime(omvk)
    if move.is_jump(omvk):
      udk = move.ud_penalty(omvk)
      # Add the applicable move/jump transition penalties:
      if imv > 0 and not move.is_jump(OMVS[imv-1]):
        tfink += udk
      if imv < nmv-1 and not move.is_jump(OMVS[imv+1]):
        tfink += udk
    cumtex.append(tfink)
    p_prev = qk
    t_prev = tfink
  assert p_prev == q
  return path.Path(p, q, tuple(OMVS), tuple(cumtex))
  # ----------------------------------------------------------------------
  

def from_points(PTS, mp_trace, mp_jump):
  assert type(PTS) is list or type(PTS) is tuple
  m = len(PTS) # Number of points or lists.
  assert m >= 1, "must have at least one point"
  p_prev = None
  OMVS = []
  cumtex = []
  jmp = None # True if previous element of {PTS} was list, false if it was point.
  for ptj in PTS:
    if hacks.is_point(ptj):
      # Append a single move with the given width:
      if p_prev != None: 
        mv = move.make(p_prev, ptj, mp_trace)
        OMVS.append(mv)
      p_prev = ptj
      jmp = False
    else:
      # {ptj}  must be a list of points. Appends a sequence of moves from {p_prev}
      # to those points.  The first move, if any, will be a jump, the rest 
      # will be traces.  Then forces a jump to whatever comes next (if any).
      assert type(ptj) is list or type(ptj) is tuple
      mpjk = mp_jump # Parameters of of next move.
      for ptjk in ptj:
        assert hacks.is_point(ptjk)
        if p_prev != None: 
          mv = move.make(p_prev, ptjk, mpjk)
          OMVS.append(mv)
        p_prev = ptjk
        mpjk = mp_trace
      mpjk = mp_jump
  assert p_prev != None # Since there must have been at least one point.
  if len(OMVS) == 0:
    return make_trivial(p_prev)
  else:
    return from_moves(OMVS)
  # ----------------------------------------------------------------------

def concat(OPHS, use_links, mp_jump):
  assert type(OPHS) is list or type(OPHS) is tuple
  assert type(use_links) is bool
  assert mp_jump == None or isinstance(mp_jump, move_parms.Move_Parms)
  
  m = len(OPHS) # Number of paths.
  assert m >= 1, "must have at least one path"
  p = pini(OPHS[0])   # Initial point of final path.
  q = pfin(OPHS[m-1]) # Final point of final path.
  oph_prev = None
  p_prev = None
  OMVS = []
  for ophj in OPHS:
    phj, drj = unpack(ophj) # For type checking.
    nmvj = nelems(ophj)
    pj = pini(ophj)
    qj = pfin(ophj)
    if p_prev != None and p_prev != pj:
      # Get the connector {olk}, possibly {None} or a jump:
      cnj = get_connector(oph_prev, ophj, use_links, mp_jump)
      if cnj == None:
        # No connector needed:
        pass
      elif isinstance(cnj, move.Move):
        # Connector is a jump, inset it:
        assert move.pini(cnj) == p_prev
        assert move.pfin(cnj) == pj
        OMVS.append(cnj)
        p_prev = pj
      else:
        # Connector must be a link path, insert its moves:
        assert pini(cnj) == p_prev
        assert pfin(cnj) == pj
        # Insert a copy of the link path 
        for kmv in range(nelems(cnj)):
          omvk = elem(cnj,kmv)
          assert move.pini(omvk) == p_prev
          OMVS.append(omvk)
          p_prev = move.pfin(omvk)

    # Append the moves of {ophj}:
    for kmv in range(nmvj):
      omvk = elem(ophj, kmv)
      pk, qk = move.endpoints(omvk)
      assert p_prev == None or p_prev == pk
      OMVS.append(omvk)
      p_prev = qk
    oph_prev = ophj

  assert p_prev == q
  
  # Create a path from all moves:
  ph = from_moves(OMVS)
  
  # Copy links from paths of {OPHS} to {ph}:
  set_concat_links(ph, OPHS)
  set_concat_links(rev(ph), [rev(ophi) for ophi in reversed(OPHS)])
  
  # Copy the contact sets from {OPHS}:
  for isd in range(2):
    for ophi in OPHS:
      for ct in get_contacts(ophi, isd):
        add_contact(ph, isd, ct)
  return ph
  # ----------------------------------------------------------------------
  
def get_connector(oph0, oph1, use_links, mp_jump):
  p0 = pfin(oph0)
  p1 = pini(oph1)
  if p0 == p1:
    conn = None
  else:
    olk = get_connecting_link(oph0, oph1) if use_links else None
    if olk != None:
      assert pini(olk) == p0
      assert pfin(olk) == p1
      conn = olk
    else:
      # Create and insert a jump:
      jmp = move.make(p0, p1, mp_jump)
      conn = jmp
  return conn
  # ----------------------------------------------------------------------

def set_concat_links(oph, OPHS):
  # Assumes that {oph} was obtained by concatenating the paths of
  # {OPHS}. Sets the input links of {oph} to the appropriate input links
  # of the paths of {OPHS}. Basically all input links of all paths in
  # {OPHS} up to the first non-trivial one, excluding duplicates and
  # those links that start on endpoints of moves of {oph} itself.
  
  LKS = [] # List of {Path} objects of links already added.
  nph = len(OPHS)
  iph = 0
  while iph < nph:
    ophi = OPHS[iph]
    OLKSi = get_links(ophi)
    for olkij in OLKSi:
      # Check if already added:
      ok = True
      lkij, drij = unpack(olkij)
      for lkr in LKS:
        if lkij == lkr: ok = False; break
      if ok:
        add_link(oph, olkij) # No-op if {olkij} starts on {oph}.
        LKS.append(lkij)
    # Check if no more links after these:
    if nelems(ophi) > 0: break
    iph += 1
  return
  # ----------------------------------------------------------------------

def displace(oph, ang, v):
  nmv = nelems(oph)
  if (nmv == 0):
    p = pini(oph)
    return make_trivial(p)
  else:
    OMVS = []
    for kmv in range(nmv):
      omvk = elem(oph, kmv)
      mpk = move.parameters(omvk)
      OMVS.append(move.displace(omvk, ang, v, mpk))
    opd = from_moves(OMVS)
    return opd
  # ----------------------------------------------------------------------

# ORIENTATION

def rev(oph):
  ph, dr = unpack(oph)
  return (ph, 1-dr)
  # ----------------------------------------------------------------------

def spin(oph, dr):
  ph1, dr1 = unpack(oph)
  return (ph1, (dr1 + dr) % 2)
  # ----------------------------------------------------------------------

def unpack(oph):
  assert oph != None
  # sys.stderr.write("oph = %s\n" % str(oph))
  if isinstance(oph, path.Path):
    # sys.stderr.write("returning %s, 0\n" % str(oph))
    return oph, 0
  else:
    assert type(oph) is tuple
    assert len(oph) == 2
    ph, dr = oph
    # sys.stderr.write("ph = %s dr = %s\n" % (str(ph), str(dr)))
    assert isinstance(ph, path.Path)
    assert dr == 0 or dr == 1
    return ph, dr
  # ----------------------------------------------------------------------

def obj_in_list(ph, OPHS):
  # Returns {True} iff the oriented path {OPHS} has a an element with the {Path}
  # object {ph}, in any orientationa.
  for ophi in OPHS:
    phi, dri = unpack(ophi)
    if phi == ph: return True
  return False
  # ----------------------------------------------------------------------

# TIMING

def fabtime(oph):
  ph, dr = unpack(oph)
  nmv = len(ph.OMVS)
  return 0 if nmv == 0 else ph.cumtex[nmv-1]
  # ----------------------------------------------------------------------
  
def tini(oph, imv):
  ph, dr = unpack(oph)
  nmv = len(ph.OMVS)
  assert imv >= 0 and imv <= nmv
  if imv == 0: return 0
  if imv == nmv: return fabtime(oph)
  if dr == 0:
    # Native direction:
    return ph.cumtex[imv-1]
  else:
    # Reverse direction
    return ph.cumtex[nmv-1] - ph.cumtex[nmv-1-imv] 
  # ----------------------------------------------------------------------
  
def tfin(oph, imv):
  return tini(oph, imv+1)
  # ----------------------------------------------------------------------
  
# INCREMENTAL PATH BUILDING

def move_to(MVS, p, q, mp, tag):
  if p != None: 
    mv = move.make(p, q, mp)
    move.set_name(mv, "M%s%d" % (tag,len(MVS)))
    MVS.append(mv)
  return q
  # ......................................................................

def finish(PHS, MVS, tag):
  assert len(MVS) != 0
  ph = from_moves(MVS)
  set_name(ph, "P%s" % tag, False)
  PHS.append(ph)
  MVS.clear()
  return
  # ......................................................................

# LINK PATHS

# The list {ph.links[0]} has all link paths that end
# at {pini(ph)}.  The list {ph.links[1]} has all link paths
# that start at {pfin(ph)}.

def clear_links(oph):
  ph, dr = unpack(oph)
  ph.links[dr] = [].copy()
  return
  # ----------------------------------------------------------------------
  
def add_link(oph, olk):

  debug = False

  # The link must end where {oph} starts:
  assert tuple(pfin(olk)) == tuple(pini(oph))

  # The link must not be a loop:
  p = tuple(pini(olk))
  assert tuple(pfin(olk)) != p

  # Ignore the link if it starts on {oph} itself:
  for kmv in range(nelems(oph)):
    if p == tuple(move.pfin(elem(oph,kmv))):
      if debug:
        sys.stderr.write("link ignored: "); show(sys.stderr, olk, False, 0,0,0); sys.stderr.write("\n")
      return
  if debug:
    sys.stderr.write("link added: "); show(sys.stderr, olk, False, 0,0,0); sys.stderr.write("\n")
          
  # Seems OK, add the link:
  ph, dr = unpack(oph)
  if dr == 0:
    # The link ends at {pini(ph)}:
    ph.links[0].append(olk)
  else:
    # The link ends at {pfin(ph)}
    ph.links[1].append(rev(olk))
  return
  # ----------------------------------------------------------------------
  
def set_links(oph, OLKS):
  clear_links(oph);
  for olk in OLKS: add_link(oph, olk)
  return
  # ----------------------------------------------------------------------

def get_links(oph):
  ph, dr = unpack(oph)
  if dr == 0:
    OLKS = ph.links[0]
  else:
    OLKS = [rev(olk) for olk in ph.links[1]]
  return OLKS
  # ----------------------------------------------------------------------
  
def get_connecting_link(oph0, oph1):
  OLKS0 = get_links(rev(oph0))
  OLKS1 = get_links(oph1)
  for olk0 in OLKS0:
    for olk1 in OLKS1:
      if (unpack(rev(olk0))) == (unpack(olk1)):
        olk = olk1
        assert pini(olk) == pfin(oph0)
        assert pfin(olk) == pini(oph1)
        return olk
  return None
  # ----------------------------------------------------------------------

def get_all_connecting_links(oph0, oph1):
  OLKS = []
  for oph0r in oph0, rev(oph0):
    for oph1r in oph1, rev(oph1):
      olk = get_connecting_link(oph0r, oph1r)
      if olk != None:
        OLKS.append(olk)
  return OLKS
  # ----------------------------------------------------------------------

def connection_time(oph0, oph1, use_links, mp_jump):
  assert oph0 != None and oph1 != None
  # Get the moves that need to be joined:
  n0 = nelems(oph0); n1 = nelems(oph1)
  omv0 = None if n0 == 0 else elem(oph0,n0-1)
  omv1 = None if n1 == 0 else elem(oph1,0)
  p = pfin(oph0); q = pini(oph1)
  if p == q:
    tconn = move.transition_penalty(omv0, omv1)
  else:
    if use_links:
      olk = get_connecting_link(oph0, oph1)
    else:
      olk = None
    if olk != None:
      tconn = fabtime(olk)
      # Assume that the link starts and ends with trace:
      if move.is_jump(omv0): tconn += move.ud_penalty(omv0)
      if move.is_jump(omv1): tconn += move.ud_penalty(omv1)
    else:
      tconn = move_parms.nozzle_travel_time(rn.dist(p, q), None, mp_jump)
      udp = move_parms.ud_penalty(mp_jump)
      if move.is_trace(omv0): tconn += udp
      if move.is_trace(omv1): tconn += udp
  return tconn
  # ----------------------------------------------------------------------

# CONTACT GRAPH

def clear_contacts(oph ,isd):
  ph, dr = unpack(oph)
  ph.contacts[isd] = set()
  return 
  # ----------------------------------------------------------------------

def add_contact(oph, isd, ct):
  ph, dr = unpack(oph)
  ph.contacts[isd].add(ct)
  return  
  # ----------------------------------------------------------------------

def get_contacts(oph, isd):
  ph, dr = unpack(oph)
  return ph.contacts[isd].copy()
  # ----------------------------------------------------------------------

# GROUPS

def set_group(oph, igr):
  assert igr == None or (type(igr) is int and igr >= 0)
  ph, dr = unpack(oph)
  ph.group = igr
  return
  # ----------------------------------------------------------------------
  
def get_group(oph):
  ph, dr = unpack(oph)
  return ph.group
  # ----------------------------------------------------------------------

def separate_by_group(OPHS):
  ngr = 0       # Number of distinct groups seen so far.
  maxgr = None  # Largest group index seen so far.
  GRPHS = []      # {GRPHS[igr]} is list of indices of all elems of group {igr}.
  for oph in OPHS:
    # Obtain and validate {igr}
    igr = get_group(oph) 
    if igr != None:
      assert type(igr) is int and  igr >= 0
      # Extend {GRPHS} to include {GRPHS[igr]}:
      while igr >= len(GRPHS): GRPHS.append([])
      # Add element to orig group list:
      if len(GRPHS[igr]) == 0:
        ngr += 1
      GRPHS[igr].append(oph)
  return GRPHS, ngr
  # ----------------------------------------------------------------------

# CONTOUR PATH NESTING

def compute_contour_nesting(OCRS):

  ncr = len(OCRS)
  POLS = [None]*ncr # Contour paths converted to {shapely} polygons.
  
  # Clear all nesting lists and convert all contour paths to polygons:
  for icr in range(ncr):
    ocr = OCRS[icr]
    ph, dr = unpack(ocr)
    assert pini(ph) == pfin(ph)
    ph.cr_children = []
    ph.cr_parent = None
    PTS = [ move.pini(elem(ocr,kmv)) for kmv in range(nelems(ocr)) ]
    POLS[icr] = Polygon(PTS)
  
  # Find the parent of each contour:
  for icr1 in range(ncr):
    # Find the index {ipar} of the parent of {OCRS[icr1]}:
    ipar = None
    pg1 = POLS[icr1]
    for icr2 in range(ncr):
      pg2 = POLS[icr2]
      if icr1 != icr2:
        # See if {icr2} can be the parent of {icr1}:
        if pg2.contains(pg1):
          if ipar == None:
            # First container of {icr1}, keep it for now:
            ipar = icr2
          else:
            # We already have a container {ipar}. Check which is closer:
            pg3 =  POLS[ipar] 
            if pg3.contains(pg2):
              # {icr2} is a closer container:
              ipar = icr2
            else:
              assert pg2.contains(pg3), "contours are not properly nested"
    if ipar != None:
      ocr1 = OCRS[icr1]
      ph1, dr1 = unpack(ocr1)
      ocr3 = OCRS[ipar]
      ph1.cr_parent = ocr3
    
  # Now define the children of each contour:
  for ocr1 in OCRS:
    ph1, dr1 = unpack(ocr1)
    
    ocr2 = ph1.cr_parent
    if ocr2 != None:
      ph2, dr2 = unpack(ocr2)
      ph2.cr_children.append(ocr1)
  # ----------------------------------------------------------------------

def inner_contours(ocr):
  ph, dr = unpack(ocr)
  return ph.cr_children
  # ----------------------------------------------------------------------

def outer_contour(ocr):
  ph, dr = unpack(ocr)
  return ph.cr_parent
  # ----------------------------------------------------------------------

def contour_nesting(ocr0, ocr1):
  # sys.stderr.write("ocr0 = %s\n" % str(ocr0))
  # sys.stderr.write("ocr1 = %s\n" % str(ocr1))

  ph0, dr0 = unpack(ocr0)
  ph1, dr1 = unpack(ocr1)
  if ph0 == ph1: return 0
  
  ph0a = ph0
  ph1a = ph1
  while ph0a != None or ph1a != None:
    if ph0a != None:
      # Walk {ph0a} up:
      ocr0a = ph0a.cr_parent
      if ocr0a == None:
        ph0a = None
      else:
        ph0a, dr0a = unpack(ocr0a)
        if ph0a == ph1:
          # {ph0} is contained in {ph1}
          return -1
    if ph1a != None:
      # Walk {ph1a} up:
      ocr1a = ph1a.cr_parent
      if ocr1a == None:
        ph1a = None
      else:
        ph1a, dr1a = unpack(ocr1a)
        if ph1a == ph0:
          # {ph1} is contained in {ph0}
          return +1
  return 0
  # ----------------------------------------------------------------------

def shift_contour(ocr, kmv):
  assert pini(ocr) == pfin(ocr)
  nmv = nelems(ocr)
  OMVS = [ elem(ocr, (imv + kmv) % nmv) for imv in range(nmv) ]
  ocrk = from_moves(OMVS)
  assert pini(ocrk) == pfin(ocrk)
  return ocrk
  # ----------------------------------------------------------------------

# PLOTTING

def plot_to_files(fname, OPHS, CLRS, rwd, wd_axes, grid, deco):

  assert type(OPHS) is list or type(OPHS) is tuple
  assert len(OPHS) > 0

  # Compute the plot's bounding box:
  B = bbox(OPHS)
  sys.stderr.write("bounding box = [ %8.2f _ %8.2f ]" % (B[0][0],B[1][0]))
  sys.stderr.write(" x [ %8.2f _ %8.2f ]\n" % (B[0][1],B[1][1]))
  
  dp = rn.sub((5,5), B[0])

  c, szx,szy = hacks.make_canvas(hacks.round_box(B,0.5), dp, None, False, grid, 1, 1)

  axes = deco
  dots = deco
  arrows = True # ??? Separate parameter? Also {deco}? ???
  matter = deco

  plot_standard(c, OPHS, dp, None, CLRS, rwd, wd_axes, axes, dots, arrows, matter)
  hacks.write_plot(c, fname)
  return
  # ----------------------------------------------------------------------

def plot_standard(c, OPHS, dp, layer, CLRS, rwd, wd_axes, axes, dots, arrows, matter, arrow_sz):

  assert type(OPHS) is list or type(OPHS) is tuple
  nph = len(OPHS)
  if nph == 0: return

  assert wd_axes != None and wd_axes > 0, "invalid wd_axes"
  
  if CLRS == None: 
    CLRS = [ pyx.color.rgb(0.050, 0.800, 0.000), ] # Default trace color.
  else:
    assert type(CLRS) is list or type(CLRS) is tuple
  nclr = len(CLRS)
  assert nclr == 1 or nclr == nph

  def pick_colors(kph):
    # Returns the colors for trace sausages and axes of path {OPHS[kph]}.
    if nclr == 1:
      ctraces = CLRS[0]
    else:
      ctraces = CLRS[kph]
    caxes =   pyx.color.rgb(0.6*ctraces.r, 0.6*ctraces.g, 0.6*ctraces.b) # Color of trace axis, dots, arrow.
    return ctraces, caxes

  # Colors:
  cmatter = hacks.matter_color(None)  # Material footprint color.
  cjumps =  pyx.color.rgb.black  # Color of jump axes, dots, arrows.

  # Dimensions relative to nominal trace widths:
  rwd_matter = 1.13; # Material footprint.

  # Absolute dimensions (mm):
  wd_dots =   2.5*wd_axes;  # Dots at ends of moves.
  sz_arrows = 5.0*wd_axes     # Size of arrows. 

  # Get list {lys} of layers to plot:
  if layer == None:
    # Plot all four layers.
    lys = (range(4))
  else:
    # Plot only the selected layer.
    assert type (layer) is int
    assert layer >= 0 and layer < 4
    lys = (layer,)

  # Plot the layers:
  for ly in lys:
    for kph in range(nph):
      oph = OPHS[kph]
      ctraces, caxes = pick_colors(kph)
      if ly == 0 and matter:
        # plots the estimate of actual material:
        plot_layer\
          ( c, oph, dp, jmp=False, \
            clr=cmatter, rwd=rwd_matter, wd=0, 
            dashed=False, wd_dots=0, sz_arrows=None
          )
      elif ly == 1:
        # plots the nominal trace material:
        plot_layer \
          ( c, oph, dp, jmp=False,
            clr=ctraces, rwd=rwd, wd=0, 
            dashed=False, wd_dots=0, sz_arrows=None
          )
      elif ly == 2 and (axes or dots or arrows):
        # Plots the trace axis:
        t_wd = wd_axes if axes else 0
        t_wd_dots = wd_dots if dots else 0
        t_sz_arrows = sz_arrows if arrows else 0
        plot_layer \
          ( c, oph, dp, jmp=False, 
            clr=caxes, rwd = 0, wd=t_wd, 
            dashed=False, wd_dots=0, sz_arrows=None
          )
      elif ly == 3:
        # Plots the jump:
        j_wd = 2.5*wd_axes;
        j_wd_dots = 3.0*wd_dots;
        j_sz_arrows = 2.5*sz_arrows
        plot_layer \
          ( c, oph, dp, jmp=True, 
            clr=cjumps, rwd = 0, wd=j_wd, 
            dashed=True, wd_dots=j_wd_dots, sz_arrows=j_sz_arrows*arrow_sz
          )
  # ----------------------------------------------------------------------

def plot_layer(c, oph, dp, jmp, clr, rwd, wd, dashed, wd_dots, sz_arrows):
  # Get the {Path} object {ph} and the direction {dr} ({None} if un-oriented):
  ph, dr = unpack(oph)

  if clr == None: return

  # Simplfications:
  if rwd == None: rwd = 0
  if wd == None: wd = 0
  if wd_dots == None: wd_dots = 0
  if sz_arrows == None: sz_arrows = 0

  assert rwd >= 0 and wd >= 0 and wd_dots >= 0 and sz_arrows >= 0

  for imv in range(nelems(oph)):
    omvk = elem(oph, imv)
    if move.is_jump(omvk) == jmp:
      wdk = rwd*move.width(omvk) + wd
      move.plot_layer(c, omvk, dp, clr, wdk, dashed, wd_dots, sz_arrows)
  return
  # ----------------------------------------------------------------------

def plot_single(c, oph, dp, split, clr):
  ph, d = unpack(oph)  # For typechecking.

  if split:
    assert clr == None
    OPHS, OJMS = split_at_jumps(oph)
    CLRS = hacks.trace_colors(len(OPHS), None)
  else:
    CLRS = [clr,] if clr != None else None 
    OPHS = [oph,]
    OJMS = []

  nph = len(OPHS)
  njm = len(OJMS)
  assert nph > 0

  wd_axes = 0.05 * move.width(elem(OPHS[0], 0))

  if nph > 0:
    # Plot the continuous raster sequences:
    if CLRS == None: CLRS = hacks.trace_colors(nph, None)
    rwd = 0.80; # Relative width of trace sausages.
    axes = True
    dots = True
    arrows = True
    matter = True
    plot_standard(c, OPHS, dp, None, CLRS, rwd, wd_axes, axes, dots, arrows, matter) 

  if njm > 0:
    # Plot the jumps:
    clr_jmp = pyx.color.rgb.black
    rwd = 0
    axes = True
    dots = True
    arrows = True
    matter = False
    move.plot_standard(c, OJMS, dp, None, [clr_jmp,], rwd, wd_axes, axes, dots, arrows, matter) 
  return
  # ----------------------------------------------------------------------

def split_at_jumps(oph):
  
  OPHS = [] # Complete trace-only paths found so far.
  OJMS = [] # Jumps found so far.
  OTRS = [] # Traces of last (possibly incomplete) trace-only path found so far.
  
  nmv = nelems(oph)
  for kmv in range(nmv + 1):
    omvk = elem(oph, kmv) if kmv < nmv else None
    if omvk == None or move.is_jump(omvk): 
      if len(OTRS) != 0:
        ophi = from_moves(OTRS)
        OPHS.append(ophi)
        OTRS = []
      if omvk != None: OJMS.append(omvk)
    else:
      OTRS.append(omvk)
  assert len(OTRS) == 0
  
  return OPHS, OJMS
  # ----------------------------------------------------------------------

def simplify(oph, tol, mp_jump):
  
  nmv = nelems(oph)
  
  # Extracts the points of the path, for convenience:

  # Alocate the dynamic programming tableau {P,J}:
  P = [ None ]*(nmv+1)
  J = [ None ]*(nmv+1)
  # For each {i} in {0..nmv}, element {P[i]} will be the best simplified
  # path for the path {ophi} consisting of the first first {i} moves of
  # {oph}. If {i>0}, that path is path {P[J[i]]} plus one last move from
  # {pfin(P[J[i]])} {pfin(ophi)}; or just {P[J[i]]} if those two points
  # coincide.

  q = [None ]*(nmv+1) # Convenience: {q[j]} is the same as {pfin(P[j])},
  T = [None ]*(nmv+1) # Convenience: {T[j]} is {fabtime(P[j])}.
  
  q[0] = pini(oph)
  T[0] = 0
  P[0] = from_points((q[0],), None, None)
  J[0] = -1
  for i in range(1,nmv+1):
    assert i >= 1 and i <= nmv
    omvi = elem(oph,i-1)
    mpi = move.parameters(omvi)
    q[i] = move.pfin(omvi)
    
    def extensible(j):
      # Checks whether the path {P[j]} can be extended with zero or one
      # moves to satisfactorily approximate all the moves {elem(oph,k)}
      # for {k} in {j..i-1}. All those moves must have the same
      # parameter vector {mpi} and, if {mpi} is not a jump, any vertices
      # other than {q[j]} and {q[i]} must be within {tol} of the segment
      # {q[j]--q[i]}.
      assert j < i
      for dk in range(i - j):
        k = j + dk
        assert j <= k and k < i
        omvk = elem(oph, k)
        # Check whether the move parameters are compatible:
        mpk = move.parameters(omvk)
        if move_parms.is_jump(mpi):
          if not move_parms.is_jump(mpk): return False
        else:
          if mpk != mpi: return False
        if k > j:
          # Check whether the point {q[k]} is close enough to the 
          # segment from {q[j]--q[i]}:
          dk, qk = rn.pt_seg_dist(q[k], q[j], q[i])
          # sys.stderr.write("    j = %d k = %d i = %d dk = %.4f qk = ( %.6f %.6f )\n" % (j,k,i,dk,qk[0],qk[1]))
          if dk > tol: return False
      return True
      # ....................................................................
    
    # Compute {P[i]} that ends at {qi} by adding one move to {P[0..i-1]}:
    Jbest = None # Index of path {P[j]} that is to be extended
    Pbest = None # The extended path.
    Tbest = +inf # Fabtime of the extended path.
    for j in range(i):
      if extensible(j):
        if q[j] == q[i]:
          # No need to add another move:
          Pji = P[j]; Tji = T[j]; 
        else:
          # Add another move:
          if move_parms.is_jump(mpi):
            # Replace all original jumps by new jumps with the given jump parameter
            omvji = move.make(q[j], q[i], mp_jump)
          elif j == i-1:
            # Reuse the original trace:
            omvji = omvi
          else:
            # Needs a new trace, but with the same params as all merged ones:
            omvji = move.make(q[j], q[i], mpi)
            move.set_name(omvji, "T%d.%d" % (j, i))
          Pji = concat([P[j], from_moves([omvji,])], False, None) # No links or jumps needed.
          Tji = fabtime(Pji)
        if Tji < Tbest:
          # sys.stderr.write("    replaced best for i = %d: j = %d (T = %.3f)\n" % (i, j, Tji))
          Pbest = Pji; Tbest = Tji; Jbest = j
      else:
        # sys.stderr.write("  j = %d not extensible to i = %d\n" % (j,i))
        pass
    # sys.stderr.write("  best for i = %d: j = %d (T = %.3f)\n" % (i, Jbest, Tbest))
    assert Pbest != None # Since {ophi} is the worst case.
    P[i] = Pbest
    J[i] = Jbest
    T[j] = Tbest
  osp = P[nmv]
  assert pini(osp) == pini(oph)
  assert pfin(osp) == pfin(oph)
  
  # Copy path name:
  set_name(osp, get_name(oph), moves=False)
  
  # Copy groups index:
  set_group(osp, get_group(oph))
  
  # Copy link attachments:
  for idr in 0, 1:
    ophr = spin(oph, idr)
    ospr = spin(osp, idr)
    set_links(ospr, get_links(ophr))
  
  return osp
  # ----------------------------------------------------------------------

# DEBUGGING AND TESTING

def check_links(OPHS,OLKS):
  OLKS_att = set() # Link {Path} objects attached to the paths.
  for oph in OPHS:
    ph, dr = unpack(oph)
    # Check incoming and outgoing links:
    for drk in range(2):
      ophk = spin(oph,drk)
      for olk in get_links(ophk):
        assert pfin(olk) == pini(ophk)
        lki, dri = unpack(olk)
        OLKS_att.add(lki)
  assert OLKS_att == set(OLKS)
  return
  # ----------------------------------------------------------------------

def validate(oph):
  nmv = nelems(oph)
  ph, dr = unpack(oph)
  assert len(ph.OMVS) == nmv
  assert len(ph.cumtex) == nmv
  # sys.stderr.write("dr = %d cumtex = %s\n" % (dr, str(ph.cumtex)))
  # Validate the list of moves:
  p_prev = pini(oph)  # Final point of previous move.
  t_prev = 0 # Cumulative execution time.
  for imv in range(nmv):
    omvi = elem(oph, imv)
    # sys.stderr.write("omvi = %s\n" % str(omvi))
    mvi, dri = move.unpack(omvi)
    # sys.stderr.write("mvi = %s dri = %s\n" % (str(mvi), str(dri)))

    # Check for duplicate moves:
    for jmv in range(imv):
      mvj, drj = move.unpack(elem(oph,jmv))
      assert mvi != mvj, "duplicate {Move} objects in path"

    # Check for gaps:
    mv_pini, mv_pfin = move.endpoints(omvi)
    # sys.stderr.write("p_prev = %s pini = %s\n" % (str(p_prev), str(mv_pini)))
    assert p_prev == mv_pini

    # Check fabtimes table:
    mv_tini = tini(oph, imv)
    mv_tfin = tfin(oph, imv)
    mv_tex = move.fabtime(omvi)
    # Account for up/down penalties:
    if move.is_jump(omvi):
      # Add the applicable move/jump transition penalties:
      udi = move.ud_penalty(omvi)
      if imv > 0 and not move.is_jump(elem(oph,imv-1)):
        mv_tex += udi
      if imv < nmv-1 and not move.is_jump(elem(oph,imv+1)):
        mv_tex += udi
    # sys.stderr.write("imv = %d" % imv)
    # sys.stderr.write(" t_prev = %12.8f tini = %12.8f tfin = %12.8f" % (t_prev,mv_tini,mv_tfin))
    # sys.stderr.write(" tex = %12.8f\n" % mv_tex)
    assert abs(t_prev - mv_tini) < 1.0e-8
    assert abs((mv_tfin - mv_tini) - mv_tex) < 1.0e-8
    p_prev = mv_pfin
    t_prev = mv_tfin
  assert p_prev == pfin(oph)
  assert abs(t_prev - fabtime(oph)) < 1.0e-8
  # ----------------------------------------------------------------------

def compare(oph0, oph1, tol, die):
  debug = False

  if debug: 
    show(sys.stderr, "oph0: ", oph0, "\n", True, 0,0)
    show(sys.stderr, "oph1: ", oph1, "\n", True, 0,0)
  nmv = nelems(oph0)
  if nelems(oph1) != nmv: return fbomb(die)
  if not hacks.same_point(pini(oph0), pini(oph1), tol): return fbomb(die)
  for imv in range(nmv):
    omv0 = elem(oph0, imv)
    omv1 = elem(oph1, imv)
    if not hacks.same_point(move.pini(omv0), move.pini(omv1), tol): return fbomb(die)
    if not hacks.same_point(move.pfin(omv0), move.pfin(omv1), tol): return fbomb(die)
  return True
  # ----------------------------------------------------------------------

def has_name(oph):
  ph, dr = unpack(oph)
  return ph.name != None
  # ----------------------------------------------------------------------

def get_name(oph):
  ph, dr = unpack(oph)
  name = ph.name
  if name == None: name = "P?"
  if dr == 1: name = "~" + name
  return name
  # ----------------------------------------------------------------------

def set_name(oph, name, moves):
  assert type(name) is str
  ph, dr = unpack(oph)
  ph.name = name
  if moves:
    nmv = nelems(ph) # Number of moves.
    ndg = len(str(nmv-1)) # Number of digits to use.
    for kmv in range(nmv):
      omv = elem(ph, kmv)
      move.set_name(omv, "%s.%0*d" % (name,ndg,kmv))
  return
  # ----------------------------------------------------------------------

def tag_names(OPHS, tag):
  if tag != None and tag != "":
    assert type(tag) is str
    for oph in OPHS:
      ph, dr = unpack(oph)
      ph.name = tag + get_name(ph)
  return
  # ----------------------------------------------------------------------

def show(wr, pref, oph, suff, moves, wna,wnm,wgr):
  ph, dr = unpack(oph)
  if pref != None: wr.write(pref)
  wr.write(" " if dr == 0 else "~")
  wr.write("%-*s" % (wna-1,get_name(ph)))
  nmv = nelems(oph)
  wr.write(" %*d" % (wnm,nmv))
  igr = get_group(oph)
  if igr != None:
    wr.write(" %*d" % (wgr,igr))
  else:
    wr.write(" %*s" % (wgr,"-"))
  p = pini(oph)
  q = pfin(oph)
  wr.write(" (%6.1f,%6.1f)" % (p[0],p[1]))
  wr.write(" (%6.1f,%6.1f)" % (q[0],q[1]))

  # Show names of links:
  nlk = 0 # Number of links printed.
  for dr in 0, 1:
    ophr = spin(oph, dr)
    for olk in get_links(ophr):
      wr.write(" " if nlk == 0 else ",");
      wr.write("<" if dr == 0 else ">")
      wr.write(get_name(spin(olk, dr)))
      nlk += 1
    
  # Show names of contacts:
  # Get the two sets and their union:
  CTS0 = get_contacts(oph, 0)
  CTS1 = get_contacts(oph, 1)
  CTS = set.union(CTS0, CTS1)
  nct = 0 # Number of contacts printed.
  for ct in CTS:
      wr.write(" " if nct == 0 else ",")
      wr.write(contact.get_name(ct))
      if not ct in CTS1:
        wr.write(":0")
      elif not ct in CTS0:
        wr.write(":1")
      else:
        wr.write("*")
      nct += 1

  if moves:
    for imv in range(nmv):
      omv = elem(oph, imv)
      wr.write(" " if imv == 0 else ",")
      wr.write(move.get_name(omv))
  if suff != None: wr.write(suff)
  return
  # ----------------------------------------------------------------------

def show_list(wr, pref, OPHS, suff, moves, links):
  assert type(OPHS) is list or type (OPHS) is tuple
  nph = len(OPHS)
  if nph == 0: return
  wna = 3 # Width of "name" column; min 3 because of the header.
  wnm = 1 # Width of "number of moves" column.
  wgr = 1 # Width of "group" column
  for oph in OPHS: 
    ph, dr = unpack(oph)
    wna = max(wna, len(get_name(ph)))
    wnm = max(wnm, len(str(nelems(oph))))
    igr = get_group(oph)
    if igr != None: wgr = max(wgr, len(str(igr)))
    wgr 
  wix = len(str(nph-1)) # Num digits in index.
  wna += 1 # Because of possible "~"
  
  wr.write("\n")

  # Write header:
  mvtit = "links, contacts" + (", moves" if moves else "")
  if pref != None: wr.write("%*s" % (len(pref),''))
  wr.write("%*s %-*s %*s %*s" % (wix,"k",wna,"name",wnm,"m",wgr,'g'))
  wr.write(" %*s %*s %s\n" % (15,"pini",15,"pfin",mvtit))

  # Write dashes line:
  if pref != None: wr.write("%*s" % (len(pref),''))
  wr.write("%s %s %s %s" % ("-"*wix,"-"*wna,"-"*wnm,"-"*wgr))
  wr.write(" %s %s %s\n" % ("-"*15,"-"*15,"-"*len(mvtit)))
  
  # Write paths:
  for kph in range(len(OPHS)):
    oph = OPHS[kph]
    if pref != None: wr.write(pref)
    wr.write("%*d " % (wix,kph))
    show(wr, None, oph, suff, moves, wna,wnm,wgr)
    wr.write("\n")
    
  if links:
    # Collect set {LKS} of the {Path} objects of all links:
    LKS = set() 
    for oph in OPHS:
      for olkj in get_links(oph) + get_links(rev(oph)):
        lkj, drj = unpack(olkj)
        LKS.add(lkj)
    # Print the links:
    show_list(wr, pref, list(LKS), suff, moves, False)

  wr.write("\n")
  return 
  # ----------------------------------------------------------------------
  
