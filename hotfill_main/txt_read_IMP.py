# Implementation of the module {txt_read}.
# Last edited on 2021-11-25 18:25:14 by stolfi

import txt_read

import block
import move
import move_parms
import path
import contact
import hacks
import rn

import sys
from math import sqrt, sin, cos, atan2, log, exp, floor, ceil, inf, nan, pi

def read(rd, mp_cont, mp_fill, mp_link, angle, shift, fixdir, smptol):
  
  debug = False
  
  Z = None          # {Z} coordinate of slice.
  OPHS = None       # Filling elements (single-raster paths) indexed by {ife}.
  PTSS_cont = None  # List of lists of vertices of contours.
  OLKS = []         # List of link paths.
  CTS = []          # List of contacts.

  ystep = move_parms.width(mp_fill) # Expected spacing of rasters.
  xdir = (cos(angle), sin(angle))   # Dir vector parallel to rasters, "left" to "right".
  ydir = (-xdir[1], xdir[0])        # Dir vector perpendicular to rasters, "up".

  def unrotate_and_shift(p_raw):
    # Given a point as read from the file,
    # separated by the character {sep}, returns the point whith those coordinates (as pair of {float}s) 
    # after rotating by {-alpha} and shifting by {shift}
    #
    xr = rn.dot(p_raw,xdir) + shift[0]
    yr = rn.dot(p_raw,ydir) + shift[1]
    p = (xr,yr)
    return p
    # ....................................................................
  
  def parse_link(vstr):  
    # Takes a list of vertices of a link, encoded as a string {vstr},
    # and returns them a list of points, un-rotated and shifted.
    # The points should be separated by ';', and each point should be two floats separated
    # by '&'.  Alternatively {verts} may be the string 'None' to indicate that there
    # is no link between the two elements; in which case the procedure returns {None}.
    #
    if 'None' in vstr: return None

    # Parse vertices of link path and collect them in the list {VTS}:
    if debug: sys.stderr.write("      vstr = %s\n" % vstr)
    aux = vstr.split(';')
    VTS = [] # Vertices of link.
    for ipt in range(len(aux)):
      if debug: sys.stderr.write("        aux[ipt] = %s\n" % aux[ipt])
      ptx = aux[ipt].split('&') 
      p_raw = (float(ptx[0]), float(ptx[1]))
      p = unrotate_and_shift(p_raw)
      VTS.append(p)
    return VTS
    # ......................................................................

  nread = 0
  nlinks_L = 0 # Number of non-"None" link specs in "L" lines.
  for line in rd:
    nread += 1
    line = line.rstrip()
    if debug: sys.stderr.write("%08d %s\n" % (nread, line))
    code = line[0]
    rest = line[1:].replace('\n', '')

    if code == 'Z':
      # {Z} coordinate of slice.
      Z = float(rest)

    elif code == 'K':
      # Number of contours.
      ncr = int(rest)
      PTSS_cont = []
      for icr in range(ncr): PTSS_cont.append([])

    elif code == 'N':
      # Number of filling raster elements.
      nfe = int(rest) 
      OPHS = [None]*nfe

    elif code == 'C':
      # Vertex of a contour:
      line = line.replace('\n', '')
      c = rest.split(',')
      icr = int(c[0])   # Contour index.
      assert icr >= 0 and icr < ncr
      p_raw = (float(c[1]), float(c[2]))
      p = unrotate_and_shift(p_raw) # Vertex coordinates.
      PTSS_cont[icr].append(p)

    elif code == 'R':
      # Raster element.
      r = rest.split(",")
      ife = int(r[0])                  # Index of raster element.
      assert ife >= 0 and ife < nfe    

      p_raw = (float(r[1]), float(r[2])) # One endpoint.
      p = unrotate_and_shift(p_raw)

      q_raw = (float(r[3]), float(r[4])) # The other endpoint.
      q = unrotate_and_shift(q_raw)

      dr = int(r[5])    # Direction bit for the move in the path.
      igr = int(r[6])   # Group index.
      # Make sure the {Move} object is oriented from left to right:
      if p[0] > q[0]:
        p,q  = q,p
        dr = 1-dr
      if fixdir: dr = 0
      oph = create_single_raster_path(ife, p, q, dr, mp_fill, igr)
      # Save the element in {OPHS} indexed by {ife}:
      OPHS[ife] = oph 

    elif code == 'L':
      # Contact and possibly pair of links between two raster elements
      l = rest.split(",")

      # Get the indices of the two rasters connected by the links
      ife0 = int(l[0].replace('L', '')); oph0 = OPHS[ife0]  # First filing element.
      assert ife0 >= 0 and ife0 < nfe    
      ife1 = int(l[1].replace('L', '')); oph1 = OPHS[ife1]  # Second filing element.
      assert ife1 >= 0 and ife1 < nfe    

      # Checks that the two rasters are in adjacent scan-lines, in proper order:
      if debug: path.show(sys.stderr, "  oph0 = ", oph0, "\n", True, 0,0,0) 
      if debug: path.show(sys.stderr, "  oph1 = ", oph1, "\n\n", True, 0,0,0) 
      check_consecutive_scan_lines(oph0, oph1, ystep)

      # Create contact and attach the paths to it:
      ct = create_and_attach_contact(ife0, oph0, ife1, oph1)
      if ct != None: CTS.append(ct)

      # Get the two links and attack them to the filling elements:
      VTS0 = parse_link(l[2]) # Vertices of first link, or {None}.
      if VTS0 != None: nlinks_L += 1
      olk0 = create_and_attach_link(VTS0, oph0, oph1, smptol, mp_link)
      if olk0 != None: OLKS.append(olk0)

      VTS1 = parse_link(l[3]) # Vertices of second link, or {None}.
      if VTS1 != None: nlinks_L += 1
      olk1 = create_and_attach_link(VTS1, oph0, oph1, smptol, mp_link)
      if olk1 != None: OLKS.append(olk1)
  
  # Assemble the contours:
  OCRS = []   # Contours.
  for icr in range(len(PTSS_cont)):
    PTS = PTSS_cont[icr]
    if PTS != None and len(PTS) != 0:
      assert len(PTS) >= 3
      ocr = path.from_points(PTS + [PTS[0]], mp_cont, None)
      path.set_name(ocr, ("CR%d" % icr), True)
      OCRS.append(ocr)
  path.compute_contour_nesting(OCRS)

  ncr = len(OCRS)
  nfe = len(OPHS)
  nlk = len(OLKS)
  nct = len(CTS)
  sys.stderr.write("read %d lines: %d contours, %d rasters, %d links, %d contacts\n" % (nread,ncr,nfe,nlk,nct))
      
  # Check the links (paranoia):
  assert nlk == nlinks_L
  OLKS_chk = set()  # Set of all {Path} objects of links:
  for ophi in OPHS:
    for ophir in ophi, path.rev(ophi):
      OLKSir = path.get_links(ophir)
      for olkirj in OLKSir:
        lkirj, drirj = path.unpack(olkirj)
        OLKS_chk.add(lkirj)
  OLKS_chk = list(OLKS_chk) # Convert from set to list.
  assert len(OLKS) == len(OLKS_chk)

  path.compute_contour_nesting(OCRS)

  return OCRS, OPHS, OLKS, CTS, Z
  # ----------------------------------------------------------------------

def create_single_raster_path(ife, p, q, dr, mp_fill, igr):
  # Creates a path {ph} with a single raster trace. The trace wiil have
  # endpoints {p,q} and parameters {mp_fill}. The direction will be from
  # {p} to {q} if {dr} is zero, and reversed otherwise.
  # 
  # The line {p--q} must be horizontal, apart from roundoff errors.
  # The raster endpoints will be precisely horizontal.
  #
  # The link and contact information is cleared; namely,
  # {path.get_links(ph)}, {path.get_links(rev(ph))}, and
  # {path.get_contacts(ph,isd)} will be empty. The move's name is set to
  # "R{ife}" and the path's name to "P{ife}".
  # 
  # Also sets the group index of each path to {igr} with {path.set_group}.
  # 
  p_adj, q_adj = check_and_fix_raster_direction(p, q)
  if dr == 0:
    mv = move.make(p_adj, q_adj, mp_fill)
  else:
    mv = move.make(q_adj, p_adj, mp_fill)
  move.set_name(mv, "T%d" % ife)
  ph = path.from_moves((mv,))
  path.set_name(ph, "R%d" % ife, False)
  path.set_group(ph, igr)
  return ph
  # ----------------------------------------------------------------------

def check_and_fix_raster_direction(p, q):
  # Checks that the segment {p,q} is horizontal apart from roundoff
  # error, with {p} to the left of {q}. Returns the coordinates of {p}
  # and {q} adjusted so that they are precisely horizontal.
  
  debug = False
  
  if debug:
    sys.stderr.write("  p = ( %10.6f, %10.6f )" % (p[0], p[1]))
    sys.stderr.write("  q = ( %10.6f, %10.6f )\n" % (q[0], q[1]))
  
  assert abs(p[1] - q[1]) <= 0.0005
  assert p[0] < q[0]
  ym = (p[1] + q[1])/2 # Mean {Y} coordinate.
  
  # sys.stderr.write("  p =     ( %20.16f %20.16f )\n" % (p[0],p[1]))
  # sys.stderr.write("  q =     ( %20.16f %20.16f )\n" % (q[0],q[1]))
  p_adj = (p[0], ym)
  q_adj = (q[0], ym)
  
  if debug:
    sys.stderr.write("  p_adj = ( %20.16f %20.16f )\n" % (p_adj[0],p_adj[1]))
    sys.stderr.write("  q_adj = ( %20.16f %20.16f )\n" % (q_adj[0],q_adj[1]))

  return p_adj, q_adj
  # ----------------------------------------------------------------------

def create_and_attach_link(VTS, oph0, oph1, smptol, mp_link):
  # Creates the link path {olk} between endpoints of two rasters {oph0,oph1}
  # and attaches it to those paths. Returns that link path, or {None}
  # the link was not created.
  # 
  # The parameters {oph0} and {oph1} must be oriented paths that
  # represent adjacent filling elements. The {VTS} parameter should be 
  # the list of vertices of a link path connecting two endpoins, or {None}
  #
  # If {VTS} is {None}, the procedure does nothing.  If {VS} is not {None}, the procedure 
  # makes a path {olk]} from those points, whose traces have parameters
  # {mp_link}. It then attaches {olk]} to the lists of links of {oph0}
  # and {oph1}, using {path.add_link} with the proper orientations.
  # Any links previously attached to {oph0} or {oph1} are preserved.
  
  debug = False
  debug_smp = True
  
  if VTS == None or len(VTS) == 0: return None
  
  # Snap vertices of the link path to the raster element endpoints:
  tol = 0.05  # Snap vertices that are this close to the endpoint (in mm).
  for oph in oph0, oph1:
    for t in path.pini(oph), path.pfin(oph):
      for ipt in range(len(VTS)):
        if hacks.same_point(VTS[ipt], t, 0.05): VTS[ipt] = t

  # Create the moves of the link path:
  TRS = [] # Traces of link.
  for ipt in range(len(VTS) - 1):
    p = VTS[ipt]
    q = VTS[ipt + 1]
    if p != q:
      mv = move.make(p, q, mp_link)
      TRS.append(mv)

  assert len(TRS) > 0, "zero-length link"
  olk_raw = path.from_moves(TRS)
  nmv_raw = path.nelems(olk_raw)

  # Simplify the link:
  if smptol > 0:
    olk = path.simplify(olk_raw, smptol, None)
  else:
    olk = olk_raw
  nmv = path.nelems(olk)
  if nmv != nmv_raw:
    T_raw = path.fabtime(olk_raw)
    T = path.fabtime(olk)
    if debug_smp: 
      sys.stderr.write("link reduced from %d moves (%.6f s) to %d  (%.6f s)\n" % (nmv_raw,T_raw,nmv,T))
  elif nmv != 1:
    T = path.fabtime(olk)
    if debug_smp:
      sys.stderr.write("link has %d moves (%.6f s) - not simplified\n" % (nmv,T))

  # Attach the link path {olk} to the paths:
  nat = 0 # Number of ends that were attached.
  for oph in oph0, path.rev(oph0), oph1, path.rev(oph1):
    for olkr in olk, path.rev(olk):
      if path.pfin(olkr) == path.pini(oph):
        if debug: 
          sys.stderr.write("    attaching link:\n")
          path.show(sys.stderr, "      link: ", olkr, "\n", False, 0,0,0)
          path.show(sys.stderr, "      path: ", oph,  "\n",  False, 0,0,0)
        path.add_link(oph, olkr)
        nat += 1

  assert nat >= 2, "link path does not connect to raster(s)"
  assert nat <= 2, "link path connects to more than 2 raster ends"
  return olk
  # ----------------------------------------------------------------------

def create_and_attach_contact(ife0, oph0, ife1, oph1):
  # Assumes that the oriented paths {oph0} and {oph1} are sngle-trace
  # filling elements with indices {ief0} and {ife1}, in consecutive scan
  # lines, with {oph0} "below" {oph1}.
  # 
  # If the two traces have a non-trivial shared border, creates a
  # {Contact} object between them, and attaches the two paths to it, in both orientations, with
  # {contact.add_side_path} and {path.add_contact}.
  
  assert path.nelems(oph0) == 1
  assert path.nelems(oph1) == 1
  mv0, dr0 = move.unpack(path.elem(oph0,0))
  mv1, dr1 = move.unpack(path.elem(oph1,0))
  mvdir = (1,0)
  
  wd0 = move.width(mv0)
  wd1 = move.width(mv1)
  tol = 0.20*min(wd0,wd1)  # Tolerance for overlaps, tilts, etc.
  ct = contact.from_moves(mv0, mv1, mvdir, 0.1, 0.05, tol)
  if ct == None: 
    sys.stderr.write("!! contact creation failed\n")
    move.show(sys.stderr, "  mv0 = ", mv0, "\n", 0)
    move.show(sys.stderr, "  mv1 = ", mv1, "\n", 0)
  else:
    contact.set_name(ct, "C%d:%d" % (ife0, ife1))

    contact.add_side_path(ct, 0, oph0, 0)
    contact.add_side_path(ct, 1, oph1, 0)

    path.add_contact(oph0, 0, ct)
    path.add_contact(oph1, 1, ct)
  return ct
  # ----------------------------------------------------------------------

def check_consecutive_scan_lines(oph0, oph1, ystep):
  # Assumes that {oph0} and {oph1} are oriented filling raster paths read from the file,
  # already unrotated, each consisting of a single horizontal trace.
  # 
  # Checks whether the two traces are in consecutive scanlines, with {oph0} lower than {oph1},
  # assuming scanlin spacing {ystep}.
  y = [None,None]
  x = [None,None]
  for iph in range(2):
    ophi = (oph0,oph1)[iph]
    x[iph], y[iph] = path.mean_projections(ophi, None, (0,1))
  dy = floor((y[1] - y[0])/ystep + 0.5)
  if dy != 1:
    sys.stderr.write("ystep = %6.3f dy = %20.16f\n" % (ystep,dy))
    for iph in range(2):
      ophi = (oph0,oph1)[iph]
      path.show(sys.stderr, "    ophi%d: " % iph, ophi, None, False, 0,0,0); 
      sys.stderr.write("y = %20.16f\n" % (y[iph]))
      assert False, "rasters are not in consecutive scan-lines"
  return
  # ----------------------------------------------------------------------
  
