# Implementation of module {path_example}
# Last edited on 2021-10-29 14:26:28 by stolfi

import path_example
import move
import move_example
import move_parms
import contact
import path
import raster
import rootray
import rootray_cart
import hacks
import rn
import sys
import random
from math import sqrt, hypot, sin, cos, floor, ceil, inf, nan, pi

def raster_rectangle(plo, axis, n, alt, sz, step, mp_trace,mp_jump):
  
  # Create all needed rasters, links, and jumps:
  TRS = move_example.rectangle_rasters(plo, axis, n, sz, step, mp_trace)
  if alt:
    LKS0, LKS1 = move_example.rectangle_links(plo, axis, n, sz, step, mp_trace)
    LJS0 = LKS0; LJS1 = LKS1
  else:
    JMS0, JMS1 = move_example.rectangle_jumps(plo, axis, n, sz, step, mp_jump)
    LJS0 = JMS0; LJS1 = JMS1
  
  # Select the moves for each path:
  MVS0 = []
  MVS1 = []
  for kmv in range(n):
    # Append the raster line in the proper orientation:
    trk = TRS[kmv]
    if alt:
      # Alternating directions:
      MVS0.append(move.spin(trk, kmv%2))
      MVS1.append(move.spin(trk, (kmv+1)%2))
    else:
      # Uniform directions:
      MVS0.append(trk)
      MVS1.append(move.rev(trk))
    if kmv <  n-1:
      # Append the proper link or jump:
      if alt:
        # Join with links:
        if kmv%2 == 0:
          MVS0.append(LKS1[kmv])
          MVS1.append(LKS0[kmv])
        else:
          MVS0.append(LKS0[kmv])
          MVS1.append(LKS1[kmv])
      else:
        # Join with jumps:
        MVS0.append(JMS0[kmv])
        MVS1.append(JMS1[kmv])
        
  ph0 = path.from_moves(MVS0); path.set_name(ph0, "P0", False)
  ph1 = path.from_moves(MVS1); path.set_name(ph1, "P1", False)

  PHS = [ ph0, ph1, ]

  for oph in PHS: path.validate(oph)

  return PHS, TRS, LJS0, LJS1
  # ----------------------------------------------------------------------

def spiral_rectangle(pini, szx, szy, axis, mp_trace):
  wd = move_parms.width(mp_trace)
  MVS = [] # The moves
  
  p = pini          # Current point.
  ax = axis         # Axis of next move.
  d = [ szx, szy ]  # Signed length of next moves along each axis.
  ntr = 0 # Number of traces made.
  # Decide the displacement vectors for the next moves in each direction:
  while True:
    # Is there space for another {ax} trace:
    if abs(d[1-ax]) < 0.45*wd: break
    # Take the next move along the axis {ax} with signed length {d[ax]}:
    q = ( p[0] + d[0], p[1] ) if ax == 0 else ( p[0], p[1] + d[1] )
    mv = move.make(p, q, mp_trace)
    xd = "a" if d[ax] > 0 else "b"
    move.set_name(mv, "HV"[ax] + xd + str(ntr//4))
    MVS.append(mv)
    ntr += 1
    if ntr > 2:
      # Subtract the trace from the area left to cover:
      d1 = d[1-ax]
      d[1-ax] = d1 - wd if d1 > 0 else d1 + wd
    # Reverse the direction of the next {ax} trace:
    d[ax] = -d[ax]
    # Prepare for the next iteration:
    ax = 1-ax
    p = q
  ph = path.from_moves(MVS)
  pini = path.pini(ph)

  path.validate(ph)
  
  return ph
  # ----------------------------------------------------------------------

def onion(ctr, Rc, mp_cont, Rf, mp_fill, phase, nt, mp_jump):

  wdc = move_parms.width(mp_cont)
  wdf = move_parms.width(mp_fill)
  
  def reduce_nt(nt1, R1, wd1):
    # Given the number of sides {nt1} in the previous outer element (filler
    # or contour), and the radius {R1} of the vertices of the next element,
    # returns a reduced number of sides {nt1/d} for some small
    # {d}, provided that the gaps between the polygon and the circle
    # are not too large compared to the trace width {wd1}.   If it fails,
    # returns {nd1} unchanged.
    while True:
      nt1_old = nt1
      # Find a small divisor {idiv} of {nt_next}
      for idiv in 7,5,3,2:
        nt_try = nt1//idiv
        if nt1 % idiv == 0 and nt_try >= 3:
          astep_try = 2*pi/nt_try
          # Check the rel deviation {sagg} of trace edges from circle would be acceptable:
          sagg = R1*(1 - cos(astep_try/2))/wd1
          # sys.stderr.write("nt1 = %d idiv = %d sagg/wd1= %.4f\n" % (nt1, idiv, sagg/wd1))
          if sagg <= 0.20: 
            # OK, reduce {nt1} and keep trying:
            nt1 = nt_try
            sys.stderr.write("R = %.3f sagg = %.3f nt = %d\n" % (R1,sagg,nt1))
      if nt1 == nt1_old: 
        # Tried all divisors but failed to reduce:
        return nt1
    assert False
    # ----------------------------------------------------------------------

  assert Rc > 0 and wdc > 0
  assert wdf > 0
  assert nt >= 3 
  astep = 2*pi/nt
  # Start with the contour:
  gap = wdc
  phc = circle(ctr, Rc, nt, phase, mp_cont, gap)
  CIRPHS = [ phc, ]
  Rin = Rc*cos(astep/2) - wdc/2  # Bullseye's inradius.
  Rex = Rc + wdc/2               # Bullseye's exradius.
  R_prev = Rc        # Endpoint radius of previous element.
  wd_prev = wdc      # Width of previous element
  phase_prev = phase # Phase angle of previous element.
  nt_prev = nt       # Number of sides in previous element.
  if Rf != Rc:
    while True:
      # Compute elements {R_next,nt_next,phase_next} of next circle.
      # Also update {Rin} or {Rex}.
      astep_prev = 2*pi/nt_prev # Angular increment of previous element.
      nt_next = nt_prev # May be reduced.
      ashift_prev = - ceil(2*wdf/R_prev/astep_prev)*astep_prev
      phase_next = phase_prev + ashift_prev
      if Rf < Rc:
        # Filling inwards. The radius and phase shift do not depend on {nt_next}:
        R_next = R_prev - (wd_prev + wdf)/2/cos(astep/2) # Endpoint radius of next element.
        # Try to reduce the number of segments{nt_prev}:
        nt_next = reduce_nt(nt_next, R_next, wdf)
        # Recompute the step between vertices and the phase shift:
        astep_next = 2*pi/nt_next
        # Compute inradius{Rin_next} of next element:
        Rin_next = R_next*cos(astep_next/2) - wdf/2 
        if R_next < 0 or Rin_next < Rf:
          # No space for the next element
          break
        Rin = Rin_next
      elif Rf > Rc:
        astep_next = 2*pi/nt_next
        R_next = R_prev + (wd_prev + wdf)/2/cos(astep_next/2) # Endpoint radius of next element.
        Rot_next = R_next + wdf/2  # Exradius of next element.
        if Rot_next > Rf:
          # No space for the next element
          break
        Rex = Rot_next
      else:
        assert False
      sys.stderr.write("adding R = %.3f nt = %d\n" % (R_next,nt_next))
      gap = wdf
      phk = circle(ctr, R_next, nt_next, phase_next, mp_fill, gap)
      CIRPHS.append(phk)
      R_prev = R_next
      wd_prev = wdf
      phase_prev = phase_next 
      nt_prev = nt_next
  use_links = False
  ph = path.concat(CIRPHS, use_links, mp_jump)

  path.validate(ph)

  return ph, Rin, Rex
  # ----------------------------------------------------------------------

def gearloose(R, zigzag, mp_cont, mp_fill, mp_jump):

  wdf = move_parms.width(mp_fill)
  wdc = move_parms.width(mp_cont)
  
  Rex = R        # Outer radius of gear.
  Rin = 0.80*Rex # Inner radius of gear.
  ctr = (Rex+1,Rex+1)

  ph0 = path.from_points((ctr, ctr), mp_fill, mp_jump)
  
  ph1 = gear(ctr, Rin, Rex, 11, True, pi/7, mp_cont, mp_jump)
  
  Sin = 0.20*R  # Inner radius of annulus.
  gap2 = wdc
  ph2 = circle(ctr, Sin, 33, pi/11, mp_cont, gap2)

  Sot = 0.60*R  # Outer radius of annulus.
  gap3 = wdc
  ph3 = circle(ctr, Sot, 44, -pi/11, mp_cont, gap3)
      
  ang = 3*pi/11
  Fin =  -1 if Sin < 0 else Sin + wdc/2
  Fot = Sot - wdc/2
  dmin = -Sot
  dmax = +Sot
  step = 2*wdf if zigzag else wdf
  
  side = 0
  ph4a = ring_fill(Fin, Fot, dmin, dmax, step, zigzag, side, mp_fill, mp_jump)
  assert ph4a != None
  ph4 = path.displace(ph4a, ang, ctr)  
  assert ph4 != None

  use_links = False
  ph = path.concat((ph0, ph1, ph2, ph3, ph4), use_links, mp_jump)
  
  path.validate(ph)
  
  return ph
  # ----------------------------------------------------------------------

# MISCELLANEOUS PATHS FOR THE UNIT TEST PROGRAMS

def misc_A(mp_trace1, mp_trace2, mp_jump):

  p0 = (2,2)
  p1 = (3,3)
  p2 = (2,4)
  p3 = (1,3)

  mv01 = move.make(p0, p1, mp_trace1); move.set_name(mv01, "T0")
  mv12 = move.make(p1, p2, mp_jump)  ; move.set_name(mv12, "J0")
  mv23 = move.make(p2, p3, mp_trace2); move.set_name(mv23, "T1")

  ph = path.from_moves((mv01, mv12, mv23))
  
  path.set_name(ph, "P", True)
  
  path.validate(ph)

  return ph
  # ----------------------------------------------------------------------

def misc_B(mp_trace, mp_jump):

  p01 = (1,1)
  p02 = (1,4)
  p03 = (2,4)

  p04 = (3,5)
  p05 = (4,0)

  p06 = (5,1)
  p07 = (7,1)
  p08 = (7,2)
  p09 = (6,2)
  p10 = (5,3)
  p11 = (7,3)
  p12 = (7,4)
  p13 = (5,4)

  ptss = ((p01, p02, p03), (p04, p05), (p06, p07, p08, p09, p10, p11, p12, p13))
  ph = path.from_points(ptss, mp_trace, mp_jump)
  
  path.set_name(ph, "P", True)
  
  path.validate(ph)

  return ph
  # ----------------------------------------------------------------------

def misc_C(mp_trace, mp_link, mp_jump):

  # Path trace endpoints:
  
  p0a = (1, 2)
  p1a = (2, 1)
  p2a = (3, 1)
  p3a = (4, 2)

  p0b = (4, 4)
  p1b = (3, 3)

  p0c = (2, 4)

  p0d = (1, 3)
  p1d = p0d

  p0e = p0c
  
  p0f = (5, 2)
  p1f = (5, 1)
  
  # Main paths:

  pha = path.from_points( ((p0a, p1a,), (p2a, p3a,)), mp_trace, mp_jump)
  phb = path.from_points( ((p0b, p1b,),), mp_trace, mp_jump)
  phc = path.from_points( (p0c,), mp_trace, mp_jump)
  phd = path.from_points( ((p0d, p1d,),), mp_trace, mp_jump)
  phe = path.from_points( ((p0e,),), mp_trace, mp_jump)
  phf = path.from_points( (p0f, p1f,), mp_trace, mp_jump)

  PHS = [ pha, phb, phc, phd, phe, phf, ]
  
  for iph in range(len(PHS)):  
    path.set_name(PHS[iph], ("P" + "abcdefgkijklm"[iph]), True)

  # Create and attach link paths:
  
  for ph1 in PHS:
    for ph2 in PHS:
      if ph1 != ph2:
        for oph2 in ph2, path.rev(ph2):
          p = path.pfin(ph1)
          q = path.pini(oph2)
          if p != q and rn.dist(p, q) <= 2.5:
            # Create a link path {lk} from {pfin(ph1)} to {pini(oph2)}:
            m = rn.mix(0.5, p, 0.5, q) # Midpoint of segment {p--q}.
            lk = path.from_points([p, m, q], mp_link, None)
            path.set_name(lk, "[" + path.get_name(ph1) + ":" + path.get_name(oph2) + "]", True)
            # Attach it to the paths:
            path.add_link(oph2, lk)
            path.add_link(path.rev(ph1), path.rev(lk))

  for ph in PHS: path.validate(ph)
  
  return PHS
  # ----------------------------------------------------------------------

def misc_D(mp_trace, mp_jump):

  TRS, JMS = move_example.misc_B(mp_trace, mp_jump)

  pha = path.from_moves((TRS[0], TRS[1], TRS[2],)) # Bottom is LR
  phb = path.from_moves((move.rev(TRS[0]), move.rev(TRS[3]), move.rev(TRS[2]),)) # Bottom is RL
  phc = path.from_moves((TRS[0], TRS[1], TRS[2], TRS[3],)) # Closed.

  phd = path.from_moves((TRS[4],JMS[0],TRS[5],))
  phe = path.from_moves((TRS[5],JMS[1],TRS[4],))

  phf = path.from_moves((TRS[6],TRS[7]))

  path.set_name(pha, "Pa", False)
  path.set_name(phb, "Pb", False)
  path.set_name(phc, "Pc", False)
  path.set_name(phd, "Pd", False)
  path.set_name(phe, "Pe", False)
  path.set_name(phf, "Pf", False)

  OPHS = [ pha, phb, path.rev(phc), phd, phe, path.rev(phf), ]

  for oph in OPHS: path.validate(oph)
  
  return OPHS, TRS, JMS
  # ----------------------------------------------------------------------

def misc_E(mp_trace, mp_jump):

  TRS, JMS = move_example.misc_C(mp_trace, mp_jump)

  # Join the moves in paths:

  opha = path.from_moves((TRS[0],JMS[0],TRS[6],))
  ophb = path.from_moves((move.rev(TRS[2]),JMS[1],TRS[5],JMS[2],TRS[4],))
  ophc = path.rev(path.from_moves((TRS[1],JMS[3],TRS[3])))
  ophd = path.from_moves((TRS[7],JMS[4],TRS[8],))
  ophe = path.from_moves((TRS[9],))

  path.set_name(opha, "Pa", False)
  path.set_name(ophb, "Pb", False)
  path.set_name(ophc, "Pc", False)
  path.set_name(ophd, "Pd", False)
  path.set_name(ophe, "Pe", False)

  OPHS = [ opha, ophb, ophc, ophd, ophe, ]

  for oph in OPHS: path.validate(oph)
  
  return OPHS, TRS, JMS
  # ----------------------------------------------------------------------

def misc_F(alt, mp_trace, mp_jump):

  wd = move_parms.width(mp_trace)

  p11 = (1, 1 + 0*wd) # Starting nozzle position.

  p21 = (4, 1 + 1*wd)
  p22 = (7, 1 + 1*wd)

  p31 = (4, 1 + 2*wd)
  p32 = (8, 1 + 2*wd)

  p41 = (1, 1 + 3*wd)
  p42 = (2, 1 + 3*wd)
  p43 = (4, 1 + 3*wd)
  p44 = (5, 1 + 3*wd)
  p45 = (7, 1 + 3*wd)
  p46 = (9, 1 + 3*wd)

  p51 = (1, 1 + 4*wd)
  p52 = (2, 1 + 4*wd)
  p53 = (4, 1 + 4*wd)
  p54 = (5, 1 + 4*wd)
  p55 = (7, 1 + 4*wd)
  p56 = (8, 1 + 4*wd)

  p61 = (1, 1 + 5*wd)
  p62 = (5, 1 + 5*wd)

  if alt:
    pts= \
      ( ( p11, ),
        ( p21, p22, p32, p31, ),
        ( p41, p42, ), ( p43, p44, ), ( p45, p46, p56, p55, ), ( p54, p53, ), ( p52, p51, p61, p62, ),
        # ( p41, p42, ), ( p43, p44, ), ( p45, p46, ), ( p54, p53, ), ( p52, p51, p61, p62, ),
      )
  else:
    pts = \
       ( ( p11, ),
        ( p21, p22, ),
        ( p31, p32, ),
        ( p41, p42, ), ( p43, p44, ), ( p45, p46, ),
        ( p51, p52, ), ( p53, p54, ),
        ( p55, p56, ),
        # ( p55, p56, ),
        ( p61, p62, ),
      )
  ph = path.from_points(pts, mp_trace, mp_jump)
 
  path.set_name(ph, "P", True)
 
  path.validate(ph)

  return ph
  # ----------------------------------------------------------------------

def misc_G(mp_cont, mp_fill, mp_jump):
  wdf = move_parms.width(mp_fill)
  wdc = move_parms.width(mp_cont)
  
  ya = 1
  yb = ya + (wdf+wdc)/2
  yc = yb + wdc

  pa0 = (  0, ya )
  qa0 = (  2, ya )
  pa1 = (  4, ya )
  qa1 = ( 11, ya )
  
  pb0 = (  1, yb )
  qb0 = (  8, yb )
  pb1 = ( 10, yb )
  qb1 = ( 11, yb )
  
  pc0 = (  0, yc )
  qc0 = (  2, yc )
  pc2 = (  7, yc )
  qc2 = ( 11, yc )
  
  tra0 = move.make(pa0, qa0, mp_fill); move.set_name(tra0, "TGa0")
  tra1 = move.make(pa1, qa1, mp_fill); move.set_name(tra1, "TGa1")
  trb0 = move.make(pb0, qb0, mp_cont); move.set_name(trb0, "TGb0")
  trb1 = move.make(pb1, qb1, mp_cont); move.set_name(trb1, "TGb1")
  trc0 = move.make(pc0, qc0, mp_cont); move.set_name(trc0, "TGc0")
  trc2 = move.make(pc2, qc2, mp_cont); move.set_name(trc2, "TGc2")
  
  TRS02 = [ tra0, tra1, trc0, trc2, ]
  TRS1 = [ trb0, trb1, ]
  
  jm_qa0_pc0 = move.make(qa0, pc0, mp_jump); move.set_name(jm_qa0_pc0, "JGqa0pc0")
  jm_qc0_pa1 = move.make(qc0, pa1, mp_jump); move.set_name(jm_qc0_pa1, "JGqc0pa1")
  jm_qa1_pc2 = move.make(qa1, pc2, mp_jump); move.set_name(jm_qa1_pc2, "JGqa1pc2")
  
  jm_qb0_pb1 = move.make(qb0, pb1, mp_jump); move.set_name(jm_qb0_pb1, "JGqb0pb1")

  jm_pa1_qa0 = move.make(pa1, qa0, mp_jump); move.set_name(jm_pa1_qa0, "JGpa1qa0")
  jm_pa0_pc0 = move.make(pa0, pc0, mp_jump); move.set_name(jm_pa0_pc0, "JGpa0pc0")
  jm_qc0_pc2 = move.make(qc0, pc2, mp_jump); move.set_name(jm_qc0_pc2, "JGqc0pc2")
  
  MVS0 = [ tra0, jm_qa0_pc0, trc0, jm_qc0_pa1, tra1, jm_qa1_pc2, trc2, ]
  MVS1 = [ trb0, jm_qb0_pb1, trb1, ]
  MVS2 = [ move.rev(tra1), jm_pa1_qa0, move.rev(tra0), jm_pa0_pc0, trc0, jm_qc0_pc2, trc2, ]
  
  ph0 = path.from_moves(MVS0)
  ph1 = path.from_moves(MVS1)
  ph2 = path.from_moves(MVS2)
  
  PHS = [ ph0, ph1, ph2 ]

  for iph in range(len(PHS)): path.set_name(PHS[iph], "PG%d" % iph, False)

  for ph in PHS: path.validate(ph)
 
  return PHS, TRS02, TRS1
  # ----------------------------------------------------------------------

def misc_H(mp_trace):
  OPHS = []
  xA = -3
  xB = +3
  R = sqrt(2)
  for iph in range(4):
    t = 2*pi*(iph/4 + 0.125)
    s = 2*pi*(iph/4 + 0.0625)
    ct = cos(t); st = sin(t)
    pA = ( xA + R*ct, R*st)
    pB = ( xB - R*ct, R*st)
    h = 2-ct
    pM = (0, h*(1 + 0.3*h)*R*st)
    if iph % 2 == 0:
      oph = path.from_points((pA,pM,pB,),mp_trace,None)
    else:
      oph = path.rev(path.from_points((pB,pM,pA,), mp_trace,None))
    path.set_name(oph, "P%d" % iph, False)
    for kmv in range(path.nelems(oph)): move.set_name(path.elem(oph,kmv), "T%d%d" % (iph,kmv))
    OPHS.append(oph)

  for oph in OPHS: path.validate(oph)
 
  return OPHS
  # ----------------------------------------------------------------------

def misc_J(mp_cont, mp_fill):
  ctr = (5, 5)
  
  R0 = 4
  nt0 = 8
  phase0 = 0.1
  ph0 = circle(ctr, R0, nt0, phase0, mp_cont, 0.0)
  assert path.pini(ph0) == path.pfin(ph0)
  path.set_name(ph0, "P0", True)
  
  R1 = 2
  nt1 = 6
  phase1 = -0.5
  gap1 = 2*move_parms.width(mp_fill)
  ph1 = circle(ctr, R1, nt1, phase1, mp_fill, gap1)
  path.set_name(ph1, "P1", True)
  
  OPHS = [ ph0, ph1, ]

  for oph in OPHS: path.validate(oph)
 
  return OPHS
  # ----------------------------------------------------------------------

def misc_K(nmv,step,pert,mp_cont,mp_fill,mp_link,mp_jump):

  # Build a square with a perturbed small moves:
  p0 = (1,1) # Initial point.
  u = (1,0) # Direction of initial side.
  random.seed(4615+417)
  MVS = []
  p = p0 # Previous unperturbed point.
  pp = p # Previous perturbed point.
  for ksd in range(4):
    # Append side {ksd}:
    tagk = "%d" % ksd
    v = (-u[1], +u[0])
    for imv in range(nmv):
      q = rn.mix(1.0, p, step, u) # Unperturbed next point.
      # Compute the perturbed point {qq}:
      if imv < nmv-1:
        du = pert*(2*random.random() - 1)
        dv = pert*(2*random.random() - 1)
        qq = rn.mix3(1.0,q, du, u, dv, v) # Perturbed point
      else:
        qq = q
      if ksd == 0:
        # Side 0 -- all moves are {mp_fill}:
        mp = mp_fill
      elif ksd == 1:
        # Side 1: half the moves are {mp_fill}, half are {mp_cont}
        mp = mp_fill if imv < nmv//2 else mp_cont
      elif ksd == 2:
        # Side 2: moves are {mp_fill} except one is {m_jump}:
        mp = mp_fill if imv != nmv//2 else mp_jump
      elif ksd == 3:
        # Side 3: all moves are jumps with different parameters:
        mp = mp_jump if imv == nmv-1 else move_parms.make(0, 1000+imv, 200, 0.2+0.01*imv)
      else:
        assert False
      path.move_to(MVS, pp, qq, mp, tagk)
      p = q; pp = qq
    u = v
  oph = path.from_moves(MVS)
  path.set_name(oph, "P", False)
  
  # Add some links:
  qini = path.pini(oph); eini = rn.add(qini, (-2, 0)); gini = rn.add(qini, (-3, 1))
  
  path.add_link(oph, path.from_points((eini, qini), mp_link, None))
  path.add_link(oph, path.from_points((gini, qini), mp_link, None))
  
  qfin = path.pfin(oph); efin = rn.add(qini, (0, -2)); gfin = rn.add(qini, (1, -3))

  path.add_link(path.rev(oph), path.from_points((efin, qfin), mp_link, None))
  path.add_link(path.rev(oph), path.from_points((gfin, qfin), mp_link, None))
  
  # Set a random index group:
  path.set_group(oph, 418)
  
  return oph
  # ----------------------------------------------------------------------

# MISCELLANEOUS CONTOUR PATHS FOR TESTING

def contours_A(mp_cont):
  
  pA0 = ( 1, 1)
  pA1 = (12, 1)
  pA2 = (12,10)
  pA3 = ( 1,10)
  ptsA = (pA0,pA1,pA2,pA3,)
  
  pC0 = ( 2, 2)
  pC1 = ( 5, 2)
  pC2 = ( 5, 9)
  pC3 = ( 2, 9)
  ptsC = (pC3,pC2,pC1,pC0,)
  
  pF0 = ( 3, 3)
  pF1 = ( 4, 3)
  pF2 = ( 4, 5)
  pF3 = ( 3, 5)
  ptsF = (pF0,pF1,pF2,pF3,)
  
  pG0 = ( 3, 6)
  pG1 = ( 4, 6)
  pG2 = ( 4, 8)
  pG3 = ( 3, 8)
  ptsG = (pG0,pG1,pG2,pG3,)
  
  pD0 = ( 6, 2)
  pD1 = (11, 2)
  pD2 = (11, 9)
  pD3 = ( 6, 9)
  ptsD = (pD3,pD2,pD1,pD0,)
  
  pH0 = ( 7, 3)
  pH1 = (10, 3)
  pH2 = (10, 8)
  pH3 = ( 7, 8)
  ptsH = (pH0,pH1,pH2,pH3,)
  
  pI0 = ( 8, 4)
  pI1 = ( 9, 4)
  pI2 = ( 9, 5)
  pI3 = ( 8, 5)
  ptsI = (pI3,pI2,pI1,pI0,)
  
  pJ0 = ( 8, 6)
  pJ1 = ( 9, 6)
  pJ2 = ( 9, 7)
  pJ3 = ( 8, 7)
  ptsJ = (pJ3,pJ2,pJ1,pJ0,)
  
  pB0 = ( 3, 11)
  pB1 = ( 7, 11)
  pB2 = ( 7, 15)
  pB3 = ( 2, 15)
  ptsB = (pB0,pB1,pB2,pB3,)
  
  pE0 = ( 5, 12)
  pE1 = ( 6, 13)
  pE2 = ( 5, 14)
  pE3 = ( 4, 13)
  ptsE = (pE3,pE2,pE1,pE0,)
  
  PTSS = ( ptsA, ptsB, ptsC, ptsD, ptsE, ptsF, ptsG, ptsH, ptsI, ptsJ, )
  
  CRS = []
  for ipo in range(len(PTSS)):
    PTS = PTSS[ipo]
    cr = path.from_points(PTS + (PTS[0],), mp_cont, None)
    path.set_name(cr, ("Q%d" % ipo), True)
    CRS.append(cr)
  
  path.compute_contour_nesting(CRS)

  for cr in CRS: path.validate(cr)
   
  return CRS, PTSS
  # ----------------------------------------------------------------------

def contours_B(mp_cont):
  ctr = (10,10)
  nt0 = 10; R0 = 5
  nt1 = 20; R1 = 7

  phc0 = path_example.circle(ctr, R0, nt0, -pi/3, mp_cont, 0.0)
  path.set_name(phc0, "P0", True)

  phc1 = path_example.circle(ctr, R1, nt1, +pi/3, mp_cont, 0.0)
  path.set_name(phc1, "P1", True)

  CRS = [ phc0, phc1 ]

  for cr in CRS: path.validate(cr)
   
  path.compute_contour_nesting(CRS)

  return CRS
  # ----------------------------------------------------------------------

# PATH COMPONENTS

def circle(ctr, R, nt, phase, mp_trace, gap):

  wd = move_parms.width(mp_trace)
  mp_jump = None # Should be just traces.

  # Leave a small gap at the end for visual clarity.
  # Also it may avoid a bump on the surface.
  pts = []
  astep = 2*pi/nt
  agap = gap/2/R # Angular size of gap between the endpoints.
  aini = phase + agap/2        # Start extruding at this angle.
  afin = phase + 2*pi - agap/2 # Stop extruding at this angle.
  aant = None # Previous value of angle {a}.
  pant = None # Corresponding vertex of unclipped polygon.
  eps = 0.01*astep
  for kt in range(nt+1):
    a = phase + kt*astep
    ared = phase + (kt % nt)*astep # To ensure last {p} same as first {p}.
    p = rn.add(ctr, (R*cos(ared), R*sin(ared)))
    if gap == 0:
      # Closed circle:
      pts.append(p)
    else:
      # Circle with a gap:
      if kt == 0:
        # Starts inside the gap.  Omit:
        pass
      elif aant < aini - eps and a > aini + eps:
        # Trace crosses out of the gap.
        # Add only the final part of the trace:
        r = (aini - aant)/astep # Not accurate, but will do.
        q = rn.mix(1-r, pant, r, p)
        pts.append(q)
        # sys.stderr.write("(0) a = %20.16f\n" % aini)
        pts.append(p)
        # sys.stderr.write("(1) a = %20.16f\n" % ared)
      elif aant < afin - eps and a > afin + eps:
        # Trace crosses into the gap.
        # Add only the initial part of the trace:
        r = (afin - aant)/astep # Not accurate, but will do.
        q = rn.mix(1-r, pant, r, p)
        pts.append(q)
        # sys.stderr.write("(3) a = %20.16f\n" % afin)
      elif a >= aini and a <= afin:
        # Add vertexto path:
        pts.append(p)
        # sys.stderr.write("(4) a = %20.16f\n" % ared)
      else:
        # Omit vertex and previous side, if any:
        pass
    pant = p
    aant = a
  ph = path.from_points(pts, mp_trace, mp_jump)
  return ph

def gear(ctr, Rin, Rex, nt, split, phase, mp_trace, mp_jump):

  wd = move_parms.width(mp_trace)

  ptss = [] # List of lists of points, each being a separate subpath.
  pts = None # Points of current subpath, or {None} between subpaths.
  # Incomplete gear:
  for kt in range(2*nt):
    # Generate the half-tooth number {kt}:
    if kt % 2 == 0:
      R0 = Rin; R1 = Rex
    else:
      R0 = Rex; R1 = Rin
    # Compute points of gear half-tooth:
    a = (kt+0.00)*pi/nt + phase
    p1 = rn.add(ctr, (R0*cos(a), R0*sin(a)))
    b = (kt+0.50)*pi/nt + phase
    p2 = rn.add(ctr, (R0*cos(b), R0*sin(b)))
    p3 = rn.add(ctr, (R1*cos(b), R1*sin(b)))
    c = (kt+1.00)*pi/nt + phase
    p4 = rn.add(ctr, (R1*cos(c), R1*sin(c)))
    if split and (kt == 0 or kt == 1 or kt == 4 or kt == 5):
      # Omit the half-tooth, ending the current path if any::
      if pts != None:
        ptss.append(pts)
        pts = None
    else:
      # Add the half-tooth:
      if pts == None:
        # First point of subpath:
        pts = [ p1 ]
      else:
        # Must be joined:
        assert p1 == pts[-1]
      pts.append(p2)
      pts.append(p3)
      pts.append(p4)
  if pts != None: ptss.append(pts)
  ph = path.from_points(ptss, mp_trace, mp_jump)
  return ph

# FILLING COMPONENTS

def circle_rootray(p, q, R):
  # Returns a root-ray for the straight line path {r(t) = p + t (q-p)}
  # relative to a circle with center at the origin and radius {R}.

  if R <= 0: return rootray_vacuum()
  # Corresponding points for the unit ball:
  pu = rn.scale(1/R, p)
  qu = rn.scale(1/R, q)
  K = rootray_cart.ball(pu,qu)
  # sys.stderr.write("R = %10.6f y = %10.6f H = %10.6f K = %s\n" % (R, y, H, str(K)))
  return K

def ring_clip_segment(Rin, Rex, p, q, side):
  # Clips the segment {p--q} to the ring with center at the origin,
  # inner radius {Rin} and outer radius {Rex}. The result is a tuple with
  # zero, one, or two pairs of points. Each pair is the endpoints of one
  # part of the clipped segment.
  #
  # When this procedure is used for filling, the caller must adjust the
  # radii to account for the thickness of contours and filling traces.
  # sys.stderr.write("\n")
  Kot = circle_rootray(p, q, Rex)
  Kin = circle_rootray(p, q, Rin)
  Ksg = (+1, (0, 1,),)
  K = rootray.intersection(rootray.difference(Kot, Kin), Ksg)
  if K == None:
    # Scan line does not meet the ring:
    return []
  assert type(K) is list or type(K) is tuple and len(K) == 2
  sK,tK = K # Starting status and root list.
  # sys.stderr.write("tK = %s\n" % str(tK))
  assert sK == +1
  nt = len(tK)
  # Decide which pair(s) of roots to use. The low indices are {ka,kb}..
  if nt == 0:
    # Segment is entirely outside the ring.
    return []
  elif nt == 2:
    # Segment crosses the outer circle but not the inner one:
    ka = 0; kb = None
  elif nt == 4:
    # Segment crosses both circles:
    if side < 0:
      ka = 0; kb = None
    elif side > 0:
      ka = 2; kb = None
    else:
      ka = 0; kb = 2
  else:
    assert False, "inconsistent number of roots"
  assert ka != None

  # Now generate the point pairs list:
  segs = [] # List of pairs of points.
  for k in ka, kb:
    if k != None:
      ptpair = []
      for ipt in k,k+1:
        t = tK[ipt]
        pt = rn.mix(1-t, p, t, q)
        ptpair.append(pt)
      segs.append(ptpair)
  return segs

def ring_raster(Rin, Rex, y, side, mp_fill, mp_jump):

  wdf = move_parms.width(mp_fill)

  # Adjust the radii accounting for the round caps of the raster lines:
  Sin = Rin if Rin < 0 else Rin + wdf/2
  Sot = Rex - wdf/2

  # Quick tests for emptiness:
  if Sot <= 0: return None
  if Sin > 0 and Sot <= Sin: return None
  if abs(y) >= Sot: return None

  # Choose {p,q} that define the ray:
  h = 1.1*Rex    # Sufficient {X} distance to be outside the ring.
  p = (-Rex, y)  # On the {Y}-axis.
  q = (+Rex, y)  # On the same scan-line, far away from the axis.

  # Find the ray intersection with the ring,
  segs = ring_clip_segment(Sin, Sot, p, q, side)
  if len(segs) == 0: return None
  return path.from_points(segs, mp_fill, mp_jump)

def ring_append_clipped_seg(pts, Rin, Rex, p, q):
  # The parameter {pts} must be a list of lists of points in the
  # format suitable for {path.from_points}.
  # Clips the segment {p--q} to the ring  with center at the origin,
  # inner radius {Rin} and outer radius {Rex}, and appends
  # the pieces to {pts}, with breaks at the gaps.
  assert type(pts) is list
  # Get last chain {cant} of {pts} and its final point {pant}:
  if len(pts) == 0:
    cant = None; pant = None
  else:
    cant = pts[-1] # Last chain in {pts}.
    assert len(cant) >= 2
    pant = cant[-1]
  # Clip {p--q} to the ring:
  # sys.stderr.write("p = %s q = %s\n" % (str(p), str(q)))
  segs = ring_clip_segment(Rin, Rex, p, q, 0)
  # sys.stderr.write("segs = %s\n" % str(segs))
  assert type(segs) is list or type(segs) is tuple
  # Now {segs} should be a list of pairs of points.
  for seg in segs:
    # Get first point of next segment {seg}:
    assert type(seg) is list or type(seg) is tuple
    assert len(seg) == 2
    ppos = seg[0]
    if pant == None or ppos != pant:
      # Start a new chain:
      cnew = [ seg[0], seg[1], ]
      pts.append(cnew)
      cant = cnew
    else:
      # Append point to last chain:
      cant.append(seg[1])
    pant = seg[1]

def ring_clip_polyline(pts_unc, Rin, Rex):
  # Given a list {pts} of vertices of a polyline,
  # returns its intersection with the ring with center at
  # origin, inner radius {Rin}, outer radius {Rex}.
  #
  # The output is a list of lists of points suitable
  # for {path;from_points}.
  assert type(pts_unc) is list or type (pts) is tuple
  # sys.stderr.write("pts_unc = %s\n" % str(pts_unc))
  pts_clp = []
  pant = None
  for p in pts_unc:
    assert hacks.is_point(p)
    if pant != None:
      ring_append_clipped_seg(pts_clp, Rin, Rex, pant, p)
    pant = p
  # sys.stderr.write("pts_clp = %s\n" % str(pts_clp))
  return pts_clp

def zigzag_segs(w, y, h, ytrim, phase):
  # Returns a list of points that are the vertices of the axes
  # of a zigzag path with trimmed 60 degree teeth. The path
  # is not clipped against anything.
  #
  # The medial line of the zigzag is the scanline at height {y}.
  # The tips of the untrimmed teeth lie at ordinates {y+h} and
  # {y-h}.  The depth of trimming is {ytrim}.
  # The polygonal line extends from at least {-w} to {+w} in {X}.
  #
  # There is the tip of a tooth at (phase*stride, y + h) where {stride} is the
  # horizontal pitch (distance between successive tooth tips).
  assert ytrim < h

  # Some geometric parameters:
  stride = 4*h/sqrt(3) # Distance between teeth on the same side of {L}.
  xtrim = ytrim/sqrt(3) # Half-width of flat part of tooth.

  # Start on {L}, on the {-X} side, well away from the circle:
  n = int(ceil(w/stride + 0.00001)) + 1 # Sufficient number of teeth.
  assert n >= 1
  phase = (phase % 1.0) # Reduce the phase to {[0 _ 1)}
  assert phase >= 0 and phase < 1.0
  xprev = (- n + phase)*stride # Ideal toot point abscissa.

  sprev = +1 # Tooth tip direction:  {+1} or {-1} for above or below {L}.
  ymin = y - h; ymax = y + h # {Y} range for points on the zigzag axis

  pts = []    # List of points.

  # Now generate the teeth:
  while xprev <= +w:
    # The point {(xprev,y+sprev*h} is the previous ideal tooth tip.
    # We generate the sloping part of the tooth after {xprev},
    # plus the flat part at the top of the next toot.

    # Compute the position {(xnext,y+snext*h)} of the next ideal tooth tip:
    xnext = xprev + stride/2
    snext = -sprev

    # Get endpoints of sloping part:
    p = (xprev + xtrim, y + sprev*(h - ytrim))
    pts.append(p)

    q = (xnext - xtrim, y + snext*(h - ytrim))
    pts.append(q)

    # Advance to next tooth:
    xprev = xnext
    sprev = snext

  assert len(pts) >= 2
  return pts

def ring_zigzag(Rin, Rex, y, step, phase, mp_fill, mp_jump):

  wdf = move_parms.width(mp_fill)
  assert step > wdf

  # Adjust radii to account for trace round caps:
  Sin = Rin if Rin < 0 else Rin + wdf/2
  Sot = Rex - wdf/2

  # Quick tests for emptiness:
  if Sot <= 0: return None
  if Sin > 0 and Sot <= Sin: return None
  if abs(y) >= Sot: return None

  # Some geometric parameters:
  height = step       # Displacement of ideal tooth tips from midline.
  ytrim = wdf # Height of chip removed from tooth tips.

  pts_unc = zigzag_segs(Rex, y, height, ytrim, phase)
  pts_clp = ring_clip_polyline(pts_unc, Sin, Sot)
  if len(pts_clp) == 0:
    ph = None
  else:
    ph = path.from_points(pts_clp, mp_fill, mp_jump)
  return ph

def ring_fill(Rin, Rex, dmin, dmax, step, zigzag, side, mp_fill, mp_jump):

  if zigzag: side = 0 # Ignore {side} in zigzag fill.

  # Find range of {kph}:
  kphmin = int(floor(dmin/step - 0.0001))-1 # Start below {dmin}, just to be sure.
  kphmax = int(ceil(dmax/step + 0.0001))+1  # Stop above {dmax}, just to be sure
  kph = kphmin
  OPHS = [ ] # List of paths
  while (kph <= kphmax):
    y = kph*step
    if dmin <= y and y <= dmax:
      if not zigzag or (kph % 2) == 0:
        phk = ring_raster(Rin, Rex, y, side, mp_fill, mp_jump)
      else:
        phase = ((kph + 3) % 4)/4
        assert phase == 0.0 or phase == 0.5
        phk = ring_zigzag(Rin, Rex, y, step, phase, mp_fill, mp_jump)
      if phk != None:
        ophk = phk if (kph % 2) == 0 else path.rev(phk)
        # sys.stderr.write("ophk =%s\n" % str(ophk))
        OPHS.append(ophk)
    kph = kph + 1
  if len(OPHS) == 0:
    sys.stderr.write("path would be trivial!\n")
    ph = None
  else:
    # Concatenate the elements, using links if appropriate:
    use_links = True
    ph = path.concat(OPHS, use_links, mp_jump)
  return ph

