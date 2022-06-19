#! /usr/bin/python3
# Test program for module {move}

import move
import move_example
import move_parms
import hacks
import palette
import job_parms
import rn
import pyx
import sys
from math import sqrt, sin, cos, floor, ceil, inf, nan, pi

parms = job_parms.slow()

mp_fill = move_parms.make_for_fillings(parms)
mp_cont = move_parms.make_for_contours(parms)
mp_jump = move_parms.make_for_jumps(parms)

wd_fill = move_parms.width(mp_fill)
wd_cont = move_parms.width(mp_cont)

def test_make_width_is_jump():
  sys.stderr.write("--- testing {make,width,is_jump,is_trace} ---\n")

  p1 = (1,1)
  p2 = (3,1)
  p3 = (2,3)
  p4 = (4,4)
  p5 = (1,4)

  assert not move.is_trace(None); 
  assert not move.is_jump(None); 
  
  mv12 = move.make(p1,p2, mp_fill); 
  assert move.is_trace(mv12); 
  assert not move.is_jump(mv12); 
  assert move.width(mv12) == wd_fill

  mv23 = move.make(p2,p3,   mp_jump); 
  assert not move.is_trace(mv23);     
  assert move.is_jump(mv23);     
  assert move.width(mv23) == 0

  mv34 = move.make(p3,p4, mp_cont); 
  assert move.is_trace(mv34); 
  assert not move.is_jump(mv34); 
  assert move.width(mv34) == wd_cont

  mv55 = move.make(p5,p5, mp_fill); 
  assert move.is_trace(mv55); 
  assert not move.is_jump(mv55); 
  assert move.width(mv55) == wd_fill

def test_move_orient_unpack_rev():
  sys.stderr.write("--- testing {spin,unpack,rev,pini,pfin,endpoints,length} ---\n")

  p1 = (1,1)
  p2 = (3,1)

  mv12 = move.make(p1,p2, mp_fill); 

  def check_pini_pfin(omvx):
    mvx, drx = move.unpack(omvx)
    assert mvx == mv12
    if drx == 0:
      assert move.pini(omvx) == p1
      assert move.pfin(omvx) == p2
      assert move.endpoints(omvx) == (p1, p2)
    else:
      assert move.pini(omvx) == p2
      assert move.pfin(omvx) == p1
      assert move.endpoints(omvx) == (p2, p1)
    assert move.length(omvx) == rn.dist(p1,p2)
  
  for dr in range(4):
    omva = move.spin(mv12, dr)
    mva, dra = move.unpack(omva)
    assert dra == dr % 2
    check_pini_pfin(omva)
    
    omvb = move.spin(move.rev(mv12), dr)
    mvb, drb = move.unpack(omvb)
    assert drb == (dr + 1) % 2
    check_pini_pfin(omvb)

def test_displace():
  sys.stderr.write("--- testing {displace} ---\n")

  p1 = (1,1)
  p2 = (3,1)

  mv12 = move.make(p1,p2, mp_fill); 

  angr = pi/6
  vr = (2,3)
  mv12r = move.displace(mv12, angr, vr, mp_fill)
  p1r = rn.add(rn.rotate2(p1, angr), vr)
  p2r = rn.add(rn.rotate2(p2, angr), vr)
  assert rn,dist(move.pini(mv12r), p1r) < 1.0e-8
  assert rn,dist(move.pfin(mv12r), p2r) < 1.0e-8

def test_plot():

  sys.stderr.write("--- testing {plot_standard,plot_layer,bbox} ---\n")

  OMVS = move_example.misc_A(mp_fill, mp_cont, mp_jump)
  nmv = len(OMVS)
  CLRS = hacks.trace_colors(nmv, None)
  nclr = len(CLRS)
 
  # Testing all combinations of options:
  
  # Compute the plot's bounding box:
  B = move.bbox(OMVS)

  dp = (0, 0)

  c, szx,szy = hacks.make_canvas(hacks.round_box(B,0.5), dp, None, True, True, 4, 8)

  # Colors:
  cmatter = pyx.color.rgb(0.750, 0.200, 0.600) # Est. material footprint.
  ctraces = pyx.color.rgb(0.850, 0.050, 0.000)  # Sausages of traces.
  caxes = pyx.color.rgb(0.450, 0.050, 0.000)   # Axis lines of traces.
  cdots = caxes                                # End dots of traces.
  cjumps = pyx.color.rgb.black                  # Axis lines, dots, and arrowheads of jumps.

  # Relative to the nominal width:
  rwd_matter = 1.13; 
  rwd =  0.80; 

  # Absolute (in mm):
  wd_ref = 0.15*min(wd_fill,wd_cont) # Reference line width.
  wd_axes = wd_ref
  wd_dots = 2.5*wd_ref
  wd_jump = wd_ref
  sz_arrows = 6*wd_ref
  wd_bbox = 0.25*wd_ref

  def test_plot_do_bbox(omv, dp, wd_bbox):
    # Plots the bounding box of move {omv} displaced by {dp}.
    B = move.bbox([omv,])
    hacks.plot_frame(c, B, dp, pyx.color.rgb.blue, wd_bbox, None, -wd_bbox/2)
    return
    # ----------------------------------------------------------------------

  def test_plot_jump(dp):
    # Plot a jump that was giving problems in one program:
    p0 = (1,1)
    p1 = (5,2)
    jmp = move.make(p0, p1, mp_jump)
    clr = pyx.color.rgb( 1.000, 0.000, 1.000 )
    wd_axis = 0.10*wd_fill
    dashed = True
    wd_dots = 2.5*wd_axis
    sz_arrow = 6*wd_axis
    move.plot_layer(c, jmp, dp, clr, wd_axis, dashed, wd_dots, sz_arrow)
    return
    # ----------------------------------------------------------------------

  def test_plot_do_layer(dp, dr, axes, dots, arrows):
    # Tests {move.plot_layer} for one combination of direction and options. 
    # sys.stderr.write("\n")
    # sys.stderr.write("  {test_plot_do_layer}: dr = %d axes = %d dots = %d arrows = %d\n" % (dr,axes,dots,arrows))

    def pick_colors(kmv):
      # Returns the colors for trace sausages and axes of move {OMVS[kmv]}.
      if nclr == 1:
        ctrace = CLRS[0]
      else:
        ctrace = CLRS[kmv]
      caxis =   pyx.color.rgb(0.6*ctrace.r, 0.6*ctrace.g, 0.6*ctrace.b) # Color of trace axis, dots, arrow.
      return ctrace, caxis

    OMVSdr = [ move.spin(omv, dr) for omv in OMVS ]

    wd_axes_f = wd_axes if axes else 0 
    wd_dots_f = wd_dots if dots else 0
    sz_arrows_f = sz_arrows if arrows else 0

    tp = (dp[0]+0.5, dp[1]+szy-0.5,)
    c.text(tp[0], tp[1], "BYL dr:%d ax:%d dt:%d ar:%d" % (dr,int(axes),int(dots),int(arrows)))
    for layer in range(4):
      # sys.stderr.write("  ... layer %d ...\n" % layer)
      for kmv in range(nmv):
        omv = OMVSdr[kmv]
        # sys.stderr.write("\n")
        jmp = move.is_jump(omv)
        wd = move.width(omv)
        if layer == 0 and not jmp:
          # plots the estimate of actual material:
          wd_matter = rwd_matter*wd
          move.plot_layer \
            (c, omv, dp, clr=cmatter, wd=wd_matter, dashed=False, wd_dots=0, sz_arrow=None)
        elif layer == 1 and not jmp:
          # plots the nominal trace material:
          wd_trace = rwd*wd
          ctrace, caxis = pick_colors(kmv)
          move.plot_layer \
            (c, omv, dp, clr=ctrace, wd=wd_trace, dashed=False, wd_dots=0, sz_arrow=None)
        elif layer == 2 and not jmp:
          # Plots the trace axis:
          ctrace, caxis = pick_colors(kmv)
          move.plot_layer \
            (c, omv, dp, clr=caxis, wd=wd_axes_f, dashed=False, wd_dots=wd_dots_f, sz_arrow=sz_arrows_f)
        elif layer == 3 and jmp:
          # Plots the jump:
          move.plot_layer \
            (c, omv, dp, clr=cjumps, wd=wd_jump, dashed=True, wd_dots=wd_dots_f, sz_arrow=sz_arrows_f)
        if layer == 0:
          test_plot_do_bbox(omv, dp, wd_bbox)
    return
    # ----------------------------------------------------------------------

  def test_plot_do_standard(dp, dr, axes, dots, arrows):
    # Tests {move.plot_standard} for one combination of direction and options. 
    # sys.stderr.write("\n")
    # sys.stderr.write("  {test_plot_do_standard}: dr = %d axes = %d dots = %d arrows = %d\n" % (dr,axes,dots,arrows))

    OMVSdr = [ move.spin(omv, dr) for omv in OMVS ]

    tp = (dp[0]+0.5, dp[1]+szy-0.5,)
    c.text(tp[0], tp[1], "STD dr:%d ax:%d dt:%d ar:%d" % (dr,int(axes),int(dots),int(arrows)))
    if dr == 0:
      # Plot all layers for each move. Note that {mv23} is a jump:
      move.plot_standard \
        ( c, OMVSdr, dp, None, CLRS, rwd=rwd, wd_axes=wd_axes, 
          axes=axes, dots=dots, arrows=arrows, matter=True
        ) 
      for omv in OMVSdr:
        test_plot_do_bbox(omv, dp, wd_bbox)
    else:
      # Plot layer by layer:
      for layer in range(4):
        # sys.stderr.write("  ... layer %d ...\n" % layer)
        move.plot_standard \
          ( c, OMVSdr, dp, layer, CLRS, rwd=rwd, wd_axes=wd_axes, 
            axes=axes, dots=dots, arrows=arrows, matter=True
          )
    return
    # ----------------------------------------------------------------------
    
  iy = 0  # Plot row index.
  for axes in False,True:
    for dots in False,True:
      for arrows in False,True:
        for dr in range(2):
          dpk = rn.add(dp, ((2*dr+0)*szx, iy*szy))
          test_plot_do_layer(dpk, dr, axes, dots, arrows)
          dpk = rn.add(dp, ((2*dr+1)*szx, iy*szy))
          test_plot_do_standard(dpk, dr, axes, dots, arrows)
        iy += 1

  dpn = rn.add(dp, (0, iy*szy))
  test_plot_jump(dpn)

  hacks.write_plot(c, "tests/out/move_TST_plot")
  return
  # ----------------------------------------------------------------------

def test_fabtime():
  sys.stderr.write("--- testing {fabtime,ud_penalty,transition_penalty} ---\n")

  p1 = (1,1)
  p2 = (3,1)
  p3 = (2,3) 
  p4 = (4,4)
  p5 = (1,4)

  mv12 = move.make(p1,p2, mp_cont); 
  mv23 = move.make(p2,p3, mp_jump); 
  mv34 = move.make(p3,p4, mp_cont); 
  mv45 = move.make(p4,p5, mp_fill); 
  
  ud12 = move.ud_penalty(mv12); assert ud12 == 0
  ud23 = move.ud_penalty(mv23); assert ud23 == parms['extrusion_on_off_time']
  ud34 = move.ud_penalty(mv34); assert ud34 == 0
  ud45 = move.ud_penalty(mv45); assert ud45 == 0

  tp123 = move.transition_penalty(mv12, mv23); assert tp123 == ud23
  tp234 = move.transition_penalty(mv23, mv34); assert tp234 == ud23
  tp345 = move.transition_penalty(mv34, mv45); assert tp345 == 0
  
  assert move.transition_penalty(mv12, None) == 0
  assert move.transition_penalty(mv23, None) == 0
  
  assert move.transition_penalty(None, mv12) == 0
  assert move.transition_penalty(None, mv23) == 0
  
  assert move.transition_penalty(None, None) == 0

  d23 = rn.dist(p2,p3)
  t23a = move_parms.nozzle_travel_time(d23, None, mp_jump)
  t23b = move.fabtime(mv23)
  assert t23a == t23b

  u34 = rn.sub(p4,p3)
  v34 = (-u34[1], +u34[0]) # A vector perpendicular to the segment {p3--p4}.
  d34 = rn.dist(p3,p4)
  r34a = 0.75
  p34a = rn.mix(r34a, p3, 1-r34a, p4) # A point on the segment {p3--p4}, 3/4 of the way.
  q34a = rn.mix(0.01, v34, 1.0, p34a) # Displace p34a off the segment a bit.
  t34a = move_parms.nozzle_travel_time(d34, r34a*d34, mp_cont)
  t34b = move_parms.nozzle_travel_time(d34, (1-r34a)*d34, mp_cont)
  assert abs(t34a + t34b - move.fabtime(mv34)) < 1.0e-8

def test_cover_time():
  sys.stderr.write("--- testing {cover_time,plot_to_files} ---\n")

  # Make two moves:

  p11 = (1,1)
  p12 = (6,1)
  mv1 = move.make(p11, p12, mp_fill)
  
  p21 = (3, 1-wd_fill)
  p22 = (7, 1-wd_fill)
  mv2 = move.make(p21, p22, mp_fill)

  # Pick a contact point {ctm} between them:
  v0 = (4, 1-wd_fill/2) # Start of contact.
  v1 = (6, 1-wd_fill/2) # End of contact.
  ctm = rn.mix(0.5, v0, 0.5, v1) # Midpoint of contact.

  OMVS = [mv1, mv2,]
  rwd = 0.80
  wd_axes = 0.05*wd_fill
  ctraces = pyx.color.rgb(1.000, 0.200, 1.000)
  move.plot_to_files("tests/out/move_TST_cover_time", OMVS, [ctraces,], rwd, wd_axes)

  # Check contact cover time on moves:
  
  for omv in mv1, mv2, move.rev(mv1), move.rev(mv2):
    mp = move.parameters(omv)
    p, q = move.endpoints(omv)
    dpq = rn.dist(p, q)
    rmq = abs(ctm[0] - p[0])/abs(q[0] - p[0]) # Rel pos of contact on {mv1}
    tca = move_parms.nozzle_travel_time(dpq, rmq*dpq, mp)
    tcb = move.cover_time(omv, ctm)
    tcc = move.cover_time(move.rev(omv), ctm)
    assert abs(tca - tcb) < 1.0e-8
    assert abs((tca + tcc) - move.fabtime(omv)) < 1.0e-8

  return
  # ----------------------------------------------------------------------

def test_shared_border():
  sys.stderr.write("--- testing {shared_border} ---\n");

  # Approximate line of the contact:
  e0 = (3,1)
  e1 = (7,3)

  xdir, um = rn.dir(rn.sub(e1, e0))  # Unit vector of {e0-->e1}.
  ydir = (-xdir[1],+xdir[0])         # Unit vector orthogonal to {xdir}.
  
  xe0 = rn.dot(e0, xdir)
  xe1 = rn.dot(e1, xdir)

  mp0 = mp_fill; wd0 = move_parms.width(mp0)
  mp1 = mp_cont; wd1 = move_parms.width(mp1)

  eps = 0.01 # Relative perturbation size
  tol = 0.20*max(wd0,wd1)  # Tolerance for overlaps, tilts, etc.
  
  ysep = (wd0+wd1)/2 # Ideal separation.

  def do_test_shared_border(r0, r1, rsep, res_exp):

    sys.stderr.write("testing with r0 = %+.7f r1 = %+.7f rsep = %.7f expects %s\n" % (r0, r1, rsep, res_exp))

    # Endpoints of first move:
    p0 = rn.mix(1.0,e0, -0.5*wd0*rsep + 1.0*eps*ysep, ydir)
    p1 = rn.mix(1.0,e1, -0.5*wd0*rsep - 2.0*eps*ysep, ydir)

    # Endpoints of second move:
    q0 = rn.mix3(1-r0, e0, r0, e1, +0.5*wd1*rsep - 2.0*eps*ysep, ydir)
    q1 = rn.mix3(1-r1, e0, r1, e1, +0.5*wd1*rsep + 2.0*eps*ysep, ydir)

    sep = rsep*ysep # Tests separation.

    # Tests with moves parallel and antiparallel:
    for kdr in range(2):
      mv0 = move.make(p0, p1, mp0)
      mv1 = move.make(q0, q1, mp1)
      c0, c1 = move.shared_border(mv0, mv1, xdir, tol)
      assert (c0 == None) == (c1 == None)
      if c0 != None:
        xc0 = rn.dot(c0, xdir); rc0_cmp = (xc0 - xe0)/(xe1 - xe0)
        xc1 = rn.dot(c1, xdir); rc1_cmp = (xc1 - xe0)/(xe1 - xe0)
        sys.stderr.write("c0  = ( %20.16f %20.16f ) @ %20.16f" % (c0[0],c0[1],rc0_cmp))
        sys.stderr.write("  c1  = ( %20.16f %20.16f ) @ %20.16f\n" % (c1[0],c1[1],rc1_cmp))
        assert ok_exp, "expected {None,None}"
      else:
        assert not ok_exp, "expected success"
      # Test indifference to orientation and order:
      for omva, omvb in (mv0, move.rev(mv1)), (mv1,mv0), (mv0, move.rev(mv1)), (mv1, move.rev(mv0)): 
        c0r,c1r = move.shared_border(omva, omvb, xdir, tol)
        if c0r != None and c1r != None:
          sys.stderr.write("c0r = ( %20.16f %20.16f ) c1r = ( %20.16f %20.16f )\n" % (c0r[0],c0r[1],c1r[0],c1r[1]))
          assert \
            ( hacks.same_point(c0r, c0, 1.0e-8) and hacks.same_point(c1r, c1, 1.0e-8) ) or \
            ( hacks.same_point(c0r, c1, 1.0e-8) and hacks.same_point(c1r, c0, 1.0e-8) )
        else:
          assert  c0r == None and c1r == None
      # Flip segment {q0--q1}:
      q0,q1 = q1,q0

    # Tests with vectors not parallel:
    q0 = rn.mix3(1-r0, p0, r0, p1, 0.5, ydir)
    q1 = rn.mix3(1-r1, p0, r1, p1, 0.1, ydir)
    for kdr in range(2):
      mv1 = move.make(q0, q1, mp1)
      c0, c1 = move.shared_border(mv0, mv1, xdir, tol)
      assert c0 == None and c1 == None
      # Flip segment {q0--q1}:
      q0,q1 = q1,q0

    return
    # ....................................................................

  for r0,r1 in (
      (-0.04, -0.02), 
      (-0.02, -0.02), 
      (-0.02, +0.02), 
      (-0.02, +0.50), 
      (-0.02, +0.98), 
      (-0.02, +1.02), 
      (+0.02, +0.02), 
      (+0.02, +0.04), 
      (+0.02, +0.50), 
      (+0.02, +0.98), 
      (+0.02, +1.02), 
      (+0.50, +0.50), 
      (+0.50, +0.98), 
      (+0.50, +1.02), 
      (+0.96, +0.98), 
      (+0.98, +0.98), 
      (+0.98, +1.02), 
      (+1.02, +1.02), 
      (+1.02, +1.04), 
    ):
    assert r0 <= r1
    for rsep in (0.50, 0.99, 1.01, 2.0):
      ok_exp = not (r1 <= 0 or r0 >= 1 or r1 - r0 < 0.001 or rsep < 0.8 or rsep > 1.2)
      res_exp = "success" if ok_exp else "failure"
      do_test_shared_border(r0, r1, rsep, res_exp)

  return
  # ----------------------------------------------------------------------

def test_name():
  sys.stderr.write("--- testing {has_name,set_name,get_name,tag_names} ---\n");
  
  def gname(var, omv):
    # Returns {get_name(omv)}, printing to {stderr}.
    xms = move.get_name(omv)
    sys.stderr.write("  get_name(%s) = %s\n" % (var,str(xms)))
    return xms
    # ............................................................

  tra = move.make((1,1), (3,2), mp_fill)
  assert not move.has_name(tra)
  xms = gname("tra", tra)
  assert xms == "T?"
  xms = gname("~tra", move.rev(tra))
  assert xms == "~T?"

  move.set_name(tra, "BB")
  assert move.has_name(tra)
  xms = gname("tra", tra)
  assert xms == "BB"
  xms = gname("~tra", move.rev(tra))
  assert xms == "~BB"

  jma = move.make((1,1), (3,2), mp_jump)
  assert not move.has_name(jma)
  xms = gname("jma", jma)
  assert xms == "J?"
  xms = gname("~jma", move.rev(jma))
  assert xms == "~J?"
  
  move.set_name(move.rev(jma), "XY")
  assert move.has_name(jma)
  xms = gname("jma", jma)
  assert xms == "XY"
  xms = gname("~jma", move.rev(jma))
  assert xms == "~XY"

  TRS, JMS = move_example.misc_C(mp_fill, mp_jump)
  
  xms = gname("TRS[2]", TRS[2])
  assert xms == "Tb0"
  
  xms = gname("~JMS[2]", move.rev(JMS[2]))
  assert xms == "~Jd0"
  
  sys.stderr.write("  ... applying {move.tag_names} ...\n")
  move.tag_names([TRS[2], JMS[2]], "Tag.")
  
  xms = gname("TRS[2]", TRS[2])
  assert xms == "Tag.Tb0"
  
  xms = gname("~TRS[2]", move.rev(TRS[2]))
  assert xms == "~Tag.Tb0"
  
  xms = gname("~JMS[2]", move.rev(JMS[2]))
  assert xms == "~Tag.Jd0"

  return
  # ----------------------------------------------------------------------

def test_show():
  sys.stderr.write("--- testing {show,show_list} ---\n");
  TRS, JMS = move_example.misc_C(mp_fill, mp_jump)
 
  sys.stderr.write("  ... {show} ...\n")
  wna = 5
  sys.stderr.write("\n")
  for omv in TRS[0], move.rev(TRS[0]), JMS[0], move.rev(JMS[0]):
    move.show(sys.stderr, "    [", omv, "]\n", wna)
    wna = wna + 5
  sys.stderr.write("\n")
  
  sys.stderr.write("  ... {show_list} ...\n")
  move.show_list(sys.stderr, "    ", TRS, None)
  sys.stderr.write("\n")
  move.show_list(sys.stderr, "    ", JMS, None)
  return
  # ----------------------------------------------------------------------

# Run the tests:

test_shared_border()
test_name()
test_show()
test_make_width_is_jump()
test_move_orient_unpack_rev()
test_displace()
test_plot()
test_fabtime()
test_cover_time()
 
