#! /usr/bin/python3
# Test program for module {hacks}
# Last edited on 2021-10-29 03:07:30 by stolfi

import hacks
import color

import rn

import pyx
import sys
from math import sqrt, log, exp, sin, cos, floor, ceil, inf, nan, pi

def test_adjust_dash_pattern():
  sys.stderr.write("--- testing {adjust_dash_pattern} ---\n")
  rlen = 0.3 # Dash length relative to {wd}.
  rgap = 0.2 # Gap length relative to {wd}.
  sys.stderr.write("ideal ratios: rlen = %.6f rgap = %.6f\n" % (rlen, rgap))
  D1 = tuple( 0.1*k for k in range(11) )
  D2 = tuple( 10+0.1*k for k in range(31) )
  for rdist in D1 + D2:
    sys.stderr.write("rdist = %.6f" % rdist)
    rdsh = hacks.adjust_dash_pattern(rdist, rlen, rgap)
    sys.stderr.write(" rdsh adjusted = %s" % str(rdsh))
    if rdsh != None:
      assert type(rdsh) is list or type(rdsh) is tuple
      assert len(rdsh) == 2
      rlen1, rgap1 = rdsh
      n = int(floor((rdist-rlen1)/(rlen1+rgap1) + 0.5)) # Number of gaps.
      rdist1 = n*(rlen1+rgap1) + rlen1
      sys.stderr.write(" n = %d rdist1 = %.6f\n" % (n, rdist1))
      assert abs(rdist - rdist1) < 1.0e-8*rdist
    else:
      sys.stderr.write(" solid\n")
  return
  # ----------------------------------------------------------------------


def test_make_canvas():
  sys.stderr.write("--- testing {make_canvas} ---\n")
  
  dp = (2, 2)
  B = ((2, 3), (5, 7))
  nx = 2
  ny = 3
  for frame in False, True:
    for grid in False, True:
      tag = "make_canvas_fr%s_gr%s" % ("FT"[frame], "FT"[grid])
      c,szx,szy = hacks.make_canvas(hacks.round_box(B,0.5), dp, None, frame, grid, 2, 3)
      hacks.write_plot(c, "tests/out/hacks_TST_plot_" + tag)
  return 
  # ----------------------------------------------------------------------
  
def test_basic_plot():
  sys.stderr.write("--- testing {pyx,round_box,plot_box,plot_frame,plot_grid} ---\n")
  tag = "basic_plot"

  # Bounding box of one drawing:
  B = ((0.5, 0.5), (10.5, 6.5))
  Bplot = hacks.round_box(B, 2)

  dp = (3,1)

  scale = 1.0
  sys.stderr.write("scale = %8.4f\n" % scale)
  c,szx,szy = hacks.make_canvas(Bplot, dp, None, False, True, 1, 1)

  tan = pyx.color.rgb(0.800, 0.700, 0.600)
  glo = pyx.color.rgb(0.400, 0.900, 1.000)

  sys.stderr.write("  ... plot_box ...\n")
  BA = ((3.1,1.8), (3.9,2.2))

  clr = pyx.color.rgb.green
  hacks.plot_box(c, BA, dp, clr)

  sys.stderr.write("  ... plot_frame ...\n")
  wd_frame_a = 0.06
  wd_frame_b = 0.02
  dash_frame = ( 3.0*wd_frame_a, 2.0*wd_frame_a )
  hacks.plot_frame(c, B, dp, pyx.color.rgb.red,   wd_frame_a, dash_frame, +0.20)
  hacks.plot_frame(c, B, dp, pyx.color.rgb.black, wd_frame_b, None,       00.00)
  hacks.plot_frame(c, B, dp, pyx.color.rgb.blue,  wd_frame_a, dash_frame, -0.20)

  sys.stderr.write("  ... plot_grid ...\n")
  dp = (dp[0], dp[1] + szy)
  wd_grid_a = 0.01
  wd_grid_b = 0.03
  dash_grid = ( 2*wd_grid_a, 1*wd_grid_a )
  hacks.plot_grid(c, B, dp, tan, wd_grid_a, dash_grid, -0.50, 0.50,0.20)
  hacks.plot_grid(c, B, dp, glo, wd_grid_b, None,      -0.35, 1.00,1.00)

  hacks.write_plot(c, "tests/out/hacks_TST_plot_" + tag)
  return 
  # ----------------------------------------------------------------------
   
def do_test_plot_line(c, X, Y, Xsz, dp, wd, dsh):
  clr = pyx.color.rgb.blue
  p = ( X, Y )
  q = rn.add(p, ( Xsz, 0 ))
  hacks.plot_line(c, p, q, dp, clr, wd, dsh)
  
  if dsh != None:
    # Simulate the dashes by hand:
    dsh = hacks.adjust_dash_pattern(rn.dist(p,q), dsh[0], dsh[1])
    if dsh != None:
      u = p
      for k in range(4):
        v = rn.add(u, (dsh[0], 0))
        hacks.plot_line(c, u, v, dp, None, 0.05, None)
        u = rn.add(v, (dsh[1], 0))

  fsize = 12

  r = rn.add(p, ( 0, -0.5*wd - 0.036*fsize))
  tx0 = (r"wd = %.2f xsz = %.2f" % (wd, Xsz))
  if wd > 0:
    tx0 += (r" rsz = %.3f" % (Xsz/wd))
  tx0 = r"\texttt{" + tx0 + "}"
  hacks.plot_text(c, tx0, r, dp, fsize, None)

  if dsh != None:
    s = rn.add(r, ( 0, -0.036*fsize))
    tx1 = (r"dsh = (%.2f,%.2f)" % (dsh[0],dsh[1]))
    if wd > 0:
      rsh = rn.scale(1/wd, dsh)
      tx1 += ( r" rsh = (%.3f,%.3f)" % (rsh[0],rsh[1]))
    tx1 = r"\texttt{" + tx1 + "}"
    hacks.plot_text(c, tx1, s, dp, fsize, None)
  return
  # ----------------------------------------------------------------------
 
def test_plot_line_solid():
  sys.stderr.write("--- testing {plot_line} (solid) ---\n")
  tag = "plot_line_solid"

  # Bounding box of one drawing:
  B = ((0.5, 0.5), (8.5, 12.5))
  Bplot = hacks.round_box(B, 1)

  dp = (1,1)

  c,szx,szy = hacks.make_canvas(Bplot, dp, None, False, True, 1, 1)

  nwd = 10 # Number of line widths including zero
  wdmax = 0.50 # Max line width (mm).
  wdmin = 0.03 # Min nonzero line width (mm).
  Xlo = B[0][0] + 0.5
  Ylo = B[0][1] + 0.5
  Ytot = B[1][1] - B[0][1] - 1.0 # Total space for the lines.
  Ystep = Ytot/(nwd-1) # Y disp between different widths.

  ll = 10 # Length of lines relative to {wd}.
  for kwd in range(nwd):
    wd = 0 if kwd == 0 else wdmin*exp(log(wdmax/wdmin)*(kwd-1)/(nwd-2))
    clr = pyx.color.rgb.blue
    X = Xlo
    Y = Ylo + kwd*Ystep
    Xsz = ll*wd
    do_test_plot_line(c, X, Y, Xsz, dp, wd, None) 

  hacks.write_plot(c, "tests/out/hacks_TST_plot_" + tag)

  return 
  # ----------------------------------------------------------------------

def test_plot_line_dashed(wd):
  sys.stderr.write("--- testing {plot_line} (dashed) ---\n")
  tag = "plot_line_dashed_%03d" % int(floor(100*wd + 0.5))

  # Parameters for the left-hand plot: 
  llmax = 15                # Max length of lines relative to {wd} for left half.
  llmin = 0.75*llmax        # Min length of lines relative to {wd} for left half.

  # Parameters for the right-hand plot:
  ndash = 4                 # Target number of dashes for right half.
  rsh_len_max = 3.0         # Max dash length rel to {wd}.
  rsh_len_min = 1.0         # Min dash length rel to {wd}.
  rsh_gap_std = 2.0         # Standard gap width rel to {wd}

  # Bounding box of one drawing:
  Xszmax_wd = wd*llmax # Est width of left half of plot.
  Xszmax_sh = wd*(ndash*rsh_len_max + (ndash-1)*rsh_gap_std)
  
  B = ((0,0), ( Xszmax_wd + Xszmax_sh + 8, 15))
  Bplot = hacks.round_box(B, 1.0)

  dp = (3,1)
  tan = pyx.color.rgb(0.800, 0.700, 0.600)

  scale = 1.0
  sys.stderr.write("scale = %8.4f\n" % scale)
  c,szx,szy = hacks.make_canvas(Bplot, dp, None, False, False, 1, 1)
  hacks.plot_grid(c, B, dp, tan, 0.01, None, -0.5, 0.20, 0.20)
  hacks.plot_grid(c, B, dp, tan, 0.02, None, -0.5, 1.00, 1.00)

  Xlo = ceil(B[0][0] + 0.5)
  Ylo = ceil(B[0][1] + 0.5)
  Xtot = floor(B[1][0] - 0.5) - Xlo # Total X space for the lines.
  Ytot = floor(B[1][1] - 0.5) - Ylo # Total Y space for the lines.
  
  # Varying the line length:
  nll = 5 # Number of line lengths including non-dashed line.
  Xlo_ll = Xlo
  Ystep_ll = Ytot/(nll-1) # Y disp between different line lengths.
  for kll in range(nll):
    X = Xlo_ll
    Y = Ylo + kll*Ystep_ll
    if kll == 0:
      # Fixed length:
      ll = llmax
      rsh = None
    else:
      # Vary length of line, keep ideal dash pattern:
      ll = llmax if kll == 0 else llmax*exp(log(llmin/llmax)*(kll-1)/(nll-2))
      rsh = (2.0, rsh_gap_std)
    Xsz = ll*wd
    sys.stderr.write("\nout: rdashpat = %s\n" % str(rsh))
    dsh = rn.scale(wd, rsh) if kll > 0 else None
    do_test_plot_line(c, X, Y, Xsz, dp, wd, dsh) 

  # Varying the dash pattern:
  nsh = 7 # Count of dash paterns including solid line.
  Ystep_sh = Ytot/(nsh-1) # Y disp between different line lengths.
  Xlo_sh = ceil(Xlo + Xtot/2) + 1
  lldef = 15 # Length of solid line, rel to {wd}
  for ksh in range(nsh):
    X = Xlo_sh
    Y = Ylo + ksh*Ystep_sh
    if ksh == 0:
      # Solid line of max len:
      ll = llmax
      rsh = None
    else:
      # Vary dash length, try adjust length of line to exact dashes:
      rlen = rsh_len_min + (rsh_len_max - rsh_len_min)*(ksh-1)/(nsh-2)
      rgap = rsh_gap_std
      rsh = (rlen, rgap) 
      ll = (ndash*rlen + (ndash-1)*rgap) + 0.001
    Xsz = ll*wd
    sys.stderr.write("\nout: rdashpat = %s\n" % str(rsh))
    dsh = rn.scale(wd, rsh) if ksh > 0 else None
    do_test_plot_line(c, X, Y, Xsz, dp, wd, dsh) 


  hacks.write_plot(c, "tests/out/hacks_TST_plot_" + tag)

  return 
  # ----------------------------------------------------------------------
  
def test_plot_text():
  sys.stderr.write("--- testing {plot_text} ---\n")
  tag = "plot_text"

  # Bounding box of one drawing:
  B = ((0.5, 0.5), (10.5, 6.5))
  Bplot = hacks.round_box(B, 1)

  dp = (3,1)

  scale = 1.0
  sys.stderr.write("scale = %8.4f\n" % scale)
  c,szx,szy = hacks.make_canvas(Bplot, dp, None, False, True, 1, 1)

  Xlo = B[0][0] + 0.5
  Ylo = B[0][1] + 0.5
  Ytot = B[1][1] - B[0][1] - 1.0 # Total space for the lines.

  nfs = 4 # Number of font sizes
  fsmin = 5
  fsmax = 20
  Ystep = Ytot/(nfs-1)
  for kfs in range(nfs):
    X = Xlo
    Y = Ylo + kfs*Ystep
    fsize = floor(fsmin*exp(log(fsmax/fsmin)*kfs/(nfs-1)) + 0.5)
    p = (X, Y)
    tx = (r"\textsf{A}")
    hacks.plot_text(c, tx, p, dp, fsize, 'red')

    p = rn.add(p, (0.03*fsize, 0))
    tx = (r"$\mathcal{B}$")
    hacks.plot_text(c, tx, p, dp, fsize, None)

    p = rn.add(p, (0.04*fsize, 0))
    tx = (r"$M^{\su{SUP}}_{\su{SUB}}$")
    hacks.plot_text(c, tx, p, dp, fsize, None)

    p = rn.add(p, (0.10*fsize, 0))
    tx = (r"\texttt{%.1f}" % fsize)
    hacks.plot_text(c, tx, p, dp, fsize, None)

  hacks.write_plot(c, "tests/out/hacks_TST_plot_" + tag)
  
  return 
  # ----------------------------------------------------------------------

def test_filet_center():
  p = (30, 50)
  q = (70, 40)
  ro = 20
  rt = 30
  for iord in range(2):
    if iord == 0:
      t = hacks.filet_center(p, q, ro, rt)
    else:
      t = hacks.filet_center(q, p, ro, rt)
    sys.stderr.write("t = ( %10.6f, %10.6f )\n" % (t[0],t[1]))
    assert abs(rn.dist(t,p) - (rt + ro)) < 1.0e-8
    assert abs(rn.dist(t,q) - (rt + ro)) < 1.0e-8
  return
  # ----------------------------------------------------------------------

def test_circle_x():
  ny = 10
  rad = 20
  yctr = 50
  for iy in range(ny+1):
    y = yctr = 1.2*rad*(2*iy/ny - 1)
    x = hacks.circle_x(y, yctr, rad)
    sys.stderr.write("y = %10.6f\n" % y)
    if abs(y - yctr) > rad:
      assert x == 0
    else:
      assert abs(rad - sqrt(x*x + (y-yctr)*(y-yctr))) < 1.0e-8
  return
  # ----------------------------------------------------------------------

def test_poly_orientation():
  sys.stderr.write("--- testing {poly_orientation} ---\n")

  C = (( 3, 1), ( 3, 9), (10, 9), (10, 1), ( 3, 1))
  assert hacks.poly_orientation(C) == -1
  
  F = (( 4, 6), ( 9, 6), ( 9, 8), ( 4, 8), ( 4, 6))
  assert hacks.poly_orientation(F) == +1
  
  K = (
    ( 0, 0), ( 9, 0), ( 9, 9), ( 0, 9), 
    ( 0, 2), ( 7, 2), ( 7, 7), ( 2, 7),
    ( 2, 4), ( 5, 4), ( 5, 5), ( 3, 5),
    ( 3, 6), ( 6, 6), ( 6, 3), ( 2, 3),
    ( 2, 8), ( 8, 8), ( 8, 1), ( 0, 1),
    ( 0, 0),
  )
  assert hacks.poly_orientation(K) == +1
  return
  # ----------------------------------------------------------------------

def test_quatratic_roots():
  sys.stderr.write("--- testing {real_quadratic_roots} ---\n")

  n = 5
  for iA in range(2*n+1):
    for iB in range(2*n+1):
      for iC in range(2*n+1):
        A = (iA-n)/n
        B = (iB-n)/n
        C = (iC-n)/n
        if (A != 0):
          t0, t1 = hacks.real_quadratic_roots(A, B, C)
          for t in t0, t1:
            if t != None:
              v = A*t*t + B*t + C
              if abs(v) > 1.0e-8:
                sys.stderr.write("** error? A = %10.7f B = %10.7f C = %10.7f" % (A, B, C))
                sys.stderr.write("  t = %21.15e v = %21.15e\n" % (t,v))
  return
  # ----------------------------------------------------------------------
  
def test_choose_grid_steps():
  sys.stderr.write("--- testing {choose_grid_steps} ---\n")

  sz = 3.14
  while (sz < 314):
    gstep_big, gstep_sma = hacks.choose_grid_steps(sz)
    sys.stderr.write("  %8.3f %8.3f %8.3f\n" % (sz, gstep_big, gstep_sma))
    sz = 1.2345*sz
    
  return
  # ----------------------------------------------------------------------
  
def test_colors(nc):
  sys.stderr.write("--- testing {trace_colors,link_colors,matter_color} nc = %d ---\n" % nc)

  tag = "colors"

  nY = 4
  nkind = 3 # Traces, links, matter.

  szpatch = (2.0, 3.0) # Size of color patches
  Bpatch = ((0,0), szpatch)
  Bplot = rn.box_expand(Bpatch, (0.25,0.25), (0.25,0.25))

  frame = False
  grid = False
  dp = (0,0)
  c,szx,szy = hacks.make_canvas(Bplot, dp, None, frame, grid, nc,nkind*nY)

  for ikind in range(nkind):
    for iY in range(nY):
      Y = ( None, 0.300, 0.600, 0.900 )[iY]
      # Get Color list {CLRS} as {pyx.color,rgb}, and a two-letter tag {xn}:
      if ikind == 0:
        CLRS = hacks.trace_colors(nc, Y); xn = "tr"
      elif ikind == 1:
        CLRS = hacks.link_colors(nc, Y); xn = "lk"
      elif ikind == 2:
        CLRS = [ hacks.matter_color(Y), ]; xn = "mt"
      else:
        assert False
      # Convert color list to {RGB} float triples:
      RGB = [ ( clr.r, clr.g, clr.b ) for clr in CLRS ]
      # Now plot the list and compute the min pair distance:
      sys.stderr.write("  ... list = %s Y = %s ...\n" % (xn, str(Y)))
      nck = len(CLRS) # Num of colors on this list.
      for k1 in range(nck):
        clrk1 = CLRS[k1]
        RGBk1 = RGB[k1] 
        sys.stderr.write("    RGB = ( %5.3f %5.3f %5.3f )" % (RGBk1[0],RGBk1[1],RGBk1[2]))
        sys.stderr.write(" Y = %5.3f" % color.Y_from_RGB(RGBk1))
        if nck >= 2:
          # Find the nearest neighbor in the table to {RGBk1}:
          RGBmin = None # The nearest neighbor.
          dmin = 2.000  # Its distance from {RGBk1}
          for k2 in range(nc):
            if k1 != k2:
              RGBk2 = RGB[k2]
              dk2 = rn.dist(RGBk1, RGBk2)
              if dk2 < dmin:
                dmin = dk2
                RGBmin = RGBk2
          sys.stderr.write(" dmin = %6.4f  RGBmin = ( %5.3f %5.3f %5.3f )" % (dmin,RGBmin[0],RGBmin[1],RGBmin[2]))
        # Plot a color patch:
        org = (k1*szx, (3*iY + ikind)*szy)
        hacks.plot_box(c, Bpatch, org, clrk1)
        sys.stderr.write("\n")

  hacks.write_plot(c, "tests/out/hacks_TST_plot_" + tag)

  return
  # ----------------------------------------------------------------------

def test_clip_seg_to_strip():
  sys.stderr.write("--- testing {clip_seg_to_strip} ---\n")
  
  xdir,xlen = rn.dir((1,2))
  ydir = (-xdir[1], xdir[0])
  
  i0 = 1
  i1 = 4
  
  def my_clip_seg(ip, iq, i0, i1):
    # Returns rel positions {rp,rq} given indices {ip,iq} of points and {i0,i1} of lines.
    if max(ip,iq) <= i0 or min(ip,iq) >= i1:
      # Should be outside:
      rp, ploc, rq, qloc = None, None, None, None
    elif min(ip,iq) >= 2 and max(ip,iq) <= 3:
      # Should be fully inside:
      rp, ploc, rq, qloc = 0, None, 1, None
    else:
      # Compute new endpoint near {p}:
      if ip <= i0:
        # Enters through line {i0}
        assert iq > i0
        rp, ploc = (i0 - ip)/(iq - ip), 0
      elif ip >= i1:
        # Enters through line {i1}
        assert iq < i1
        rp, ploc = (i1 - ip)/(iq - ip), 1
      else:
        # Starts inside:
        rp, ploc = 0, None
      # Compute new endpoint near {q}:
      if iq <= i0:
        # Exits through line {i0}
        assert ip > i0
        rq, qloc = (i0 - ip)/(iq - ip), 0
      elif iq >= i1:
        # Exits through line {i1}
        assert ip < i1
        rq, qloc = (i1 - ip)/(iq - ip), 1
      else:
        # Ensh inside:
        rq, qloc = 1, None
    return rp, ploc, rq, qloc
    # ....................................................................
  
  def do_test(ip, iq):

    p = rn.mix(12, xdir, ip, ydir)
    q = rn.mix(15, xdir, iq, ydir)

    sys.stderr.write("  ip = %d (%12.8f,%12.8f)" % (ip,p[0],p[1]))
    sys.stderr.write("  iq = %d (%12.8f,%12.8f)" % (iq,q[0],q[1]))
    sys.stderr.write("\n")

    # Compute by module:
    ps_cmp, ploc_cmp, qs_cmp, qloc_cmp = hacks.clip_seg_to_strip(p,q, ydir, i0, i1)

    if ps_cmp != None:
      sys.stderr.write("  ps_cmp = (%12.8f,%12.8f) %s" % (ps_cmp[0],ps_cmp[1],str(ploc_cmp)))
      sys.stderr.write("  qs_cmp = (%12.8f,%12.8f) %s" % (qs_cmp[0],qs_cmp[1],str(qloc_cmp)))
      sys.stderr.write("\n")

    # Compute by indices:
    rp_exp, ploc_exp, rq_exp, qloc_exp = my_clip_seg(ip, iq, i0, i1)
    ps_exp = None if rp_exp == None else rn.mix(1-rp_exp, p, rp_exp, q)
    qs_exp = None if rq_exp == None else rn.mix(1-rq_exp, p, rq_exp, q)

    if ps_exp != None:
      sys.stderr.write("  ps_exp = (%12.8f,%12.8f) %s" % (ps_exp[0],ps_exp[1],str(ploc_exp)))
      sys.stderr.write("  qs_exp = (%12.8f,%12.8f) %s" % (qs_exp[0],qs_exp[1],str(qloc_exp)))
      sys.stderr.write("\n")

    # Check it:
    if ps_cmp == None:
      assert ps_cmp == None and qs_cmp == None
      assert ploc_cmp == None and qloc_cmp == None
      assert ps_exp == None and qs_exp == None
      assert ploc_exp == None and qloc_exp == None
    else:
      assert ps_cmp != None and qs_cmp != None
      assert ps_exp != None and qs_exp != None
      assert rn.dist(ps_cmp, ps_exp) < 1.0e-8
      assert ploc_cmp == ploc_exp

      assert rn.dist(qs_cmp, qs_exp) < 1.0e-8
      assert qloc_cmp == qloc_exp
    return
    # ....................................................................

  for ip in range(6):
    for iq in range(6):
      # Avoid segs that just touch the border:
      if ip != i0 and ip != i1 and iq != i0 and iq != i1:
        do_test(ip, iq)

  return
  # ----------------------------------------------------------------------

def test_trim_poly():
  sys.stderr.write("--- testing {find_point_on_poly,trim_poly} ---\n")
  
  CS = [ (1,1), (5,1), (5,5), (3,5), (2,3), (1,5), (1,3) ]
  
  ia, qa = hacks.find_point_on_poly(CS, 0, 6)
  sys.stderr.write("{find_point_on_poly} i = %d q = ( %10.8f %10.8f )\n" % (ia, qa[0], qa[1]))
  assert ia == 2
  assert rn.dist(qa, (5,3)) < 1.0e-8
  
  ib, qb = hacks.find_point_on_poly(CS, 1, 0.5)
  sys.stderr.write("{find_point_on_poly} i = %d q = ( %10.8f %10.8f )\n" % (ib, qb[0], qb[1]))
  assert ib == 5
  assert rn.dist(qb, (1,3.5)) < 1.0e-8
  
  TS_cmp = hacks.trim_poly(CS, 9, 0.5)
  TS_exp = [ (4,5), (3,5), (2,3), (1,5), (1,3.5) ]
  assert len(TS_cmp) == len(TS_exp)
  for i in range(len(TS_cmp)):
    sys.stderr.write("{trim_poly} i = %d" % i)
    sys.stderr.write(" cmp = ( %10.8f %10.8f )" % (TS_cmp[i][0], TS_cmp[i][1]))
    sys.stderr.write(" exp = ( %10.8f %10.8f )\n" % (TS_exp[i][0], TS_exp[i][1]))
    assert rn.dist(TS_cmp[i], TS_exp[i]) < 1.0e-8

  return
  # ----------------------------------------------------------------------

sys.stderr.write("??? Needs more tests ???\n")

test_trim_poly()

test_clip_seg_to_strip()
test_circle_x()
test_filet_center()

test_quatratic_roots()
test_poly_orientation()

test_adjust_dash_pattern()

test_choose_grid_steps()

test_make_canvas()
test_basic_plot()

test_plot_line_solid()
test_plot_line_dashed(1.00)
test_plot_line_dashed(0.50)

test_plot_text()

test_colors(10)
test_colors(50)




