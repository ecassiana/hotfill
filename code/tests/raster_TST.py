#! /usr/bin/python3
# Test program for module {raster}

import raster
import path
import rootray_shape
import raster_example
import path_example
import block_example
import move 
import move_parms
import input_data
import contact
import hacks
import palette

import rn
import rmxn

import pyx
import sys
from math import sqrt, sin, cos, floor, ceil, inf, nan, pi

# Some dynamics parameters:
ac = 3000
sp_cont = 20
sp_fill = 40
sp_jump = 130
ud_jump = 0.05

sp_link = sp_fill

# Move parameters:
wd_cont = 0.50; mp_cont = move_parms.make(wd_cont, ac, sp_cont, 0.00)
wd_fill = 1.00; mp_fill = move_parms.make(wd_fill, ac, sp_fill, 0.00)
wd_link = 0.70; mp_link = move_parms.make(wd_link, ac, sp_link, 0.00)
wd_jump = 0.00; mp_jump = move_parms.make(wd_jump, ac, sp_jump, ud_jump)

wd_axes = 0.15*wd_fill

def plot(OPHS, LKS, CTS, B, tag):
  # Plots the paths {OPHS}, the links {LKS}, and the contacts {CTS}.
  # Plot files are "tests/out/raster_regroup_TST_{tag}.{ext}"
  # where {ext} is "png", "jpg", "eps".
  
  nph = len(OPHS)
  nlk = len(LKS)
  nct = len(CTS)
  
  assert nph + nlk + nct > 0, "nothing to plot"
  
  # Grow the box as needed to include everything:
  B = rn.box_join(B, path.bbox(OPHS))
  B = rn.box_join(B, path.bbox(LKS))
  B = rn.box_join(B, contact.bbox(CTS))
  
  dp = (0,0)
  
  frame = False
  grid = True
  c,szx,szy = hacks.make_canvas(hacks.round_box(B,0.5), dp, None, frame, grid, 1, 1)
  
  if nph > 0:
    CLRS_ph = hacks.trace_colors(nph, None)
    rwd_ph = 0.80
    axes_ph = False
    dots_ph = False
    arrows_ph = True
    matter_ph = True
    path.plot_standard(c, OPHS, dp, None, CLRS_ph, rwd_ph, wd_axes, axes_ph, dots_ph, arrows_ph, matter_ph)

  if nlk > 0:
    CLRS_lk = hacks.link_colors(nph, None)
    rwd_lk = 0.60
    axes_lk = True
    dots_lk = True
    arrows_lk = False
    matter_lk = False
    path.plot_standard(c, LKS, dp, None, CLRS_lk, rwd_lk, wd_axes, axes_lk, dots_lk, arrows_lk, matter_lk)
    
  if nct > 0:
    clr_ct = pyx.color.rgb(1.000, 0.100, 0.000)
    wd_ct = 2*wd_axes
    sz_tics_ct = 0
    arrows_ct = False
    for ct in CTS:
      contact.plot_single(c, ct, dp, clr_ct, None, 0, wd_ct, sz_tics_ct, arrows_ct)

  fname = "tests/out/raster_TST" "_" + tag
  hacks.write_plot(c, fname)
  return
  # ----------------------------------------------------------------------

def test_basic():
  sys.stderr.write("--- testing {get_spacing_and_phase,path_scanline_index,sort_by_scanline,separate_by_scanline} ---\n")
  
  xdir = ( cos(pi/6), sin(pi/6) )
  ydir = ( -xdir[1], +xdir[0] )
  ystep_exp = wd_fill
  yphase_exp = 0.48*wd_fill
  eps = 0.15*wd_fill
  OPHS = raster_example.rasters_A(mp_fill, xdir, ydir, yphase_exp, eps)
  nph = len(OPHS)
  assert nph == 8
  
  ystep, yphase = raster.get_spacing_and_phase(OPHS, xdir, ydir)
  sys.stderr.write("ystep = %20.16f yphase = %20.16f\n" % (ystep, yphase))
  assert abs(ystep - ystep_exp) <= 0.5*eps
  assert abs(yphase - yphase_exp) <= 0.5*eps
  
  assert raster.path_scanline_index(OPHS[0], xdir, ydir, ystep, yphase) == 0
  assert raster.path_scanline_index(OPHS[1], xdir, ydir, ystep, yphase) == 0
  assert raster.path_scanline_index(OPHS[2], xdir, ydir, ystep, yphase) == 1
  assert raster.path_scanline_index(OPHS[3], xdir, ydir, ystep, yphase) == 2
  assert raster.path_scanline_index(OPHS[4], xdir, ydir, ystep, yphase) == 2
  assert raster.path_scanline_index(OPHS[5], xdir, ydir, ystep, yphase) == 2
  assert raster.path_scanline_index(OPHS[6], xdir, ydir, ystep, yphase) == 3
  assert raster.path_scanline_index(OPHS[7], xdir, ydir, ystep, yphase) == 3
  
  # Scramble the paths:
  OPHS_scr = [ OPHS[3], OPHS[4], OPHS[5], OPHS[2], OPHS[6], OPHS[0], OPHS[1], OPHS[7], ]
  
  # Sort them:
  OPHS_srt = raster.sort_by_scanline(OPHS_scr, xdir, ydir, ystep, yphase)
  assert len(OPHS_srt) == nph
  for k in range(nph):
    assert OPHS_srt[k] == OPHS[k]

  # Separate by scanline:
  SCS = raster.separate_by_scanline(OPHS, xdir, ydir, ystep, yphase)
  sys.stderr.write("len(SCS) = %d\n" % len(SCS))
  sys.stderr.write("SCS = %s\n" % SCS)
  
  assert len(SCS) == 4
  assert len(SCS[0]) == 2
  assert SCS[0][0] == 0
  assert SCS[0][1] == 1
  assert len(SCS[1]) == 1
  assert SCS[1][0] == 2
  assert len(SCS[2]) == 3
  assert SCS[2][0] == 3
  assert SCS[2][1] == 4
  assert SCS[2][2] == 5
  assert len(SCS[3]) == 2
  assert SCS[3][0] == 6
  assert SCS[3][1] == 7
  
  return
  # ----------------------------------------------------------------------
  
def test_links_and_contacts():
  sys.stderr.write("--- testing {make_raster_raster_link,make_raster_raster_contact,create_all_raster_raster_links_and_contacts} ---\n")
  
  debug = False

  xdir = ( cos(pi/6), sin(pi/6) )
  ydir = ( -xdir[1], +xdir[0] )
  
  dmax = 3.0*wd_fill # max length of link path.
  
  yphase = 0.48*wd_fill
  eps = 0.03*wd_fill
  OPHS = raster_example.rasters_A(mp_fill, xdir, ydir, yphase, eps)
  assert len(OPHS) == 8
  
  sys.stderr.write("  ... raster elements ...\n")
  path.show_list(sys.stderr, "    ", OPHS, None, True, True)
  
  mp_link = mp_fill
 
  for oph in OPHS: 
    path.clear_links(oph)
    path.clear_links(path.rev(oph))
    for isd in range(2):
      path.clear_contacts(oph, isd)

  def do_test_single(opha, ophb, dxmid, name, LKS, CTS):
    if debug:
      sys.stderr.write("  ... testing single ...\n")
      path.show(sys.stderr, "    opha = ", opha, "\n", True, 0,0,0)
      path.show(sys.stderr, "    ophb = ", ophb, "\n", True, 0,0,0)
    
    mva,dra = move.unpack(path.elem(opha,0))
    mvb,drb = move.unpack(path.elem(ophb,0))
    lk = raster.make_raster_raster_link(opha, ophb, xdir, dmax, dxmid, mp_link)
    path.set_name(lk, "L" + name, True)
    if debug: path.show(sys.stderr, "    link = ", lk, "\n", True, 0,0,0)
    
    assert lk != None
    assert path.pini(lk) == path.pini(opha)
    assert path.pfin(lk) == path.pini(ophb)
    LKSa = path.get_links(opha)
    LKSb = path.get_links(ophb)
    if debug:
      sys.stderr.write("get_links(opha) = %s\n" % LKSa)
      sys.stderr.write("get_links(ophb) = %s\n" % LKSb)
    assert len(LKSa) == 1
    assert len(LKSb) == 1
    path.compare(LKSa[0], path.rev(LKSb[0]), 1.0e-6, True)
    LKS.append(lk)
    
    szmin = 1.1*wd_fill
    rszmin = 0.01
    tol = 0.20*wd_fill
    ct = raster.make_raster_raster_contact(opha, ophb, xdir, szmin, rszmin, tol)
    contact.set_name(lk, "C" + name)

    contact.show(sys.stderr, "contact = ", ct, "\n", 0)

    assert ct != None
    assert contact.side_move(ct,0) == mva
    assert contact.side_move(ct,1) == mvb
    assert ct in path.get_contacts(opha, 0)
    assert ct in path.get_contacts(ophb, 1)
    CTS.append(ct)
    
    return 
    # ......................................................................
  
  LKS = []
  CTS = []
  
  do_test_single(OPHS[0], OPHS[2], None, "00", LKS, CTS)
  do_test_single(path.rev(OPHS[1]), path.rev(OPHS[2]), None, "20", LKS, CTS)
  do_test_single(path.rev(OPHS[4]), path.rev(OPHS[7]), None, "23", LKS, CTS)

  # Clear all connections:
  for oph in OPHS: 
    for ophr in oph, path.rev(oph):
      path.clear_links(ophr)
      path.clear_contacts(ophr, 0)
      path.clear_contacts(ophr, 1)
  for ct in CTS:
    for isd in range(2):
      contact.clear_side_paths(ct, isd)

  dmax = 3.0*wd_fill
  dxmid = 0.25*wd_fill
  LKS,CTS = raster.create_all_raster_raster_links_and_contacts(OPHS, xdir, ydir, dmax, dxmid, mp_link)
  
  sys.stderr.write("  ... link paths ...\n")
  path.show_list(sys.stderr, "    ", LKS, None, False, False)
  
  sys.stderr.write("  ... contacts ...\n")
  contact.show_list(sys.stderr, "    ", CTS, None)

  # ??? Check {contact.get_side_lists} etc. more ???
  
  plot(OPHS, LKS, CTS, None, "test_contacts")
  
  return
  # ----------------------------------------------------------------------

def make_rect_shape(name, plo, phi, mag, sgn):
  # Makes a closed rectangular shape with those two corners, scaled by {mag}.
  # If {sgn} is {+1} the path is CCW, if {sgn} is {-1} it is CW.
  
  S = rootray_shape.box(rn.scale(mag, plo),  rn.scale(mag, phi))
  if sgn == -1: S = rootray_shape.complement(S)
  return S
  # ----------------------------------------------------------------------
  
def make_circ_shape(name, ctr, R, mag, sgn):
  # Makes a closed circular shape with center {ctr} and radius {R}, scaled by {mag}.
  # If {sgn} is {+1} the path is CCW, if {sgn} is {-1} it is CW.
  
  C = rootray_shape.quadric(2, 2, -1)
  Rm = mag*R
  Adir = ((Rm, 0), (0, Rm));     bdir = rn.scale(mag,ctr)
  Ainv = ((1/Rm, 0), (0, 1/Rm)); binv = rn.scale(-1, rmxn.map_row(bdir, Ainv))
  S = rootray_shape.affine(C, Ainv, binv)
  if sgn == -1: S = rootray_shape.complement(S)
  return S
  # ----------------------------------------------------------------------

def test_from_shape():
  sys.stderr.write("--- testing {endpoints_from_shape} ---\n")

  xdir = ( 1, 0 )
  ydir = ( 0, 1 )
  
  # Some nested contours:
  mag = 3
  A = make_rect_shape("A", ( 0, 0), (16,10), mag, +1) 
  B = make_rect_shape("B", ( 1, 1), ( 2, 9), mag, -1) 
  C = make_rect_shape("C", ( 3, 1), (10, 9), mag, -1) 
  D = make_rect_shape("D", (11, 6), (15, 9), mag, -1) 
  E = make_circ_shape("E", (13, 3), 2, mag, -1)
  F = make_rect_shape("F", ( 4, 6), ( 9, 8), mag, +1) 
  G = make_rect_shape("G", ( 4, 2), ( 5, 5), mag, +1) 
  H = make_rect_shape("H", ( 6, 2), ( 9, 5), mag, +1) 
  I = make_rect_shape("I", ( 7, 3), ( 8, 4), mag, -1) 
  J = make_circ_shape("J", (13, 3), 1, mag, +1)
  U = make_rect_shape("U", (17, 0), (19, 4), mag, +1) 
  V = make_rect_shape("V", (17, 8), (19,10), mag, +1) 

  Box = ((mag*0,mag*0), (mag*19,mag*10))
  mrg = (2,2) # Extra box margin.
  Box = rn.box_expand(Box, mrg, mrg)

  # ----------------------------------------------------------------------
  
  SH = J
  test_from_shape_aux(SH, Box, xdir,ydir, "10_circ")

  # ----------------------------------------------------------------------
  
  SH = rootray_shape.union([ A, U ])
  test_from_shape_aux(SH, Box, xdir,ydir, "20_2sqr")

  # ----------------------------------------------------------------------
  
  SH = rootray_shape.intersection([ A, C ])
  test_from_shape_aux(SH, Box, xdir,ydir, "30_sqhole")

  # ----------------------------------------------------------------------
  
  SH = rootray_shape.intersection([ A, rootray_shape.complement(F) ])
  test_from_shape_aux(SH, Box, xdir,ydir, "40_revhole")
  
  # ----------------------------------------------------------------------
  
  xdir = ( cos(pi/6), sin(pi/6) )
  ydir = ( -xdir[1], +xdir[0] )
    
  HI = rootray_shape.intersection([ H, I ])
  EJ = rootray_shape.union([ E, J ])
  CFGHI = rootray_shape.union([ C, F, G, HI ])
  BCDEFGHIJ = rootray_shape.intersection([ B, CFGHI, D, EJ ])
  ABCDEFGHIJ = rootray_shape.intersection([ A, BCDEFGHIJ ])
  ABCDEFGHIJUV = rootray_shape.union([ ABCDEFGHIJ, U, V ])
  
  # SH = EJ
  # SH = CFGHI 
  # SH = BCDEFGHIJ
  # SH = ABCDEFGHIJ
  SH = ABCDEFGHIJUV
  
  # Rotate the shape so sides are parallel to {xdir,ydir}:
  Adir = (xdir, ydir);   bdir = (0,0)
  Ainv = ((xdir[0], ydir[0]), (xdir[1],ydir[1])); binv = (0,0)
  SH = rootray_shape.affine(SH, Ainv, binv)
  Box = rmxn.box_affine(Box, Adir, bdir)
  
  test_from_shape_aux(SH, Box, xdir,ydir, "50_bigtilt")

  return
  # ----------------------------------------------------------------------
  
def test_from_shape_aux(SH, B, xdir,ydir, tag):
  
  ystep = wd_fill
  yphase = 0.48*wd_fill

  sys.stderr.write("... sub-test %s...\n" % tag)
  
  PTRS = raster.endpoints_from_shape(SH, B, xdir, ydir, ystep, yphase)
  OPHS = raster.from_endpoints(PTRS, mp_fill)
  
  plot(OPHS, [], [], B, "test_endpoints_from_shape_" + tag)
  
  return
  # ----------------------------------------------------------------------

def make_rect_poly(plo, phi, mag, sgn):
  # Makes a closed rectangular polygon with the given corners, scaled by {mag}.
  # If {sgn} is {+1} the path is CCW, if {sgn} is {-1} it is CW.
  
  PTS = [ 
    plo,
    (phi[0], plo[1]),
    phi,
    (plo[0], phi[1]),
    plo,
  ]
  ori = hacks.poly_orientation(PTS)
  if sgn != ori: PTS = list(reversed(PTS))
  # Apply {mag} scale:
  PTS = [ rn.scale(mag, p) for p in PTS]
  return PTS
  # ----------------------------------------------------------------------
  
def make_circ_poly(ctr, R, mag, sgn):
  # Makes a closed circular shape with center {ctr} and radius {R}, scaled by {mag}.
  # If {sgn} is {+1} the path is CCW, if {sgn} is {-1} it is CW.
  
  PTS = []
  dstep = 0.25
  nsteps = max(8, int(ceil(2*pi*R/dstep)))
  for ipt in range(nsteps):
    a = sgn*2*pi*(ipt/nsteps)
    p = rn.add(ctr, (R*cos(a), R*sin(a)))
    PTS.append(p)
  # Apply {mag} scale:
  PTS = [ rn.scale(mag, p) for p in PTS]

  # Close the polygon:
  PTS.append(PTS[0])
  return PTS
  # ----------------------------------------------------------------------

def test_from_polys():
  sys.stderr.write("--- testing {endpoints_from_polys} ---\n")

  # Some nested contours:
  mag = 3
  PTCS_1 = [
    make_rect_poly(( 0, 0), (16,10), mag, +1), # "A"
    make_rect_poly(( 1, 1), ( 2, 9), mag, -1), # "B"
  ]                                            
                                               
  PTCS_2 = [                                   
    make_rect_poly(( 3, 1), (10, 9), mag, -1), # "C"
    make_rect_poly((11, 6), (15, 9), mag, -1), # "D"
    make_circ_poly((13, 3), 2, mag, -1),       # "E"
    make_rect_poly(( 4, 6), ( 9, 8), mag, +1), # "F"
    make_rect_poly(( 4, 2), ( 5, 5), mag, +1), # "G"
    make_rect_poly(( 6, 2), ( 9, 5), mag, +1), # "H"
    make_rect_poly(( 7, 3), ( 8, 4), mag, -1), # "I"
    make_circ_poly((13, 3), 1, mag, +1),       # "J"
    make_rect_poly((17, 0), (19, 4), mag, +1), # "U"
    make_rect_poly((17, 8), (19,10), mag, +1), # "V"
  ]

  Box = ((mag*0,mag*0), (mag*19,mag*10))
  mrg = (2,2) # Extra box margin.
  Box = rn.box_expand(Box, mrg,mrg)
  
  xdir = ( cos(pi/6), sin(pi/6) )
  ydir = ( -xdir[1], +xdir[0] )
  
  PTCS = PTCS_1; tag = "1"
  test_from_polys_aux(PTCS, Box, xdir,ydir, tag)
  
  PTCS = PTCS_1 + PTCS_2; tag = "2"
  test_from_polys_aux(PTCS, Box, xdir,ydir, tag)
  
  return
  # ----------------------------------------------------------------------
  
def test_from_polys_aux(PTCS, B, xdir,ydir, tag):

  sys.stderr.write("... sub-test %s...\n" % tag)
  
  ystep = wd_fill
  yphase = 0.48*wd_fill

  PTRS = raster.endpoints_from_polys(PTCS, xdir, ydir, ystep, yphase)
  OPHS = raster.from_endpoints(PTRS, mp_fill)
  
  plot(OPHS, [], [], B, "test_endpoints_from_polys_" + tag)

  return
  # ----------------------------------------------------------------------

def test_analyze_fabtime():

  sys.stderr.write("--- testing {analyze_fabtime} ---\n")
  
  tag = "test_analyze_fabtime"

  ang = pi/6
  xdir = (cos(ang), sin(ang)) 
  ydir = (-xdir[1], xdir[0]) 
  
  # Create two paths:
  nmv0 = 9
  nmv1 = 2
  pt = [[None]*(nmv0+1), [None]*(nmv1+1)]
  mpp = [[None]*nmv0, [None]*nmv1]
  
  pt[0][0] = rn.mix(1.000,xdir, 1.000,ydir); mpp[0][0] = mp_jump
  pt[0][1] = rn.mix(2.000,xdir, 1.000,ydir); mpp[0][1] = mp_fill
  pt[0][2] = rn.mix(3.500,xdir, 1.000,ydir); mpp[0][2] = mp_jump
  pt[0][3] = rn.mix(4.000,xdir, 1.000,ydir); mpp[0][3] = mp_jump
  pt[0][4] = rn.mix(5.000,xdir, 1.000,ydir); mpp[0][4] = mp_link
  pt[0][5] = rn.mix(5.001,xdir, 1.001,ydir); mpp[0][5] = mp_link
  pt[0][6] = rn.mix(6.000,xdir, 2.000,ydir); mpp[0][6] = mp_fill
  pt[0][7] = rn.mix(3.000,xdir, 2.000,ydir); mpp[0][7] = mp_fill
  pt[0][8] = rn.mix(2.200,xdir, 2.000,ydir); mpp[0][8] = mp_jump
  pt[0][9] = rn.mix(1.000,xdir, 2.000,ydir)
  
  pt[1][0] = rn.mix(1.000,xdir, 3.000,ydir); mpp[1][0] = mp_fill
  pt[1][1] = rn.mix(1.001,xdir, 3.000,ydir); mpp[1][1] = mp_fill
  pt[1][2] = rn.mix(4.000,xdir, 3.000,ydir)
  
  OPHS = []
  Trast_exp = 0
  Tlink_exp = 0
  Tjump_exp = 0
  for iph in range(len(pt)):
    nmv = len(pt[iph]) - 1
    OMVS = []
    for imv in range(nmv):
      p = pt[iph][imv]
      q = pt[iph][imv+1]
      mp = mpp[iph][imv]
      omv = move.make(p, q, mp)
      Tfab_mv = move.fabtime(omv)
      if mp == mp_jump:
        udp = move_parms.ud_penalty(mp_jump)
        if imv > 0 and not move_parms.is_jump(mpp[iph][imv-1]):
          Tjump_exp += udp
        Tjump_exp += Tfab_mv
        if imv < nmv-1 and not move_parms.is_jump(mpp[iph][imv+1]):
          Tjump_exp += udp
      else:
        if mp == mp_fill:
          Trast_exp += Tfab_mv
        elif mp == mp_link:
          Tlink_exp += Tfab_mv
        else:
          assert False
      OMVS.append(omv)
    oph = path.from_moves(OMVS)
    OPHS.append(oph)

  ytol = 1.0e-7;
  minlen = 0.1
  Trast_cmp, Tlink_cmp, Tjump_cmp = raster.analyze_fabtime(OPHS, ydir, ytol, minlen)

  assert abs(Trast_cmp - Trast_exp) < 1.0e-8
  assert abs(Tlink_cmp - Tlink_exp) < 1.0e-8
  assert abs(Tjump_cmp - Tjump_exp) < 1.0e-8
  
  B = path.bbox(OPHS)
  mrg = (0.5,0.5) # Extra box margin.
  B = rn.box_expand(B, mrg,mrg)
  
  plot(OPHS, [], [], B, tag)
  
  return
  # ----------------------------------------------------------------------

test_analyze_fabtime()
### test_links_and_contacts()
### test_from_polys()
### test_basic()
### test_from_shape()
