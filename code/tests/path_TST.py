#! /usr/bin/python3
# Test program for module {path}

import path
import path_example
import raster_example
import move 
import move_parms
import contact 
import block
import hacks
import palette
import job_parms
import rn
import pyx
import random
import sys
from math import sqrt, hypot, sin, cos, floor, ceil, inf, nan, pi

parms = job_parms.typical_js()
mp_jump = move_parms.make_for_jumps(parms)

# Some arbitrary dynamics parameters:
ac = parms['acceleration'] 
sp = parms['max_extrusion_speed'] 

# Move parameters matching the raster spacings in the file:
wd_cont = 0.75; mp_cont = move_parms.make(wd_cont, ac, sp, 0.0)
wd_fill = 1.00; mp_fill = move_parms.make(wd_fill, ac, sp, 0.0)
wd_link = 0.50; mp_link = move_parms.make(wd_link, ac, sp, 0.0)

parms['solid_raster_width'] = wd_fill
parms['contour_trace_width'] = wd_cont

def test_make_trivial():
  sys.stderr.write("--- testing {make_trivial} ---\n")

  p00 = (1,1)

  ph0 = path.make_trivial(p00)

  assert path.nelems(ph0) == 0
  assert path.nelems(path.rev(ph0)) == 0
  assert path.pini(ph0) == p00
  assert path.pfin(ph0) == p00
  
  mvx = move.make(p00, p00, mp_fill)
  assert path.find_move(ph0, mvx) == None

  path.validate(ph0)
  path.validate(path.rev(ph0))

  sys.stderr.write("\n")
  return
  # ----------------------------------------------------------------------

def test_from_move():
  sys.stderr.write("--- testing {from_move,endpoints} ---\n")

  p10 = (1,2)
  p11 = (2,3)

  ph1 = path.from_move(p10, p11, mp_fill)
  assert path.nelems(ph1) == 1
  mv10 = path.elem(ph1, 0)
  assert not move.is_jump(mv10)
  assert path.pini(ph1) == p10
  assert path.pfin(ph1) == p11
  assert path.endpoints(ph1) == (p10, p11)

  ph1r = path.rev(ph1)
  assert path.elem(ph1r, 0) == move.rev(mv10)
  assert path.pini(ph1r) == p11
  assert path.pfin(ph1r) == p10
  assert path.endpoints(ph1r) == (p11, p10)
  
  assert path.find_move(ph1, mv10) == 0
  assert path.find_move(ph1, move.rev(mv10)) == 0
  assert path.find_move(ph1r, mv10) == 0

  mvx = move.make((0,0), (1,1), mp_fill)
  assert path.find_move(ph1, mvx) == None

  path.validate(ph1)
  path.validate(ph1r)

  sys.stderr.write("\n")
  return
  # ----------------------------------------------------------------------

def test_from_moves(): 
  sys.stderr.write("--- testing {from_moves} ---\n")

  p0 = (2,2)
  p1 = (3,3)
  p2 = (2,4)
  p3 = (1,3)
  
  mv01 = move.make(p0, p1, mp_fill)
  mv12 = move.make(p1, p2, mp_jump)
  mv32 = move.make(p3, p2, mp_fill)

  ph = path.from_moves((mv01, mv12, move.rev(mv32)))

  assert path.nelems(ph) == 3
  omv0 = path.elem(ph, 0)
  omv1 = path.elem(ph, 1)
  omv = path.elem(ph, 2)
  assert (move.unpack(omv0)) == (mv01, 0)
  assert (move.unpack(omv1)) == (mv12, 0)
  assert (move.unpack(omv)) == (mv32, 1)
  assert not move.is_jump(omv0)
  assert move.is_jump(omv1)
  assert not move.is_jump(omv)
  assert path.pini(ph) == p0
  assert path.pfin(ph) == p3
  assert path.elem(path.rev(ph), 0) == move.rev(omv)
  assert path.elem(path.rev(ph), 1) == move.rev(omv1)
  assert path.elem(path.rev(ph), 2) == move.rev(omv0)
  assert path.pini(path.rev(ph)) == p3
  assert path.pfin(path.rev(ph)) == p0

  mvx = move.make((0,0), (1,1), mp_fill)
  assert path.find_move(ph, mvx) == None

  path.validate(ph)
  path.validate(path.rev(ph))

  sys.stderr.write("\n")
  return
  # ----------------------------------------------------------------------

def test_from_points():
  sys.stderr.write("--- testing {from_points,get_name,set_name,tag_names,points} ---\n")

  p30 = (4,3)
  p31 = (3,2)
  p32 = (2,3)

  ph3 = path.from_points((p30, ( p31, p32, ),), mp_fill, mp_jump)
  assert path.nelems(ph3) == 2
  mv30 = path.elem(ph3, 0)
  mv31 = path.elem(ph3, 1)
  assert move.is_jump(mv30)
  assert not move.is_jump(mv31)
  assert path.pini(ph3) == p30
  assert path.pfin(ph3) == p32
  
  ph3r = path.rev(ph3)
  assert path.elem(ph3r, 0) == move.rev(mv31)
  assert path.elem(ph3r, 1) == move.rev(mv30)
  assert path.pini(ph3r) == p32
  assert path.pfin(ph3r) == p30
  
  sys.stderr.write("points_cmp = %s\n" % str(path.points(ph3)))
  sys.stderr.write("points_exp = %s\n" % str((p30, p31, p32)))
  assert list(path.points(ph3)) == [p30, p31, p32]
  assert list(path.points(ph3r)) == [p32, p31, p30]

  path.validate(ph3)
  path.validate(ph3r)

  assert path.get_name(ph3) == "P?"
  path.set_name(ph3, "Tao", False)
  assert path.get_name(ph3) == "Tao"
  assert path.get_name(ph3r) == "~Tao"
  
  ph4 = path.from_points((p30, p32, ), mp_fill, mp_jump)

  sys.stderr.write("  ... applying {path.tag_names} ...\n")
  path.tag_names([ph3, ph4], "Tag.")
  assert path.get_name(ph3) == "Tag.Tao"
  assert path.get_name(ph3r) == "~Tag.Tao"
  assert path.get_name(ph4) == "Tag.P?"

  sys.stderr.write("\n")
  return
  # ----------------------------------------------------------------------
 
def test_links():
  sys.stderr.write("--- testing {clear_links,add_link,set_links,get_links,get_connecting_link,get_all_connecting_links} ---\n")
  
  OPHS = path_example.misc_H(mp_fill)
  sys.stderr.write("OPHS:\n")
  path.show_list(sys.stderr, "    ", OPHS, None, True, True)
  
  def makelink(p,q):
    # A slightly twisted path from {p} to {q}
    m = rn.mix(0.45, rn.add(p,q), 0.10, (sin(100*p[0]), cos(200*p[1])))
    link = path.from_points((p,m,q,), mp_cont, None)
    return link
    # ......................................................................
  
  def listeq(L0,L1):
    # True iff path lists {L0} and {L1} are the same 
    # considering {(ph,0)} and {ph} equal.
    L0r = tuple( (path.unpack(oph)) for oph in L0 )
    L1r = tuple( (path.unpack(oph)) for oph in L1 )
    return L0r == L1r
    # ......................................................................
  
  def patheq(oph0,oph1):
    # True iff paths {oph0} and {oph11} are the same 
    # considering {(ph,0)} and {ph} equal.
    return (path.unpack(oph0)) == (path.unpack(oph1))
    # ......................................................................
  
  def check_all_links(oph0, oph1, OLKS_exp):
    OLKS_cmp = path.get_all_connecting_links(oph0, oph1)
    OLKS_cmp = sorted(OLKS_cmp, key = path.pini)
    sys.stderr.write("OLKS_cmp = %s\n" % str([ path.get_name(olk) for olk in OLKS_cmp]))

    OLKS_exp = sorted(OLKS_exp, key = path.pini)
    sys.stderr.write("OLKS_exp = %s\n" % str([ path.get_name(olk) for olk in OLKS_exp]))

    assert listeq(OLKS_cmp, OLKS_exp)
    return
    # ....................................................................
  
  for oph in OPHS: 
    path.clear_links(oph)
    path.clear_links(path.rev(oph))
    
  # Links along the square, and diagonals:
  LKSA = []
  LKSB = []
  for iph, jph in (0,1), (1,2), (2,3), (3,0), (0,2), (1,3):
    
    linkA = makelink(path.pini(OPHS[iph]), path.pini(OPHS[jph]))
    path.set_name(linkA, "LA%d%d" % (iph,jph), False)
    LKSA.append(linkA)
    path.add_link(OPHS[iph], path.rev(linkA))
    path.add_link(OPHS[jph], linkA)
    
    linkB = makelink(path.pfin(OPHS[iph]), path.pfin(OPHS[jph]))
    path.set_name(linkB, "LB%d%d" % (iph,jph), False)
    LKSB.append(linkB)
    path.add_link(path.rev(OPHS[iph]), path.rev(linkB))
    path.add_link(path.rev(OPHS[jph]), linkB)
   
  sys.stderr.write("LKSA:\n")
  path.show_list(sys.stderr, "    ", LKSA, None, False, False)
  sys.stderr.write("LKSB:\n")
  path.show_list(sys.stderr, "    ", LKSB, None, False, False)

  # Check:
  L00a = path.get_links(OPHS[0])
  L00z = [ path.rev(LKSA[0]), LKSA[3], path.rev(LKSA[4]), ]
  sys.stderr.write("L00a:\n")
  path.show_list(sys.stderr, "    ", L00a, None, False, False)
  sys.stderr.write("L00z:\n")
  path.show_list(sys.stderr, "    ", L00z, None, False, False)
  assert listeq(L00a, L00z)
  
  L10a = path.get_links(OPHS[1])
  sys.stderr.write("L10a:\n")
  path.show_list(sys.stderr, "    ", L10a, None, False, False)
  L10z = [ LKSA[0], path.rev(LKSA[1]), path.rev(LKSA[5]), ]
  sys.stderr.write("L10z:\n")
  path.show_list(sys.stderr, "    ", L10z, None, False, False)
  assert listeq(L10a, L10z)

  assert listeq(path.get_links(path.rev(OPHS[0])), [ path.rev(LKSB[0]), LKSB[3], path.rev(LKSB[4]), ])
  assert listeq(path.get_links(path.rev(OPHS[1])), [ LKSB[0], path.rev(LKSB[1]), path.rev(LKSB[5]), ])
  
  connR0P1 = path.get_connecting_link(path.rev(OPHS[0]), OPHS[1])
  path.show(sys.stderr, "    connR0P1:\n", connR0P1, "\n", False, 0,0,0)
  assert patheq(connR0P1, LKSA[0])
  
  assert patheq(path.get_connecting_link(path.rev(OPHS[1]), OPHS[2]), LKSA[1])
  assert patheq(path.get_connecting_link(path.rev(OPHS[0]), OPHS[2]), LKSA[4])
  
  assert patheq(path.get_connecting_link(path.rev(OPHS[1]), OPHS[0]), path.rev(LKSA[0]))
  assert patheq(path.get_connecting_link(path.rev(OPHS[2]), OPHS[1]), path.rev(LKSA[1]))
  assert patheq(path.get_connecting_link(path.rev(OPHS[2]), OPHS[0]), path.rev(LKSA[4]))
  
  assert patheq(path.get_connecting_link(OPHS[0], path.rev(OPHS[1])), LKSB[0])
  assert patheq(path.get_connecting_link(OPHS[1], path.rev(OPHS[2])), LKSB[1])
  assert patheq(path.get_connecting_link(OPHS[0], path.rev(OPHS[2])), LKSB[4])
                                                               
  assert patheq(path.get_connecting_link(OPHS[1], path.rev(OPHS[0])), path.rev(LKSB[0]))
  assert patheq(path.get_connecting_link(OPHS[2], path.rev(OPHS[1])), path.rev(LKSB[1]))
  assert patheq(path.get_connecting_link(OPHS[2], path.rev(OPHS[0])), path.rev(LKSB[4]))
  
  check_all_links(OPHS[0], OPHS[1], [LKSA[0], LKSB[0],])
  check_all_links(OPHS[0], OPHS[2], [LKSA[4], LKSB[4],])
  check_all_links(OPHS[1], OPHS[2], [LKSA[1], LKSB[1],])
  
  # Check {set_links}:
  for ophk in OPHS:
    for ophr in ophk, path.rev(ophk):
      LKS = path.get_links(ophr)
      path.clear_links(ophr)
      assert len(path.get_links(ophr)) == 0
      path.set_links(ophr, LKS)
      assert tuple(path.get_links(ophr)) == tuple(LKS)
  
  sys.stderr.write("\n")
  return 
  # ----------------------------------------------------------------------

def test_groups():
  sys.stderr.write("--- testing {set_group,get_group} ---\n")
  
  ph = path_example.misc_B(mp_fill, mp_jump)
  
  assert path.get_group(ph) == None
  assert path.get_group(path.rev(ph)) == None

  for oph in ph, path.rev(ph):
    path.set_group(oph, 418)
    assert path.get_group(oph) == 418
    path.set_group(path.rev(oph), 4615)
    assert path.get_group(oph) == 4615

  sys.stderr.write("\n")
  return 
  # ----------------------------------------------------------------------

def test_separate_by_groups():
  sys.stderr.write("--- testing {separate_by_group} ---\n")
  
  OPHS, TRS, JMS = path_example.misc_D(mp_fill, mp_jump)
  nph = len(OPHS)

  # Assign the paths to three groups plus {None}:
  kgr = 0;
  igr_max = -1
  for oph in OPHS:
    if kgr == 3:
      path.set_group(oph, None)
    else:
      igr = 10*kgr
      path.set_group(oph, igr)
      igr_max = max(igr_max, igr)
    kgr = (kgr + 1) % 4
  assert igr_max >= 0  
  
  GRPHS, ngr = path.separate_by_group(OPHS)

  assert ngr == 3
  assert len(GRPHS) == igr_max + 1
  
  kgr = 0
  for igr in range(len(GRPHS)):
    GPS = GRPHS[igr]
    if igr == 10*kgr:
      assert len(GPS) != 0
      for kph in range(len(GPS)):
        assert path.get_group(GPS[kph]) == igr
        assert OPHS[kgr + 4*kph] == GPS[kph]
      kgr += 1
    else:
      assert len(GPS) == 0

  sys.stderr.write("\n")
  return 
  # ----------------------------------------------------------------------

def collect_links(OPHS, oph):
  # Returns a set of {Path} objects that are links leading to the paths
  # of {OPHS} and could be used as links in their concatenation {oph}.
  # Basically all links of all paths in {OPHS} up to the first
  # non-trivial one.
  OLKS = set()
  for ophi in OPHS:
    for olkij in path.get_links(ophi):
      # Add {olkij} but avoid 
      lkij, drij = path.unpack(olkij)
      OLKS.add((lkij, drij))
    if path.nelems(ophi) > 0: break

def check_concat_links(oph, OPHS):
  # The path {oph} must be the result of concatenating the paths in
  # {OPHS}. Checks if the input links of the paths {OPHS} have been
  # added as input links of {oph}, minus the applicable exceptions.
  
  p0 = path.pini(oph)
  
  # Collect the {Path} objects of input links of {OPHS}, up to first non-trivial path:
  LKS_oel = set()
  for iph in range(len(OPHS)):
    ophi = OPHS[iph]
    for olkij in path.get_links(ophi):
      assert path.pfin(olkij) == path.pini(ophi), "link does not connect to path"
      assert path.pini(ophi) == p0, "inconsistent pini(oph)"
      assert path.pini(olkij) != p0, "link is a closed loop"
      lkij, drij = path.unpack(olkij)
      LKS_oel.add(lkij)
    if path.nelems(ophi) > 0: break
  
  # Collect the {Path} objects of the input links of {oph}:
  LKS_oph = set((path.unpack(olk))[0] for olk in path.get_links(oph))
  
  # The latter must be a subset of the former:
  assert LKS_oph <= LKS_oel

  # Check whether the excluded links have been omitted for good reason:
  for lk in LKS_oel - LKS_oph:
    # Set {olk} to the oriented version of {lk} that goes into {oph}:
    if path.pfin(lk) == p0:
      olk = lk
    elif path.pini(lk) == p0:
      olk = path.rev(lk)
    else:
      assert False # Should have been caught earlier.

    # Link {lk} must have been excluded because it starts on {oph}:
    q = path.pini(olk)
    # We already checked that {q} is not {p0=pini(oph)}.  Check the other nodes of {oph}:
    ok = True
    for kmv in range(path.nelems(oph)):
       pk = move.pfin(path.elem(oph, kmv))
       if hacks.same_point(pk, q, 1.0e-5): ok = False; break
    if ok:
       sys.stderr.write("** link of factor path was improperly excluded:\n")
       path.show(sys.stderr, "  oph = ", oph, "\n", False, 0,0,0)
       path.show(sys.stderr, "  lk =  ", lk,  "\n", False, 0,0,0)
       assert False

  return
  # ----------------------------------------------------------------------
  
def test_concat():
  
  nph = 3
  TRS, PHS, LKS, CTS = raster_example.rasters_B(nph, (0.8, 0.6), mp_fill, mp_link)
  assert len(PHS) == nph
  OPHS = PHS
  tag = "rasters_B"
  do_test_concat(tag, OPHS)


  OEXS = path_example.misc_C(mp_fill, mp_link, mp_jump)
  OPHS = [ OEXS[0], OEXS[1], OEXS[2], path.rev(OEXS[4]), path.rev(OEXS[3]) ]
  tag = "path_misc_C"
  do_test_concat(tag, OPHS)

  sys.stderr.write("\n")
  return
  # ----------------------------------------------------------------------

def do_test_concat(tag, OPHS):

  sys.stderr.write("--- testing {concat} tag = %s ---\n" % tag)

  path.show_list(sys.stderr, "  ", OPHS, None, True, True)
  assert type(OPHS) is list
  assert isinstance(OPHS[0], path.Path)

  for use_links in False, True:
    phc = path.concat(OPHS, use_links, mp_jump)
    path.set_name(phc, "PP", False)
    sys.stderr.write("phc = %s\n" % str(phc))

    # General tests:

    assert path.pini(phc) == path.pini(OPHS[0])
    assert path.pfin(phc) == path.pfin(OPHS[-1])
    assert path.pini(path.rev(phc)) == path.pfin(OPHS[-1])
    assert path.pfin(path.rev(phc)) == path.pini(OPHS[0])
    nmv = path.nelems(phc)
    omvprev = None
    for imv in range(nmv):
      assert path.elem(path.rev(phc), imv) == move.rev(path.elem(phc, nmv-1-imv))
      omvk = path.elem(phc, imv)
      if omvprev != None:
        assert move.pini(omvk) == move.pfin(omvprev)
      assert path.find_move(phc, omvk) == imv
      assert path.find_move(phc, move.rev(omvk)) == imv
      assert path.find_move(path.rev(phc), omvk) == nmv-1-imv
      assert path.find_move(path.rev(phc), move.rev(omvk)) == nmv-1-imv
      omvprev = omvk

    # Check use of links in gaps:

    kmv = 0
    oph_prev = None
    for oph in OPHS:
      p_prev = None if oph_prev == None else path.pfin(oph_prev)
      p = path.pini(oph)
      if p_prev != None and p_prev != p:
        # A connector must have been inserted between {oph_prev} and {oph}:
        sys.stderr.write("checking connector %s -- %s\n" % (str(p_prev), str(p)))
        if use_links:
          olk = path.get_connecting_link(oph_prev, oph)
          # Should there be a link path here?
          should_exist = rn.dist(path.pfin(oph_prev), path.pini(oph)) <= 2.5
          assert (olk != None) == should_exist
        else:
          olk = None
        if olk != None:
          # The moves of {olk} must be in the path {phc}:
          for jmv in range(path.nelems(olk)):
            omvj = path.elem(olk, jmv)
            omvk = path.elem(phc, kmv)
            sys.stderr.write("link omvj = %s  %s -- %s\n" % (str(omvj), move.pini(omvj), move.pfin(omvj)))
            sys.stderr.write("link omvk = %s  %s -- %s\n" % (str(omvk), move.pini(omvk), move.pfin(omvk)))
            assert omvj == omvk
            kmv += 1
        else:
          # There must be a jump:
          omvk = path.elem(phc, kmv)
          sys.stderr.write("jump omvk = %s  %s -- %s\n" % (str(omvk), move.pini(omvk), move.pfin(omvk)))
          assert move.is_jump(omvk)
          assert move.pini(omvk) == path.pfin(oph_prev)
          assert move.pfin(omvk) == path.pini(oph)
          kmv += 1
      # The moves of {oph} must be in {phc}:
      for jmv in range(path.nelems(oph)):
        omvj = path.elem(oph, jmv)
        omvk = path.elem(phc, kmv)
        sys.stderr.write("path omvj = %s  %s -- %s\n" % (str(omvj), move.pini(omvj), move.pfin(omvj)))
        sys.stderr.write("path omvk = %s  %s -- %s\n" % (str(omvk), move.pini(omvk), move.pfin(omvk)))
        assert omvj == omvk 
        kmv += 1
      oph_prev = oph

    assert kmv == nmv

    # Check if links have been preserved:
    check_concat_links(phc, OPHS)
    check_concat_links(path.rev(phc), [path.rev(ophi) for ophi in reversed(OPHS)])

    path.validate(phc)
    path.validate(path.rev(phc))

    ctrace = pyx.color.rgb(1.000, 0.000, 0.700)
    rwd = 0.80
    wd_axes = 0.05*wd_fill;  # Axes of traces.
    grid = True
    deco = True

    fname = "tests/out/path_TST_concat_%s_ul%s" % (tag, ("T" if use_links else "F"))
    path.plot_to_files(fname, [phc,], [ctrace,], rwd, wd_axes, grid, deco)

  sys.stderr.write("\n")
  return
  # ----------------------------------------------------------------------

def test_find_nearest_move_endpoint():
  sys.stderr.write("--- testing {find_nearest_move_endpoint} ---\n")
  
  def test_path(oph):
    nmv = path.nelems(oph)
    PTS = [ move.pini(path.elem(oph,imv)) for imv in range(nmv) ]
    if nmv == 0:
      sys.stderr.write("  ... trivial path ...\n")
      kmax = 0
    elif path.pini(oph) == path.pfin(oph):
      sys.stderr.write("  ... closed path ...\n")
      kmax = nmv - 1
    else:
      sys.stderr.write("  ... open path ...\n")
      PTS.append(path.pfin(oph))
      kmax = nmv
    for kpt in range(len(PTS)):
      pt = rn.mix(1.0, PTS[kpt], 0.001, (sin(kpt), cos(kpt)))
      d = rn.dist(pt, PTS[kpt]) 
      kptx, dx = path.find_nearest_move_endpoint(oph, pt)
      assert kptx <= kmax
      ptx = move.pini(path.elem(oph, kptx)) if kptx < nmv else path.pfin(oph)
      # sys.stderr.write("kpt =  %d pt =  ( %20.16f %20.16f ) d =  %20.16f\n" % (kpt,pt[0],pt[1],d))
      # sys.stderr.write("kptx = %d ptx = ( %20.16f %20.16f ) dx = %20.16f\n" % (kptx,ptx[0],ptx[1],dx))
      assert kptx == kpt
      assert dx == rn.dist(pt, PTS[kpt])
    return 
    # ......................................................................
    
  OPHS = path_example.misc_J(mp_cont,mp_fill)
  
  test_path(OPHS[0])
  test_path(OPHS[1])

  sys.stderr.write("\n")
  return
  # ----------------------------------------------------------------------
  
def test_mean_projections():
  sys.stderr.write("--- testing {mean_projections} ---\n")
  
  ph = path_example.misc_A(mp_fill, mp_cont, mp_jump)
  assert isinstance(ph, path.Path)

  angle = pi/3
  xdir = rn.rotate2((1, 0), angle)
  ydir = rn.rotate2((0, 1), angle)
  
  m = rn.mix(0.5, path.pini(ph), 0.5, path.pfin(ph))
  xm_ex = rn.dot(m, xdir)
  ym_ex = rn.dot(m, ydir)
  
  for oph in ph, path.rev(ph):
    xm_cp, ym_cp = path.mean_projections(oph, xdir, ydir)
    assert abs(xm_cp - xm_ex) < 1.0e-14
    assert abs(ym_cp - ym_ex) < 1.0e-14

    xm_cp, ym_cp = path.mean_projections(oph, None, ydir)
    assert xm_cp == None
    assert abs(ym_cp - ym_ex) < 1.0e-14

    xm_cp, ym_cp = path.mean_projections(oph, xdir, None)
    assert abs(xm_cp - xm_ex) < 1.0e-14
    assert ym_cp == None
    
  sys.stderr.write("\n")
  return
  # ----------------------------------------------------------------------

def test_displace():
  sys.stderr.write("--- testing {displace} ---\n")

  ph = path_example.misc_A(mp_fill, mp_cont, mp_jump)
  assert isinstance(ph, path.Path)

  angr = pi/3
  vr = (5,7)
  phr =  path.displace(ph, angr, vr)
  nr = path.nelems(phr) 
  assert nr == 3
  for imv in range(nr):
    mvk = path.elem(ph, imv)
    pk, qk = move.endpoints(mvk)
    pkr = rn.add(rn.rotate2(pk, angr), vr)
    qkr = rn.add(rn.rotate2(qk, angr), vr)
    mvrk = path.elem(phr,imv)
    prk, qrk = move.endpoints(mvrk)
    assert rn.dist(prk, pkr) < 1.0e-8
    assert rn.dist(qrk, qkr) < 1.0e-8

  path.validate(phr)
  path.validate(path.rev(phr))

  sys.stderr.write("\n")
  return
  # ----------------------------------------------------------------------

def test_orient_unpack():
  sys.stderr.write("--- testing {spin, path.unpack} ---\n")

  def check_orient_unpack(ph):
    for dr in range(4):
      pht, drt = path.unpack(path.spin(ph, dr))
      assert pht == pht
      assert drt == dr % 2
      phb, drb = path.unpack(path.spin(path.rev(pht), dr))
      assert phb == pht
      assert drb == (dr + 1) % 2

      path.validate(pht)
      path.validate(path.rev(pht))

  PHS = path_example.misc_C(mp_fill, mp_link, mp_jump)
  assert type(PHS) is list
  assert isinstance(PHS[0], path.Path)
  pha = PHS[0]
  check_orient_unpack(pha)

  sys.stderr.write("\n")
  return
  # ----------------------------------------------------------------------

def test_name():
  sys.stderr.write("--- testing {set_name,get_name,tag_names} ---\n");
  
  def gname(var, oph):
    # Returns {get_name(oph)}, printing to {stderr}.
    xms = path.get_name(oph)
    sys.stderr.write("  get_name(%s) = %s\n" % (var,str(xms)))
    return xms
    # ............................................................

  p1 = (1,1)
  p2 = (2,3)
  p3 = (4,5)
  p4 = (6,7)
  pha = path.from_points(( (p1,), (p2,p3,p4), ), mp_fill,mp_jump) 
  xms = gname("pha", pha)
  assert xms == "P?"
  xms = gname("~pha", path.rev(pha))
  assert xms == "~P?"

  path.set_name(pha, "Tao", False)
  xms = gname("pha", pha)
  assert xms == "Tao"
  xms = gname("~pha", path.rev(pha))
  assert xms == "~Tao"

  sys.stderr.write("\n")
  return
  # ----------------------------------------------------------------------

def test_show():

  sys.stderr.write("--- testing {show,show_list} ---\n")

  OPHS, TRS, JMS = path_example.misc_D(mp_fill, mp_jump)
  nph = len(OPHS)
  for iph in range(nph): path.set_group(OPHS[iph], 8 + iph//2)
  
  sys.stderr.write("  ... {show} ...\n")
  wna = 5
  wnm = 5
  wgr = 3
  sys.stderr.write("\n")
  for kph in range(nph):
    moves = (kph % 2 == 1)
    for oph in OPHS[kph], path.rev(OPHS[kph]):
      path.show(sys.stderr, "   [ ", oph, " ]\n", moves, wna,wnm,wgr)
      wna = wna + 2
      wnm = max(wnm - 1, 0)
  sys.stderr.write("\n")
  
  sys.stderr.write("  ... {show_list} ...\n")
  path.show_list(sys.stderr, "    < ", OPHS, " >", True, True)

  sys.stderr.write("\n")
  return 
  # ----------------------------------------------------------------------

def do_plot_bbox(c, OPHS, dp, wdbox):
  # Plots the bounding box of oriented paths in the list {OPHS} displaced by {dp}.
  B = path.bbox(OPHS)
  hacks.plot_frame(c, B, dp, pyx.color.rgb.blue, wdbox, None, -wdbox/2)
  return
  # ----------------------------------------------------------------------

def test_contacts():
  sys.stderr.write("--- testing {get_contacts,set_contacts} ---\n")
  
  ph = path_example.misc_B(mp_fill, mp_jump)
  
  # ??? Should test with real contacts ???
  ct1 = "FOO"
  ct2 = "BAR"
  ct3 = "BAZ"
  
  for oph in ph, path.rev(ph):
    for isd in range(2):
      path.clear_contacts(oph, isd)
      assert tuple(path.get_contacts(oph, isd)) == ()
    path.add_contact(oph, 0, ct1)
    path.add_contact(oph, 0, ct2)
    path.add_contact(oph, 1, ct3)
    assert tuple(path.get_contacts(oph, 0)) == (ct1, ct2,)
    assert tuple(path.get_contacts(oph, 1)) == (ct3,)
  
  sys.stderr.write("\n")
  return 
  # ----------------------------------------------------------------------

def test_incremental():
  
  sys.stderr.write("--- testing {move_to, finish} ---\n")

  MVS = []

  PHS = []
  
  npt0 = 30; R0 =  4
  npt1 = 40; R1 = 10
  ctr = ( R1 + 1, R1 + 1 )
  p = None
  for ipt in range(npt0+1):
    a = 2*pi*ipt/npt0
    q = rn.mix(1, ctr, R0, ( cos(a), sin(a) + 0.15*sin(5*a) ))
    p = path.move_to(MVS, p, q, mp_cont, "C0")
  path.finish(PHS, MVS, "C0")
  assert len(MVS) == 0
  assert len(PHS) == 1
  path.validate(PHS[0])
  assert path.nelems(PHS[0]) == npt0
  
  p = None
  for ipt in range(npt1+1):
    a = 2*pi*ipt/npt1
    q = rn.mix(1, ctr, R1*(1 + 0.10*cos(5*a)), ( cos(a), sin(a) ))
    p = path.move_to(MVS, p, q, mp_cont, "C1")
  path.finish(PHS, MVS, "C1")
  assert len(MVS) == 0
  assert len(PHS) == 2
  path.validate(PHS[1])
  assert path.nelems(PHS[1]) == npt1
  
  CLRS = hacks.trace_colors(len(PHS), None)
  rwd = 0.80
  wd_axes = 0.05*wd_fill;  # Axes of traces.
  grid = True
  deco = True

  fname = "tests/out/path_TST_incremental"
  path.plot_to_files(fname, PHS, CLRS, rwd, wd_axes, grid, deco)
  
  sys.stderr.write("\n")
  return
  # ----------------------------------------------------------------------

def test_contours():
  sys.stderr.write("--- testing {compute_contour_nesting,inner_contours,outer_contours,shift_contour} ---\n")
  
  sys.stderr.write("  ... two nested circles ...\n")

  OCRS0 = path_example.contours_B(mp_cont)
  
  phc0 = OCRS0[0]
  phc1 = OCRS0[1]
  assert frozenset(path.inner_contours(phc0)) == frozenset()
  assert frozenset(path.inner_contours(phc1)) == frozenset((phc0,))
  assert path.outer_contour(phc0) == phc1
  assert path.outer_contour(phc1) == None
  
  phc0r = path.rev(phc0)
  phc0s = path.shift_contour(phc0r, 3)
  assert path.nelems(phc0r) == path.nelems(phc0)
  assert path.pini(phc0s) == move.pini(path.elem(phc0r, 3))
  assert path.pfin(phc0s) == path.pini(phc0s)
  nmv0 = path.nelems(phc0)
  for imv in range(nmv0):
    assert path.elem(phc0s, imv) == path.elem(phc0r, (imv + 3) % nmv0)

  sys.stderr.write("  ... ten paths, complex nesting ...\n")

  OCRS, PTSS = path_example.contours_A(mp_cont)
 
  ncr = len(OCRS)
  
  # Check for correct points:
  for icr in range(ncr):
    ocri = OCRS[icr]
    PTSi = PTSS[icr]
    PTSx = [ move.pini(path.elem(ocri, kmv)) for kmv in range(path.nelems(ocri)) ]
    assert tuple(PTSx) == tuple(PTSi)

  # Check containment relation:
  # Shorter names for the contours:
  assert ncr == 10
  A,B,C,D,E,F,G,H,I,J = OCRS

  assert set(path.inner_contours(A)) == set((C,D,))
  assert set(path.inner_contours(B)) == set((E,))
  assert set(path.inner_contours(C)) == set((F,G,))
  assert set(path.inner_contours(D)) == set((H,))
  assert set(path.inner_contours(E)) == set()
  assert set(path.inner_contours(F)) == set()
  assert set(path.inner_contours(G)) == set()
  assert set(path.inner_contours(H)) == set((I,J,))
  assert set(path.inner_contours(I)) == set()
  assert set(path.inner_contours(J)) == set()

  assert path.outer_contour(A) == None
  assert path.outer_contour(B) == None
  assert path.outer_contour(C) == A
  assert path.outer_contour(D) == A
  assert path.outer_contour(E) == B
  assert path.outer_contour(F) == C
  assert path.outer_contour(G) == C
  assert path.outer_contour(H) == D
  assert path.outer_contour(I) == H
  assert path.outer_contour(J) == H
  
  assert path.contour_nesting(A,A) == 0
  assert path.contour_nesting(A,B) == 0
  assert path.contour_nesting(A,C) == +1
  assert path.contour_nesting(A,D) == +1
  assert path.contour_nesting(A,E) == 0
  assert path.contour_nesting(A,F) == +1
  assert path.contour_nesting(A,G) == +1
  assert path.contour_nesting(A,H) == +1
  assert path.contour_nesting(A,I) == +1
  assert path.contour_nesting(A,J) == +1
  
  assert path.contour_nesting(C,A) == -1
  assert path.contour_nesting(C,B) == 0
  assert path.contour_nesting(C,C) == 0
  assert path.contour_nesting(C,D) == 0
  assert path.contour_nesting(C,E) == 0
  assert path.contour_nesting(C,F) == +1
  assert path.contour_nesting(C,G) == +1
  assert path.contour_nesting(C,H) == 0
  assert path.contour_nesting(C,I) == 0
  assert path.contour_nesting(C,J) == 0
  
  assert path.contour_nesting(D,A) == -1
  assert path.contour_nesting(D,B) == 0
  assert path.contour_nesting(D,C) == 0
  assert path.contour_nesting(D,D) == 0
  assert path.contour_nesting(D,E) == 0
  assert path.contour_nesting(D,F) == 0
  assert path.contour_nesting(D,G) == 0
  assert path.contour_nesting(D,H) == +1
  assert path.contour_nesting(D,I) == +1
  assert path.contour_nesting(D,J) == +1
  
  assert path.contour_nesting(H,A) == -1
  assert path.contour_nesting(H,B) == 0
  assert path.contour_nesting(H,C) == 0
  assert path.contour_nesting(H,D) == -1
  assert path.contour_nesting(H,E) == 0
  assert path.contour_nesting(H,F) == 0
  assert path.contour_nesting(H,G) == 0
  assert path.contour_nesting(H,H) == 0
  assert path.contour_nesting(H,I) == +1
  assert path.contour_nesting(H,J) == +1
  
  for ocr0 in OCRS:
    for ocr1 in OCRS:
      ne01 = path.contour_nesting(ocr0,ocr1)
      ne10 = path.contour_nesting(ocr1,ocr0)
      if ocr0 == ocr1:
        assert ne01 == 0 and ne10 == 0
      else:
        assert ne01 == -ne10
   
  CLRS = hacks.trace_colors(len(OCRS), None)
  rwd = 0.80
  wd_axes = 0.05*wd_fill;  # Axes of traces.
  grid = True
  deco = True

  fname = "tests/out/path_TST_contours"
  path.plot_to_files(fname, OCRS, CLRS, rwd, wd_axes, grid, deco)
  
  sys.stderr.write("\n")
  return
  # ----------------------------------------------------------------------
 
def test_plot_basic():
  # Generates a plots showing the same set of oriented paths with
  # various options, plotting with {plot_standard} and mutiple calls of
  # {plot_layer}. Puts all in a single EPS figure. The file name will be
  # tests/out/path_TST_plot_basic.{ext}"

  sys.stderr.write("--- testing {bbox,plot_layer,plot_standard} ---\n")

  OPHS = path_example.misc_C(mp_fill, mp_link, mp_jump)
  nph = len(OPHS)

  CLRS = hacks.trace_colors(nph, None)
  nclr = len(CLRS)

  def pick_colors(kph):
    # Returns the colors for trace sausages and axes of path {OPHS[kph]}.
    if nclr == 1:
      ctraces = CLRS[0]
    else:
      ctraces = CLRS[kph]
    caxes =   pyx.color.rgb(0.6*ctraces.r, 0.6*ctraces.g, 0.6*ctraces.b) # Color of trace axis, dots, arrow.
    return ctraces, caxes

  # Dimensions
  rwd = 0.80; # Trace sausage width.
  wd_axes =   0.05*wd_fill;  # Axes of traces.
  wd_dots =   2.5*wd_axes;  # Dots at ends of moves.
  sz_arrows = 8*wd_axes      # Size of arrows. 
  wd_box = 0.5*wd_axes

  def do_test_plot_layer(c, dp, axes, dots, arrows, matter):
    # Plots on the {pyx} canvas {c} the oriented paths in the list {OPHS}
    # displaced by {dp} using multiple calls of {plot_layer}.

    # Colors:
    cmatter = pyx.color.rgb(0.800, 0.250, 0.750)  # Material footprint color.
    cjumps =  pyx.color.rgb.black  # Color of jump axes, dots, arrows.

    # Dimensions relative to nominal trace widths:
    rwd_matter = 1.13; # Material footprint.

    for layer in range(4):
      for kph in range(nph):
        oph = OPHS[kph]
        ctraces, caxes = pick_colors(kph)
        if layer == 0:
          # plots the estimate of actual material:
          path.plot_layer\
            ( c, oph, dp, jmp=False, \
              clr=cmatter, rwd=rwd_matter, wd=0, 
              dashed=False, wd_dots=0, sz_arrows=None
            )
        elif layer == 1:
          # plots the nominal trace material:
          path.plot_layer \
            ( c, oph, dp, jmp=False,
              clr=ctraces, rwd=rwd, wd=0, 
              dashed=False, wd_dots=0, sz_arrows=None
            )
        elif layer == 2:
          # Plots the trace axis:
          t_wd = wd_axes if axes else 0
          t_wd_dots = wd_dots if dots else 0
          t_sz_arrows = sz_arrows if arrows else 0
          path.plot_layer \
            ( c, oph, dp, jmp=False, 
              clr=caxes, rwd=0, wd=t_wd, 
              dashed=False, wd_dots=t_wd_dots, sz_arrows=t_sz_arrows
            )
        elif layer == 3:
          # Plots the jump:
          j_wd = wd_axes;
          j_wd_dots = wd_dots;
          j_sz_arrows = sz_arrows
          path.plot_layer \
            ( c, oph, dp, jmp=True, 
              clr=cjumps, rwd=0, wd=j_wd, 
              dashed=True, wd_dots=j_wd_dots, sz_arrows=j_sz_arrows
            )
        if layer == 0:
          # Plot the path's bounding box over the matter shadow:
          wd_box = 0.5*wd_axes
          do_plot_bbox(c, [oph,], dp, wd_box)    
    return

  def do_test_plot_standard(c, dp, axes, dots, arrows, matter):
    # Plots the oriented paths in the list {OPHS} shifted by {dp} using
    # {path.plot_standard} with no layers.

    path.plot_standard(c, OPHS, dp, None, CLRS, rwd, wd_axes, axes, dots, arrows, matter)
    # Plot the path's bounding box over everything:
    do_plot_bbox(c, OPHS, dp, wd_box)  
    return

  B = path.bbox(OPHS)
  
  dp = None
  
  c, szx, szy = hacks.make_canvas(hacks.round_box(B,0.5), dp, None, True, True, 4, 4)

  # Plot with various options:
  Y = 0
  matter = True
  for axes in False, True:
    for dots in False, True:
      for arrows in False, True:

        title = "ax%d dt%d ar%d" % (int(axes),int(dots), int(arrows))

        # Plot with {plot_standard}:
        dp0 = ((2*arrows+0)*szx, Y )
        pt0 = rn.add(dp0, (0.5, 0.5))
        c.text(pt0[0], pt0[1], "STD " + title)
        do_test_plot_standard(c, dp0, axes, dots, arrows, matter)

        # Plot with {plot_layer}:
        dp1 = ((2*arrows+1)*szx, Y)
        pt1 = rn.add(dp1, (0.5, 0.5))
        c.text(pt1[0], pt1[1], "BYL " + title)
        do_test_plot_layer(c, dp1, axes, dots, arrows, matter)

      Y = Y + szy

  hacks.write_plot(c, "tests/out/path_TST_plot_basic")

  sys.stderr.write("\n")
  return
  # ----------------------------------------------------------------------

def test_plot_single():
  # Generates a plots showing the same set of oriented paths with
  # various options, plotting with {plot_standard} and mutiple calls of
  # {plot_layer}. Puts all in a single EPS figure. The file name will be
  # tests/out/path_TST_plot_basic.{ext}"

  sys.stderr.write("--- testing {plot_single} ---\n")

  oph = path_example.misc_B(mp_fill, mp_jump)

  B = path.bbox([oph,])
  
  dp = (0,0)

  c, szx, szy = hacks.make_canvas(hacks.round_box(B,0.5), dp, None, True, True, 1, 2)

  clr = pyx.color.rgb(0.850, 0.050, 0.000)

  dpk = dp
  path.plot_single(c, oph, dpk, False, clr); dpk = rn.add(dpk, (0, szy))
  path.plot_single(c, oph, dpk, True, None); dpk = rn.add(dpk, (0, szy))

  hacks.write_plot(c, "tests/out/path_TST_plot_single")

  sys.stderr.write("\n")
  return
  # ----------------------------------------------------------------------

def test_plot_to_files(grid, deco):
  # Generates a plots showing a set of oriented paths plotted with {plot_to_file}.
  # The file name will be tests/out/path_TST_plot_to_files_deco{deco}.{ext}"

  sys.stderr.write("--- testing {plot_to_files} deco = %s---\n" % deco)
  
  tag = "grid%s_deco%s" % ("FT"[grid],"FT"[deco])

  OPHS, TRS, JMS = path_example.misc_E(mp_fill, mp_jump)
  nph = len(OPHS)

  CLRS = hacks.trace_colors(nph, None)
  rwd = 0.80
  wd_axes = 0.05*wd_fill;  # Axes of traces.

  fname = "tests/out/path_TST_plot_to_files_" + tag
  path.plot_to_files(fname, OPHS, CLRS, rwd, wd_axes, grid, deco)

  sys.stderr.write("\n")
  return 
  # ----------------------------------------------------------------------

def test_fabtime_tini_tfin():
  sys.stderr.write("--- testing {fabtime, path.tini, path.tfin} ---\n")

  def check_path_times(label, ph):
    sys.stderr.write("  ... testing times of %s ...\n" % label)
    path.validate(ph)
    phr = path.rev(ph)
    tpha = path.fabtime(ph)
    tphar = path.fabtime(phr)
    assert abs(tpha - tphar) < 1.0e-8
    nmv = path.nelems(ph)
    for imv in range(nmv):
      omvk = path.elem(ph,imv)
      mvk, drk = move.unpack(omvk)
      tinik = path.tini(ph, imv)
      tfink = path.tfin(ph, imv)
      # sys.stderr.write("tini[%d] = %12.8f" % (imv, tinik))
      # sys.stderr.write("  tfin[%d] = %12.8f\n" % (imv, tfink))
      omvkr = path.elem(phr, nmv-1-imv)
      mvkr, drkr = move.unpack(omvkr)
      assert mvk == mvkr
      assert drk == 1-drkr
      tinikr = path.tini(path.rev(ph), nmv-1-imv)
      tfinkr = path.tfin(path.rev(ph), nmv-1-imv)
      assert abs((tfink + tinikr) - tpha) < 1.0e-8
      assert abs((tinik + tfinkr) - tpha) < 1.0e-8
      if move.is_jump(omvk):
        assert (tfink - tinik) >= move.fabtime(omvk) - 1.0e-8
      else:
        assert abs((tfink - tinik) - move.fabtime(omvk)) < 1.0e-8
    tphb = path.tfin(ph, nmv-1)
    tphbr = path.tfin(phr, nmv-1)
    # sys.stderr.write("fabtime(ph) = %12.8f tphb = %12.8f\n" % (tpha, tphb))
    assert abs(tphb - tphbr) < 1.0e-8
    assert abs(tpha - tphb) < 1.0e-8

    assert abs(path.tini(ph,nmv) - tpha) < 1.0e-8
    assert path.tfin(ph,-1) == 0

  pa1 = (1,1)
  pa2 = (2,2)
  pa3 = (3,1)
  pa4 = (4,3)
  pa5 = (5,1)
  pa6 = (6,2)
  pa7 = (7,1)

  pha = path.from_points(((pa1, pa2, pa3,), (pa4, pa5,), (pa6, pa7,)), mp_fill, mp_jump)

  check_path_times("pha", pha)
  check_path_times("rev(pha)", path.rev(pha))
  
  # Checking that two jumps take longer than one:
  q1 = (2,2)
  q2 = (3,2)
  q3 = (4,2)
  
  oph1 = path.from_points(((q1,), (q2,), (q3,)), mp_fill, mp_jump)
  oph2 = path.from_points(((q1,), (q3,)), mp_fill, mp_jump)
  
  sys.stderr.write("two jumps: fabtime = %10.6f\n" % (path.fabtime(oph1)))
  sys.stderr.write("one jump:  fabtime = %10.6f\n" % (path.fabtime(oph2)))
  
  sys.stderr.write("  ... testing (non) convexity of fabtime ...\n")
  
  pf0 = (1.0, 1.0)
  pf1 = (1.1, 1.0)
  pf2 = (1.5, 1.0)
  pf3 = (2.0, 1.0)
  
  ophfa = path.from_points((pf0, pf1, pf3), mp_fill, None);  tfa = path.fabtime(ophfa)
  ophfb = path.from_points((pf0, pf2, pf3), mp_fill, None);  tfb = path.fabtime(ophfb)
  sys.stderr.write("tfa = %.4f  tfb = %.4f\n" % (tfa,tfb))
  assert tfa < tfb

  sys.stderr.write("\n")
  return
  # ----------------------------------------------------------------------

def test_check_links_and_contacts():
  sys.stderr.write("--- testing {check_links,check_contacts} ---\n")
  cols = 5
  rows = 2
  nx = 3
  ny = 20
  islands = True
  mp_link = mp_cont
  OCRS, OPHS, OLKS, CTS = raster_example.patch_array(cols, rows, nx,ny, islands, mp_cont, mp_fill, mp_link)
  
  sys.stderr.write("### CONTOUR PATHS (%d)\n" % len(OCRS))
  path.show_list(sys.stderr, "  ", OCRS, None, True, True)
  sys.stderr.write("\n")
  
  sys.stderr.write("### RASTER PATHS (%d)\n" % len(OPHS))
  path.show_list(sys.stderr, "  ", OPHS, None, True, False)
  sys.stderr.write("\n")
  
  sys.stderr.write("### LINK PATHS (%d)\n" % len(OLKS))
  path.show_list(sys.stderr, "  ", OLKS, None, True, True)
  sys.stderr.write("\n")
  
  sys.stderr.write("### CONTACTS (%d)\n" % len(CTS))
  contact.show_list(sys.stderr, "", CTS, None)
  sys.stderr.write("\n")

  ydir = (0,1) 
  path.check_links(OPHS,OLKS)
  contact.check_side_paths(CTS,OPHS,ydir)

  sys.stderr.write("\n")
  return
  # ----------------------------------------------------------------------

def test_simplify():
  sys.stderr.write("--- testing {simplify} ---\n")
  
  tag = "simplify"

  ctrace = pyx.color.rgb(1.000, 0.700, 0.000)
  rwd = 0.80
  wd_axes = 0.05*wd_fill;  # Axes of traces.
  grid = True
  deco = True

  nmv = 10
  step = 1.5*wd_fill
  pert = 0.50*wd_fill
  oph = path_example.misc_K(nmv,step,pert,mp_cont,mp_fill,mp_link,mp_jump)
  move_parms.show(sys.stderr, "mp_jump = ", mp_jump, "\n")

  path.show(sys.stderr, "oph = ", oph, "\n", False, 0,0,0)
  assert path.nelems(oph) == 4*nmv

  fname = "tests/out/path_TST_plot_to_files_" + tag + "oph"
  path.plot_to_files(fname, [oph,], None, rwd, wd_axes, grid, deco)

  tol = 2.1*pert
  
  osp = path.simplify(oph, tol, mp_jump)
  path.show(sys.stderr, "osp = ", osp, "\n", True, 0,0,0)

  fname = "tests/out/path_TST_plot_to_files_" + tag + "osp"
  path.plot_to_files(fname, [osp,], None, rwd, wd_axes, grid, deco)
  
  assert path.nelems(osp) == 7
  assert path.pini(osp) == path.pini(oph)
  assert path.pfin(osp) == path.pfin(oph)
  
  omv0 = path.elem(osp, 0); move.show(sys.stderr, "omv0 = ", omv0, "\n", 0)
  omv1 = path.elem(osp, 1); move.show(sys.stderr, "omv1 = ", omv1, "\n", 0)
  omv2 = path.elem(osp, 2); move.show(sys.stderr, "omv2 = ", omv2, "\n", 0)
  omv3 = path.elem(osp, 3); move.show(sys.stderr, "omv3 = ", omv3, "\n", 0)
  omv4 = path.elem(osp, 4); move.show(sys.stderr, "omv4 = ", omv4, "\n", 0)
  omv5 = path.elem(osp, 5); move.show(sys.stderr, "omv5 = ", omv5, "\n", 0)
  omv6 = path.elem(osp, 6); move.show(sys.stderr, "omv6 = ", omv6, "\n", 0)

  assert move.parameters(omv0) == mp_fill
  assert move.pini(omv0) == path.pini(oph)
  assert move.pfin(omv0) == move.pini(path.elem(oph, nmv))
  
  assert move.parameters(omv1) == mp_fill
  assert move.pfin(omv1) == move.pini(path.elem(oph, nmv+nmv//2))
  
  assert move.parameters(omv2) == mp_cont
  assert move.pfin(omv2) == move.pini(path.elem(oph, 2*nmv))
  
  assert move.parameters(omv3) == mp_fill
  assert move.pfin(omv3) == move.pini(path.elem(oph, 2*nmv + nmv//2))
  
  assert move.parameters(omv4) == mp_jump
  assert move.pfin(omv4) == move.pini(path.elem(oph, 2*nmv + nmv//2 + 1))
  
  assert move.parameters(omv5) == mp_fill
  assert move.pfin(omv5) == move.pini(path.elem(oph, 3*nmv))
  
  assert move.parameters(omv6) == mp_jump
  assert move.pfin(omv6) == path.pfin(oph)

  sys.stderr.write("\n")
  return
  # ----------------------------------------------------------------------

def test_find_midline_points():
  sys.stderr.write("--- testing {find_point_by_lenth,find_nearest_midline_point} ---\n")
  
  PTS = [ (1,1), (5,1), (5,5), (3,5), (2,3), (1,5), (1,3) ]
  oph = path.from_points(PTS, mp_fill, None)
  path.set_name(oph, "P", True)
  nmv = path.nelems(oph)
  
  sys.stderr.write("  ... testing {find_midline_point_by_length,find_nearest_midline_point} ---\n")
  
  imva, qa = path.find_midline_point_by_length(oph, 6.0)
  sys.stderr.write("    a: i = %d q = ( %10.8f %10.8f )\n" % (imva, qa[0], qa[1]))
  assert imva == 1
  assert rn.dist(qa, (5,3)) < 1.0e-8
  
  imvb, qb = path.find_midline_point_by_length(path.rev(oph), 0.5)
  sys.stderr.write("    b: i = %d q = ( %10.8f %10.8f )\n" % (imvb, qb[0], qb[1]))
  assert imvb == 0
  assert rn.dist(qb, (1,3.5)) < 1.0e-8
  
  sys.stderr.write("  ... testing {find_nearest_midline_point} ---\n")
  
  def do_test_fnmp(tag, ophr, p, q_exp, t_exp, imv_exp_a, imv_exp_b):
    sys.stderr.write("    test '%s'\n" % tag)
    q_cmp, imv_cmp, t_cmp, d_cmp = path.find_nearest_midline_point(ophr, p)
    
    sys.stderr.write("    cmp: i = %d   q = ( %10.8f %10.8f ) t = %10.8f d = %10.8f\n" % (imv_exp_a, q_cmp[0], q_cmp[1], t_cmp, d_cmp))
    
    d_exp = rn.dist(p, q_exp)
    sys.stderr.write("    exp: i =")
    if imv_cmp != imv_exp_b: sys.stderr.write(" %d" % imv_exp_a)
    if imv_exp_a != imv_exp_b and imv_cmp != imv_exp_a: 
      sys.stderr.write(" %d" % imv_exp_b)
    else:
      sys.stderr.write("  ")
    sys.stderr.write(" q = ( %10.8f %10.8f ) t = %10.8f d = %10.8f\n" % (q_exp[0], q_exp[1], t_exp, d_exp))

    assert imv_cmp == imv_exp_a or imv_cmp == imv_exp_b
    assert abs(t_cmp - t_exp) < 1.0e-8
    assert rn.dist(q_cmp, q_exp) < 1.0e-8
    assert rn.dist(q_cmp, p) == d_cmp
    return 
    # ....................................................................
  
  do_test_fnmp("c", oph,           (6.0,2.0), (5.0, 2.0), 5.0,                     1,1)
  
  do_test_fnmp("d", oph,           (0.0,6.0), (1.0,5.0), 10.0 + 2*hypot(1.0,2.0),  4,5)
  
  do_test_fnmp("e", path.rev(oph), (6.0,2.0), (5.0, 2.0), 7.0 + 2*hypot(1.0,2.0),  4,4)
  
  do_test_fnmp("f", oph,           (0.0,1.99), (1.0,1.0), 0.0,                     0,0)
  
  do_test_fnmp("g", oph,           (0.0,2.01), (1.0,3.0), 12.0 + 2*hypot(1.0,2.0), 5,6)
  
  sys.stderr.write("\n")
  return
  # ----------------------------------------------------------------------

def test_extract_trim():
  sys.stderr.write("--- testing {extract_section,trim} ---\n")
  
  PTS = [ (1,1), (5,1), (5,5), (3,5), (2,3), (1,5), (1,3) ]
  oph = path.from_points(PTS, mp_fill, None)
  path.set_name(oph, "P", True)
  nmv = path.nelems(oph)
  
  sys.stderr.write("  ... testing {extract_section} ---\n")
  
  def do_test_exsec(tit, ophr, imv0,q0, imv1,q1, sec_exp):
    sys.stderr.write("    test '%s'\n" % tit)
    sec_cmp = path.extract_section(ophr, imv0,q0, imv1,q1)
     
    path.set_name(ophr, "XC", False)
    path.show(sys.stderr, "    cmp = ", sec_cmp, "\n", True, 0,0,0)
    OMVS_cmp = [ path.elem(sec_cmp, k) for k in range(path.nelems(sec_cmp)) ]
    move.show_list(sys.stderr, "      ", OMVS_cmp, None)
    
    path.set_name(sec_exp, "XE", True)
    path.show(sys.stderr, "    exp = ", sec_exp, "\n", True, 0,0,0)
    OMVS_exp = [ path.elem(sec_exp, k) for k in range(path.nelems(sec_exp)) ]
    move.show_list(sys.stderr, "      ", OMVS_exp, None)

    path.compare(sec_cmp, sec_exp, 1.0e-8, True)
    return
    # ....................................................................
    
  tita = "whole path"
  ophra = oph; imv0a = -1; q0a = path.pini(oph); imv1a = nmv; q1a = path.pfin(oph)
  sec_expa = oph
  do_test_exsec(tita, ophra, imv0a,q0a, imv1a,q1a, sec_expa)
    
  titb = "trivial path (vertex)"
  ophrb = oph; imv0b = 1; q0b = (5,5); imv1b = 2; q1b = q0b
  sec_expb = path.make_trivial(q0b)
  do_test_exsec(titb, ophrb, imv0b,q0b, imv1b,q1b, sec_expb)
    
  titc = "trivial path (mid-move)"
  ophrc = oph; imv0c = 1; q0c = (5,3); imv1c = 1; q1c = q0c
  sec_expc = path.make_trivial(q0c)
  do_test_exsec(titc, ophrc, imv0c,q0c, imv1c,q1c, sec_expc)
    
  titd = "piece of single move"
  ophrd = oph; imv0d = 1; q0d = (5,2); imv1d = 1; q1d = (5,4)
  sec_expd = path.from_points((q0d, q1d), mp_fill, None) 
  do_test_exsec(titd, ophrd, imv0d,q0d, imv1d,q1d, sec_expd)
    
  tite = "pieces of two consecutive moves"
  ophre = oph; imv0e = 1; q0e = (5,2); imv1e = 2; q1e = (4,5)
  sec_expe = path.from_points((q0e, (5,5), q1e), mp_fill, None) 
  do_test_exsec(tite, ophre, imv0e,q0e, imv1e,q1e, sec_expe)

  titf = "multiple moves"
  ophrf = oph; imv0f = 1; q0f = (5,2); imv1f = 5; q1f = (1,4)
  sec_expf = path.from_points((q0f, (5,5), (3,5), (2,3), (1,5), q1f), mp_fill, None) 
  do_test_exsec(titf, ophrf, imv0f,q0f, imv1f,q1f, sec_expf)
    
  titg = "reverse path"
  ophrg = path.rev(oph); imv0g = 0; q0g = (1,4); imv1g = 4; q1g = (5,2)
  sec_expg = path.from_points((q0g, (1,5), (2,3), (3,5), (5,5), q1g), mp_fill, None) 
  do_test_exsec(titg, ophrg, imv0g,q0g, imv1g,q1g, sec_expg)
    
  sys.stderr.write("  ... testing {trim} ---\n")

  otr_cmp = path.trim(oph, 9, 0.5)
  otr_exp = path.from_points([ (4,5), (3,5), (2,3), (1,5), (1,3.5) ], mp_fill, None)
  nmv = path.nelems(otr_exp)
  assert path.nelems(otr_cmp) == nmv
  for i in range(nmv+1):
    if i < nmv:
      omv_exp = path.elem(otr_exp, i); qi_exp = move.pini(omv_exp)
      omv_cmp = path.elem(otr_cmp, i); qi_cmp = move.pini(omv_cmp)
    else:      
      qi_exp = path.pfin(otr_exp)
      qi_cmp = path.pfin(otr_cmp)
    sys.stderr.write("    i = %d" % i)
    sys.stderr.write(" cmp = ( %10.8f %10.8f )" % (qi_cmp[0], qi_cmp[1]))
    sys.stderr.write(" exp = ( %10.8f %10.8f )\n" % (qi_exp[0], qi_exp[1]))
    assert rn.dist(qi_cmp, qi_exp) < 1.0e-8

  sys.stderr.write("\n")
  return
  # ----------------------------------------------------------------------

test_show()
test_extract_trim()
test_find_midline_points()

test_fabtime_tini_tfin()
test_simplify()
test_concat()

test_check_links_and_contacts()
test_groups()
test_separate_by_groups()

test_plot_single()
test_incremental()
test_contours()

test_make_trivial()
test_from_move()
test_from_moves() 
test_from_points()
test_links()

test_displace()
test_orient_unpack()
test_find_nearest_move_endpoint()
test_mean_projections()
test_name()
test_plot_basic()
for grid in False, True:
  for deco in False, True:
    test_plot_to_files(grid, deco)

