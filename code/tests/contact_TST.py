#! /bin/usr/python3
# Test program for module {contact}
# Last edited on 2021-11-05 02:49:07 by stolfi

import contact
import contact_example
import move
import move_parms
import path
import block
import palette
import path_example
import raster_example
import block_example
import hacks
import job_parms
import rn
import pyx
import sys
from math import sqrt, sin, cos, floor, ceil, pi, nan, inf

parms = job_parms.typical_js()
parms['solid_raster_width'] = 1.00
parms['contour_trace_width'] = 0.50

mp_jump = move_parms.make_for_jumps(parms)
mp_cont = move_parms.make_for_contours(parms)
mp_fill = move_parms.make_for_fillings(parms)

wdf = move_parms.width(mp_fill)
wdc = move_parms.width(mp_cont)

def test_basic():

  sys.stderr.write("--- testing {make,from_moves,side,tcool_limit} ---\n")
  
  ya = 1
  yb = ya + (wdf+wdc)/2
  yc = yb + wdc
  
  sys.stderr.write("traces ya = %10.8f  yb = %10.8f  yc = %10.8f\n" % (ya,yb,yc))
  
  eps = 0.02
  
  pa0 = (  0, ya )
  qa0 = (  2, ya )
  
  pb0 = (  1, yb + 0.1*eps )
  qb0 = (  8, yb - 0.3*eps )
  pb1 = ( 10, yb )
  
  pc0 = (  0, yc )
  qc0 = (  2, yc )

  tra0 = move.make(pa0, qa0, mp_fill)
  trb0 = move.make(pb0, qb0, mp_cont)
  trc0 = move.make(pc0, qc0, mp_cont)
  
  jm_qb0_pb1 = move.make(qb0, pb1, mp_jump)
 
  def check_contact(name, ct, sides_exp, ends_exp, tclim_exp):
    sys.stderr.write("checking contact %s\n" % name)
    assert ct != None
    assert isinstance(ct, contact.Contact)
    assert contact.side_move(ct, 0) == sides_exp[0]
    assert contact.side_move(ct, 1) == sides_exp[1]
    assert contact.side_moves(ct) == (sides_exp[0], sides_exp[1])
    ends_cmp = contact.endpoints(ct)
    if ends_cmp[0][0] > ends_cmp[1][0]: ends_cmp = (ends_cmp[1], ends_cmp[0])
    sys.stderr.write("  endpoints cmp = ( %.9f %.9f ) ( %.9f %.9f )" % (ends_cmp[0]+ends_cmp[1]))
    sys.stderr.write(" length = %.9f\n" % rn.dist(ends_cmp[0], ends_cmp[1]))
    sys.stderr.write("  endpoints exp = ( %.9f %.9f ) ( %.9f %.9f )\n" % (ends_exp[0]+ends_exp[1]))
    assert rn.dist(ends_cmp[0], ends_exp[0]) < eps
    assert rn.dist(ends_cmp[1], ends_exp[1]) < eps
    assert contact.tcool_limit(ct) == tclim_exp
    return 
    # ....................................................................
    
  yab = ya + wdf/2
  ybc = yb + wdc/2
  
  eab0 = (1, yab)
  eab1 = (2, yab)
  eab_exp = ( eab0, eab1 )

  tclim_exp = 17.0
  
  sys.stderr.write("contacts (exp) yab = %10.8f  ybc = %10.8f\n" % (yab,ybc))

  ctX = contact.make(eab0, eab1, tra0, trb0)
  contact.set_tcool_limit(ctX, tclim_exp)
  check_contact("ctX", ctX, (tra0,trb0), eab_exp, tclim_exp)

  e0 = move.pini(tra0)
  e1 = move.pfin(tra0)
  mvdira0, mvlena0 = rn.dir(rn.sub(e1, e0))

  tol = 0.20*wdf  # Tolerance for overlaps, tilts, etc.

  for mvdir in mvdira0, None:
    nsu = 0 # Number of times {mkcontact} succeeded.
    for otra0 in tra0, move.rev(tra0):
      for otrb0 in trb0, move.rev(trb0):
        ctA = contact.from_moves(otra0, otrb0, mvdir, 0.9, 0.49, tol)
        if ctA != None:
          tclim_exp += 1
          contact.set_tcool_limit(ctA, tclim_exp)
          check_contact("ctA", ctA, (tra0,trb0), eab_exp, tclim_exp)
          nsu += 1
    assert nsu == 0 or nsu == 4, "depends on orientations"
 
  # Check length limits:
  mvdir = None # Let the prodeure determine it.
  ctAx = contact.from_moves(tra0, trb0, mvdir, 1.1, 0.00, tol)
  if ctAx != None: 
    contact.show(sys.stderr, "ctAx = ", ctAx, "\n", 0)
  assert ctAx == None
  
  ctAy = contact.from_moves(tra0, trb0, mvdir, 0.0, 0.51, tol)
  assert ctAy == None

  ctB = contact.from_moves(tra0, trc0, mvdir, 0, 0, tol)
  assert ctB == None
 
  ctC = contact.from_moves(tra0, jm_qb0_pb1, mvdir, 0, 0, tol)
  assert ctC == None
  
  return
  # ----------------------------------------------------------------------

def test_names():

  sys.stderr.write("--- testing {has_name,set_name,get_name,tag_names} ---\n")
  
  ya = 1
  yb = ya + (wdf+wdc)/2
  yc = yb + wdc
  
  mvdir = (1,0)

  pa0 = (  0, ya )
  qa0 = (  2, ya )
  
  pb0 = (  1, yb )
  qb0 = (  8, yb )
  
  pc0 = (  0, yc )
  qc0 = (  2, yc )

  tra0 = move.make(pa0, qa0, mp_fill); move.set_name("Tra0")
  trb0 = move.make(pb0, qb0, mp_cont); move.set_name("Trb0")
  trc0 = move.make(pc0, qc0, mp_cont); move.set_name("Trc0")
  
  tol = 0.20*wdf  # Tolerance for overlaps, tilts, etc.

  ctA = contact.from_moves(tra0, trb0, mvdir, 0.9, 0.49, tol)
  assert ctA != None
  assert not contact.has_name(ctA)
  nameA_exp = "C(%s:%s)" % (move.get_name(tra0), move.get_name(trb0))
  assert contact.has_name(ctA)
  assert contact.get_name(ctA) == nameA_exp

  ctD = contact.from_moves(trb0, trc0, mvdir, 0.9, 0.49, tol)
  assert ctD != None
  assert not contact.has_name(ctD)
  nameD_exp = "C(%s:%s)" % (move.get_name(trb0), move.get_name(trc0))
  assert contact.has_name(ctD)
  assert contact.get_name(ctD) == nameD_exp

  contact.set_name(ctA, "Close")
  assert contact.get_name(ctA) == "Close"
  
  sys.stderr.write("applying {contact.tag_names}:\n")
  contact.tag_names([ctA, ctD], "Tag.")
  assert contact.get_name(ctA) == "Tag.Close"
  assert contact.get_name(ctD) == "Tag." + nameD_exp
  return
  # ----------------------------------------------------------------------

def test_more_makes():

  sys.stderr.write("--- testing {from_moves,from_move_lists,from_paths,from_blocks} ---\n")

  BCS,PHS,TRS0,TRS1 = block_example.misc_G(mp_cont, mp_fill, mp_jump)
  
  sys.stderr.write("  ... traces of block 0 ...\n")
  move.show_list(sys.stderr, "    ", TRS0, None)

  sys.stderr.write("  ... traces of block 1 ...\n")
  move.show_list(sys.stderr, "    ", TRS1, None)
  
  sys.stderr.write("  ... paths ...\n")
  path.show_list(sys.stderr, "    ", PHS, None, True, True)
  
  sys.stderr.write("  ... blocks ...\n")
  block.show_list(sys.stderr, "    ", BCS, None, True)
  
  ph0 = PHS[0]; nmv0 = path.nelems(ph0)
  ph1 = PHS[1]; nmv1 = path.nelems(ph1)
  ph2 = PHS[2]; nmv2 = path.nelems(ph2)
    
  bcA = BCS[0]
  bcB = BCS[1]

  szmin = 0.9
  rszmin = 0.19
  tol = 0.20*wdf  # Tolerance for overlaps, tilts, etc.
  ydir = (1,1)
  MVS0 = [ path.elem(ph0,kmv) for kmv in range(nmv0) ]
  MVS1 = [ path.elem(ph1,kmv) for kmv in range(nmv1) ]

  for ifun in range(6):
    if ifun == 0 or ifun == 1:
      if ifun == 0:
        sys.stderr.write("  ... {from_move_lists} ...\n")
        OMVS0 = [ mv for mv in MVS0 if not move.is_jump(mv) ]
        OMVS1 = [ mv for mv in MVS1 if not move.is_jump(mv) ]
      else:
        sys.stderr.write("  ... {from_move_lists} (reversed) ...\n")
        OMVS0 = [ move.rev(mv) for mv in MVS0 if not move.is_jump(mv) ]
        OMVS1 = [ move.rev(mv) for mv in MVS1 if not move.is_jump(mv) ]
      CTS = contact.from_move_lists(OMVS0, OMVS1, szmin, rszmin, tol, ydir)
      sys.stderr.write("OMVS0:\n")
      move.show_list(sys.stderr, " ", OMVS0, None)
      sys.stderr.write("OMVS1:\n")
      move.show_list(sys.stderr, " ", OMVS1, None)
    elif ifun == 2 or ifun == 3:
      if ifun == 2:
        sys.stderr.write("  ... {from_paths} ...\n")
        oph0 = ph0; oph1 = ph1
      else:
        sys.stderr.write("  ... {from_paths} (reversed) ...\n")
        oph0 = path.rev(ph0); oph1 = path.rev(ph1)
      CTS = contact.from_paths(oph0, oph1, szmin, rszmin, tol, ydir) 
      path.show(sys.stderr, " oph0 = ", oph0, "\n", True, 0,0,0)
      path.show(sys.stderr, " oph1 = ", oph1, "\n", True, 0,0,0)
    elif ifun == 4:
      sys.stderr.write("  ... {from_path_lists} ...\n")
      OPHS0 = [ph0,path.rev(ph2),]
      OPHS1 = [ph1,]
      CTS = contact.from_path_lists(OPHS0, OPHS1, szmin, rszmin, tol, ydir) 
      sys.stderr.write("OPHS0:\n")
      path.show_list(sys.stderr, " ", OPHS0, None, True,False)
      sys.stderr.write("OPHS1:\n")
      path.show_list(sys.stderr, " ", OPHS1, None, True,False)
    elif ifun == 5:
      sys.stderr.write("  ... {from_blocks} ...\n")
      block.show(sys.stderr, " bcA = ", bcA, "\n", True, 0,0)
      block.show(sys.stderr, " bcB = ", bcB, "\n", True, 0,0)
      CTS = contact.from_blocks(bcA, bcB, szmin, rszmin, tol, ydir) 
    else:
      assert False
    
    # Show and check:
    contact.show_list(sys.stderr, "    ", CTS, None)
    for kct in range(len(CTS)):
      ctk = CTS[kct]
      contact.show(sys.stderr, ("  contact %d: " % kct), ctk, "\n", 4)

      # Check tcool_limit:
      assert contact.tcool_limit(ctk) == +inf
      tclim_exp = 15.0
      contact.set_tcool_limit(ctk, tclim_exp)
      assert contact.tcool_limit(ctk) == tclim_exp
    
    # Check contacts by names of sides:
    CTNS_obs = [ ( move.get_name(contact.side_move(ct, 0)), move.get_name(contact.side_move(ct, 1)) ) for ct in CTS ]
    CTNS_obs = list.sort(CTNS_obs)

    CTNS_exp = [
      ("TGa0", "TGb0"),
      ("TGc0", "TGb0"),
      ("TGa1", "TGb0"),
      ("TGa1", "TGb1"),
      ("TGc2", "TGb0"),
      ("TGc2", "TGb1"),
    ]
    CTNS_exp = list.sort(CTNS_exp)
    
    assert CTNS_obs == CTNS_exp

  return 
  # ----------------------------------------------------------------------

def test_side_paths():
  
  sys.stderr.write("  ... testing {clear_side_paths,add_side_path,get_side_paths} ...\n")

  OPHS, TRS, JMS = path_example.misc_E(mp_fill, mp_jump)
  path.show_list(sys.stderr, "    ", OPHS, None, True, True)
  move.show_list(sys.stderr, "    ", TRS, None)
  
  Pd = OPHS[3]; assert path.get_name(Pd) == "Pd"
  Pe = OPHS[4]; assert path.get_name(Pe) == "Pe"
  
  szmin = 0.1
  rszmin = 0.02
  tol = 0.20*wdf  # Tolerance for overlaps, tilts, etc.
  CTS = contact.from_paths(Pd,Pe, szmin, rszmin, tol, None)
  assert CTS != None and len(CTS) == 2
  contact.show_list(sys.stderr, "    ", CTS, None)
  
  ct = CTS[0]
  
  # The {Move} objects that are sides of the contact:
  mv0 = TRS[7]; assert move.get_name(mv0) == "Td0"
  mv1 = TRS[9]; assert move.get_name(mv1) == "Te0"
  
  assert contact.side_move(ct, 0) == mv0
  assert contact.side_move(ct, 1) == mv1
  
  # Some paths that contain the sides of the contact:
  assert path.find_move(Pd, mv0) == 0
  assert path.find_move(Pe, mv1) == 0

  for round in range(2):
    
    for isd in 0, 1:
      for P in Pd, Pe:
        path.clear_contacts(P, isd)
        assert len(path.get_contacts(P, isd)) == 0
      contact.clear_side_paths(ct,isd)
      assert len(contact.get_side_paths(ct, isd)) == 0

    for isd in range(2):
      P0,dr0 = path.unpack((Pd,Pe)[isd])
      for dr in range(2):
        ophr = (P0,dr)
        imv = path.find_move(ophr,contact.side_move(ct,isd))
        assert imv != None
        contact.add_side_path(ct, isd, ophr, imv)
        assert (P0,dr,imv) in contact.get_side_paths(ct, isd)
      path.add_contact(P0, isd, ct)
      assert ct in path.get_contacts(P0, isd)
      assert ct in path.get_contacts(path.rev(P0), isd)
  
  return
  # ----------------------------------------------------------------------


def test_side_index():
  
  sys.stderr.write("  ... testing {set_side_index,get_side_index} ...\n")

  pa = (1, 1 + 0*wdf)
  qa = (3, 1 + 0*wdf)
  tra = move.make(pa, qa, mp_fill)

  pb = (2, 1 + 1*wdf)
  qb = (5, 1 + 1*wdf)
  trb = move.make(pb, qb, mp_fill)
  
  ct = contact.from_moves(tra, trb, None, 0.1, 0.0, 0.1)
  assert ct != None
  
  for kix in None, 10, 20:
    for isd in 0, 1: 
      ix = None if kix == None else kix + 7*isd
      contact.set_side_index(ct, isd, ix)
    for isd in 0, 1: 
      ix_exp = None if kix == None else kix + 7*isd
      ix_cmp = contact.get_side_index(ct, isd)
      assert ix_cmp == ix_exp
  
  return
  # ----------------------------------------------------------------------

def test_show():

  sys.stderr.write("--- testing {show,show_list} ---\n")

  CTS, OPHS, TRS = contact_example.misc_B(mp_fill, mp_jump)
  
  sys.stderr.write("  ... {show} ...\n")
  wna = 5
  sys.stderr.write("\n")
  for kct in range(len(CTS)):
    ct = CTS[kct]
    contact.show(sys.stderr, "    [", ct, "]\n", wna)
    wna = wna + 2
  sys.stderr.write("\n")

  sys.stderr.write("  ... {show_list} ...\n")
  contact.show_list(sys.stderr, "    < ", CTS, " >")
  return 
  # ----------------------------------------------------------------------

def test_plot_to_files():

  sys.stderr.write("--- testing {plot_to_files} ---\n")

  tag = "plot_to_files"
  CTS, OPHS, TRS = contact_example.misc_B(mp_fill, mp_jump)
  nph = len(OPHS)

  CLRS = [ pyx.color.rgb(0.300, 0.600, 0.000), ]
  nclr = len(CLRS)

  rwd = 0.80
  wd_axes = 0.05*wdf   
  clr_ct = pyx.color.rgb.red # Color for contact lines

  for do_tics in (False, True):
    for do_arrows in (False, True):
      dashops = (False,) if do_tics or do_arrows else (False, True,)
      for do_dashed in dashops:
        if not (do_tics and do_arrows):
          if do_dashed:
            dashpat = (2*wd_axes, 1.5*wd_axes)
            ext = 0.35*wdf
          else:
            dashpat = None
            ext = 0
          subtag = "ds%s_tc%d_ar%d" % (int(do_dashed),int(do_tics),int(do_arrows))
          fname = ("tests/out/contact_TST_%s_%s" % (tag,subtag))
          contact.plot_to_files(fname, CTS, clr_ct, dashpat, ext, OPHS, CLRS, rwd, wd_axes, do_tics, do_arrows)
  return
  # ----------------------------------------------------------------------

def test_plot_single():

  sys.stderr.write("--- testing {plot_single} ---\n")

  tag = "plot_single"
  CTS, OPHS, TRS = contact_example.misc_B(mp_fill, mp_jump)
  nph = len(OPHS)

  CLRS = hacks.trace_colors(nph, None)
  nclr = len(CLRS)

  # Get the enclosing box of the paths:
  B = path.bbox(OPHS)
  B = rn.box_join(B, contact.bbox(CTS))
  
  dp = (0,0)
  
  wd_axes = 0.05*wdf 
  rwd = 0.80
  wd_ct = 1.5*wd_axes
  clr_ct = pyx.color.rgb.red # Color for contact lines

  for do_tics in (False, True):
    for do_arrows in (False, True):
      dashops = (False,) if do_tics or do_arrows else (False, True,)
      for do_dashed in dashops:
        if not (do_tics and do_arrows):
          c, szx,szy = hacks.make_canvas(hacks.round_box(B,0.5), dp, None, True, True, 1, 1)

          axes = False
          dots = True
          ph_arrows = True
          matter = False
          path.plot_standard(c, OPHS, None, None, CLRS, rwd, wd_axes, axes, dots, ph_arrows, matter)

          sz_tics = wd_ct if do_tics else 0
          nct = len(CTS)
          if do_dashed:
            dashpat = (2*wd_ct, 1.5*wd_ct)
            ext = 0.35*wdf
          else:
            dashpat = None
            ext = 0
          for kct in range(nct):
            contact.plot_single(c, CTS[kct], None, clr_ct, dashpat, ext, wd=wd_ct, sz_tic=sz_tics, arrow=do_arrows)

          subtag = "ds%s_tc%d_ar%d" % (int(do_dashed),int(do_tics),int(do_arrows))
          fname = ("tests/out/contact_TST_%s_%s" % (tag,subtag))
          hacks.write_plot(c, fname)

  return
  # ----------------------------------------------------------------------

def test_endpoints_sides():
  sys.stderr.write("--- testing {endpoints,pmid,side,which_side,ixcovs} ---\n")
  
  CTS, OPHS, TRS = contact_example.misc_B(mp_fill, mp_jump)
  
  # Testing {pmid}
  q00, q01 = contact.endpoints(CTS[0])
  m1a = rn.mix(0.50, q00, 0.50, q01)
  m1b = contact.pmid(CTS[0])
  assert rn.dist(m1a, m1b) < 1.0e-8
 
  # This validation depends on the specific {OPHS,CTS} created above:
  for   ct,    mv0,   mv1,   oph0,   ix0, oph1,   ix1 in ( 
      ( CTS[0], TRS[0], TRS[1], OPHS[0], 0,   OPHS[2], 2 ),
      ( CTS[1], TRS[1], TRS[2], OPHS[2], 2,   OPHS[1], 0 ),
      ( CTS[2], TRS[0], TRS[2], OPHS[0], 0,   OPHS[1], 0 ),
      ( CTS[3], TRS[2], TRS[3], OPHS[1], 0,   OPHS[2], 0 ),
      ( CTS[4], TRS[1], TRS[4], OPHS[2], 2,   OPHS[1], 4 ),
      ( CTS[5], TRS[7], TRS[9], OPHS[3], 0,   OPHS[4], 0 ),
      ( CTS[6], TRS[8], TRS[9], OPHS[3], 2,   OPHS[4], 0 ),
      ( CTS[7], TRS[7], TRS[5], OPHS[3], 0,   OPHS[1], 2 ),
  ):
    contact.show(sys.stderr, "  contact = ", ct, "\n", 0)

    assert contact.side_move(ct, 0) == mv0
    assert contact.side_move(ct, 1) == mv1
    assert contact.which_side(mv0, ct) == 0
    assert contact.which_side(move.rev(mv0), ct) == 0
    assert contact.which_side(mv1, ct) == 1
    assert contact.which_side(move.rev(mv1), ct) == 1
    assert contact.which_side(TRS[6], ct) == None

    omv0 = path.elem(oph0, ix0)
    omv1 = path.elem(oph1, ix1)
    mv0a, drm0 = move.unpack(omv0)
    mv1a, drm1 = move.unpack(omv1)
    assert mv0a == mv0
    assert mv1a == mv1

    assert contact.path_ixcovs(oph0, ct) == ( ix0, None )
    assert contact.path_ixcovs(oph1, ct) == ( None, ix1 )

    ph0, drp0 = path.unpack(oph0)
    ph1, drp1 = path.unpack(oph1)
    if ph0 != ph1:
      nmv0 = path.nelems(oph0)
      nmv1 = path.nelems(oph1)
      gap = (path.pfin(oph0) != path.pini(oph1))
      use_links = False
      oph01 = path.concat((oph0, oph1), use_links, mp_jump)
      assert path.nelems(oph01) == nmv0 + int(gap) + nmv1
      assert contact.path_ixcovs(oph01, ct) == ( ix0, nmv0+int(gap)+ix1 )

  return
  # ----------------------------------------------------------------------

def do_test_side_tcov(ct):
  # Tests {side_tcov} on {t}:
  sys.stderr.write("... testing{side_tcov} ...\n")
  sys.stderr.write("ct = %s\n" % contact.get_name(ct))

  # Test {side_tcov}:
  for isd in 0, 1:
    mvi = contact.side_move(ct, isd)
    assert isinstance(mvi, move.Move)
    tcvi_cmp = contact.side_tcov(ct, isd)
    tcvi_exp = move.cover_time(mvi, contact.pmid(ct))
    assert abs(tcvi_exp - tcvi_cmp) < 1.0e-8
  return
  # ----------------------------------------------------------------------

def do_test_path_tcovs(oph, ct):
  # Tests {path_tcov,path_tcovs,path_ixcovs} on {oph,ct}:

  sys.stderr.write("... testing{path_tcov,path_tcovs,path_ixcovs} ...\n")
  sys.stderr.write("ct = %s oph = %s\n" % (contact.get_name(ct),path.get_name(oph)))

  # Test {path_tcovs}:
  ixs_cmp = contact.path_ixcovs(oph, ct)
  assert type(ixs_cmp) is list or type(ixs_cmp) is tuple
  assert len(ixs_cmp) == 2

  tcs_cmp = contact.path_tcovs(oph, ct)
  assert type(tcs_cmp) is list or type(tcs_cmp) is tuple
  assert len(tcs_cmp) == 2

  for isd in 0, 1:
    # Check side {isd} of contact
    mvi_cmp = contact.side_move(ct, isd)
    assert isinstance(mvi, move.Move)
    imv = path.find_move(oph, mvi)
    assert ixs_cmp[isd] == imv
    if imv != None:
      assert contact.path_tcov(oph, imv, ct, isd) == tcs_cmp[isd]
      # Check if indeed covers:
      omvi_exp = path.elem(oph, imv)
      mvi_exp, dri_exp = move_unpack(omvi_exp)
      assert mvi_cmp == mvi_exp
      # Check coverage time:
      tcovi_exp = path.tini(oph, imv) + move.cover_time(omvi_exp, contact.pmid(ct))
      assert abs(tcovi_exp - tcs_cmp[isd]) < 1.0e-8
    else:
      # Check if indeed does not cover:
      for jmv in range(path.nelems(oph)):
        omvj = path.elem(oph, jmv)
        mvj, drj = move_unpack(omvk,c)
        assert mvj != mvi_cmp
      assert tcs_cmp[isd] == None
  return
  # ----------------------------------------------------------------------

def do_test_est_rcool(oph0, ct):
  # Tests {est_rcool} on {oph0,ct}:
  # Assumes the paths and contacts of {contact_example.misc_K}.
  
  sys.stderr.write("... testing{tcool,rcool_max,est_rcool} ...\n")
  sys.stderr.write("oph0 = %s ct = %s\n" % (path.get_name(oph0),contact.get_name(ct)))

  sys.stderr.write("  ... testing {path_tcovs} ...\n")
  ixs0 = contact.path_ixcovs(oph0, ct)
  tcs0 = contact.path_tcovs(oph0, ct)
  
  sys.stderr.write("  ... testing {tcool} ...\n") 
  tcool_cmp = contact.tcool(oph0, ct)
  if tcs0[0] == None or tcs0[1] == None:
    tcool_exp = None
    assert tcool_cmp == None
  else:
    tcool_exp = abs(tcs0[0] - tcs0[1])
    assert abs(tcool_exp - tcool_cmp) < 1.0e-8

  # Compute {est_rcool_cmp} by the module:
  est_rcool_cmp = contact.est_rcool(oph0, ct, mp_jump, False)
  
  # Check the {quick} version:
  est_rcool_cmp_quick = contact.est_rcool(oph0, ct, mp_jump, True)
  if est_rcool_cmp_quick != est_rcool_cmp:
    assert est_rcool_cmp > 1.0
    assert est_rcool_cmp_quick == +inf
  
  # Compute {est_rcool_exp} by hand:
  if ixs0[0] == None and ixs0[1] == None: 
    # Path {oph0} does not cover {ct}:
    assert tcs0[0] == None and tcs0[1] == None
    est_rcool_exp = 0
  elif ixs0[0] != None and ixs0[1] != None:
    # Path {oph0} itself closes {ct}:
    est_rcool_exp = tcool_exp/contact.tcool_limit(ct)
  else:
    # Path {oph0} covers only one side of {ct}, find which:
    isd = 0 if ixs0[0] != None else 1
    # Compute {est_rcool_exp} as the best {rcool} obtainable with the other paths that cover {ct}
    est_rcool_exp = 0
    for ph1, dr1, imv1 in contact.get_side_paths(ct, 1-isd):
      oph1 = path.spin(ph1,dr1)
      # Now {oph1} (in that direction) is one of the paths that could cover the other side of {ct}.
      for use_links in False, True:
        oph = path.concat([oph0, oph1], use_links, mp_jump)
        ixs = contact.path_ixcovs(oph, ct)
        assert ixs[0] != None and ixs[1] != None
        sys.stderr.write("ixs0 = (%s,%s) ixs = (%d, %d)\n" % (str(ixs0[0]),str(ixs0[1]),ixs[0],ixs[1]))
        tcs = contact.path_tcovs(oph, ct)
        sys.stderr.write("tcs0 = (%s,%s) tcs = (%.8f, %.8f)\n" % (str(tcs0[0]),str(tcs0[1]),tcs[0],tcs[1]))
        assert ixs[isd] == ixs0[isd]
        assert abs(tcs[isd] - tcs0[isd]) < 1.0e-8
        rcool_u1 = abs(tcs[0] - tcs[1])/contact.tcool_limit(ct)
        est_rcool_exp = min(est_rcool_exp, rcool_u1)
  # Now compare hand and module results:
  assert est_rcool_cmp >= est_rcool_exp - 1.0e-8
  if est_rcool_cmp > est_rcool_exp + 1.0e-8:
    sys.stderr.write("!! {contact.est_rcool} too pessimistic?\n")
    sys.stderr.write("  est_rcool_cmp = %14.8f\n" % est_rcool_cmp)
    sys.stderr.write("  est_rcool_exp = %14.8f\n" % est_rcool_exp)
    assert False, "aborted"
    
  return
  # ----------------------------------------------------------------------
    
def do_test_est_max_rcool(oph, CTS):
  # Tests {est_max_rcool} on {oph,CTS}:
  # Assumes the paths and contacts of {contact_example.misc_K}.
  
  sys.stderr.write("... testing{est_max_rcool} ...\n")
  sys.stderr.write("oph0 = %s CTS = %s \n" % (path.get_name(oph), str([contact.get_name(ct) for ct in CTS])))
    
  # Discard contacts not covered by {oph}:
  CTS = [ ct for ct in CTS if contact.path_ixcovs(oph, ct) != (None,None) ]
  
  # Compute {est_max_rcool_cmp} by the module:
  est_max_rcool_cmp = contact.est_max_rcool(oph, CTS, mp_jump, False);
  
  # Check the quick version:
  est_max_rcool_cmp_quick = contact.est_max_rcool(oph, CTS, mp_jump, True);
  if est_max_rcool_cmp_quick != est_max_rcool_cmp:
    assert est_max_rcool_cmp > 1.0
    assert est_max_rcool_cmp_quick == +inf
    
  # Compute {est_max_rcool_exp} by hand:
  est_max_rcool_exp = 0
  for ct in CTS:
    est_rcool_k = contact.est_rcool(oph, ct, mp_jump, False);
    est_max_rcool_exp = max(est_max_rcool_exp, est_rcool_k)
  
  # Check results:
  assert est_max_rcool_exp == est_max_rcool_cmp
  return
  # ----------------------------------------------------------------------

def test_tcov_rcool():
  sys.stderr.write("--- testing {est_rcool,est_max_rcool} ---\n")

  CTS, OPHS = contact_example.misc_K(mp_fill, mp_jump)
  
  for ct in CTS:
    do_test_side_tcov(ct)
  
  for kph0 in range(len(OPHS)):
    for dr0 in 0, 1:
      oph = path.spin(OPHS[kph0], dr0)
      # Test cover and closing times with {oph} as the first path.
      for ct in CTS:
        sys.stderr.write("ct = %s\n" % (contact.get_name(ct)))
        do_test_est_rcool(oph, ct)
      sys.stderr.write("CTS = %s\n" % (str([contact.get_name(ct) for ct in CTS])))
      do_test_est_max_rcool(oph, CTS)

  return
  # ----------------------------------------------------------------------

def test_check_side_paths():
  sys.stderr.write("--- testing {check_side_paths} ---\n")
  cols = 5
  rows = 2
  nx = 3
  ny = 20
  islands = True
  mp_link = mp_cont
  OCRS, OPHS, OLKS, CTS = raster_example.patch_array(cols, rows, nx,ny, islands, mp_cont, mp_fill, mp_link)
  # OPHSr = [ path.rev(oph) for oph in OPHS ]
  contact.check_side_paths(CTS,OPHS,(0,1))
  return
  # ----------------------------------------------------------------------

test_show()
test_basic()
test_more_makes()
test_endpoints_sides()
test_side_index()
test_side_paths()
test_check_side_paths()
test_tcov_rcool()
test_plot_single()
test_plot_to_files()
