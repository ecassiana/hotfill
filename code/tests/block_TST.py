#! /usr/bin/python3
# Test program for module {block}

import block
import block_example
import move 
import move_parms
import path
import hacks
import job_parms
import rn
import pyx
import sys
from math import sqrt, sin, cos, floor, ceil, inf, nan, pi

parms = job_parms.typical_js()
parms['solid_raster_width'] = 1.00
parms['contour_trace_width'] = 0.50

mp_jump = move_parms.make_for_jumps(parms)
mp_cont = move_parms.make_for_contours(parms)
mp_fill = move_parms.make_for_fillings(parms)

wdf = move_parms.width(mp_fill)

def test_basics():
  sys.stderr.write("--- testing {from_paths,nchoices,choice,has_name,get_name,set_name,tag_names} ---\n")
  
  BCS, OPHS, TRS, JMS = block_example.misc_A(mp_fill, mp_jump)

  bc0 = BCS[0]; block.validate(bc0)
  bc1 = BCS[1]; block.validate(bc1)
  bc2 = BCS[2]; block.validate(bc2)
  
  mex = min(path.fabtime(OPHS[0]), min(path.fabtime(OPHS[1]), path.fabtime(OPHS[2])))
  assert block.min_fabtime(bc0) == mex

  assert block.nchoices(bc0) == 4
  assert block.nchoices(bc1) == 2
  assert block.nchoices(bc2) == 1
  avgn = block.avg_choices([bc0,bc1,bc2])
  assert abs(avgn - 2.0) <= 1.0e-12
  assert block.choice(bc0, 0) == OPHS[0]
  assert block.choice(bc0, 1) == path.rev(OPHS[0])
  assert block.choice(bc0, 2) == OPHS[1]
  assert block.choice(bc0, 3) == path.rev(OPHS[1])
  
  # Two blocks with no name:
  bce = block.from_paths(OPHS[0:2])
  bcf = block.from_paths(OPHS[1:2])
  
  assert not block.has_name(bce)
  assert block.get_name(bce) == "B?"
  block.set_name(bce, "Lego")
  assert block.has_name(bce)
  assert block.get_name(bce) == "Lego"
  
  sys.stderr.write("applying {block.tag_names}:\n")
  assert not block.has_name(bcf)
  block.tag_names([bce, bcf], "Tag.")
  assert block.get_name(bce) == "Tag.Lego"
  assert block.has_name(bcf)
  assert block.get_name(bcf) == "Tag.B?"

  return
  # ----------------------------------------------------------------------

def test_show():

  sys.stderr.write("--- testing {show,show_list} ---\n")

  BCS, OPHS, TRS, JMS = block_example.misc_A(mp_fill, mp_jump)
  
  sys.stderr.write("  ... {show} ...\n")
  wna = 5
  wnc = 5
  sys.stderr.write("\n")
  for kbc in range(len(BCS)):
    paths = (kbc % 2 == 1)
    bc = BCS[kbc]
    block.show(sys.stderr, "    [", bc, "]\n", paths, wna, wnc)
    wna = wna + 2
    wnc = max(wnc - 1, 0)
  sys.stderr.write("\n")
  
  sys.stderr.write("  ... {show_list} ...\n")
  block.show_list(sys.stderr, "  < ", BCS, " >", True)
  return 
  # ----------------------------------------------------------------------

def test_finding():
  sys.stderr.write("--- testing {moves,find_block_with_move,find_choice_with_move} ---\n")
  
  BCS, OPHS, TRS, JMS = block_example.misc_A(mp_fill, mp_jump)

  bc0 = BCS[0]
  bc1 = BCS[1]
  bc2 = BCS[2]
  
  MVS0_cp = block.moves(bc0)
  MVS0_ex = (TRS[0],JMS[0],TRS[6],TRS[2],JMS[1],TRS[5],JMS[2],TRS[4])
  assert set(MVS0_cp) == set(MVS0_ex)
  
  MVS1_cp = block.moves(bc1)
  MVS1_ex = (TRS[7],JMS[4],TRS[8],TRS[9])
  assert set(MVS1_cp) == set(MVS1_ex)
  
  MVS2_cp = block.moves(bc2)
  MVS2_ex = (TRS[3],JMS[3],TRS[1])
  assert set(MVS2_cp) == set(MVS2_ex)
  
  assert  block.find_block_with_move(BCS, TRS[5]) == 0
  assert  block.find_block_with_move(BCS, TRS[9]) == 1
  assert  block.find_block_with_move(BCS, JMS[3]) == 2
  
  assert  block.find_block_with_move(BCS, move.rev(TRS[5])) == 0
  assert  block.find_block_with_move(BCS, move.rev(TRS[9])) == 1
  assert  block.find_block_with_move(BCS, move.rev(JMS[3])) == 2

  assert  block.find_block_with_move(BCS[:2], JMS[3]) == None

  assert BCS[block.find_block_with_move(BCS, path.elem(OPHS[0],2))] == bc0
  assert BCS[block.find_block_with_move(BCS, path.elem(OPHS[1],2))] == bc0
  assert BCS[block.find_block_with_move(BCS, path.elem(OPHS[2],2))] == bc2

  assert BCS[block.find_block_with_move(BCS, JMS[0]) == bc1]
  assert BCS[block.find_block_with_move(BCS, TRS[7]) == bc2]
   
  assert block.find_choice_with_move(bc0, JMS[0]) == 0
  assert block.find_choice_with_move(bc0, TRS[6]) == 0
  assert block.find_choice_with_move(bc0, TRS[4]) == 2
  assert block.find_choice_with_move(bc0, JMS[2]) == 2
  assert block.find_choice_with_move(bc0, TRS[7]) == None

  assert block.find_choice_with_move(bc1, TRS[4]) == None
  assert block.find_choice_with_move(bc1, TRS[7]) == 0
  assert block.find_choice_with_move(bc1, JMS[4]) == 0
  assert block.find_choice_with_move(bc1, TRS[8]) == 0
  assert block.find_choice_with_move(bc1, TRS[9]) == 1
  
  assert block.find_choice_with_move(bc2, TRS[4]) == None
  assert block.find_choice_with_move(bc2, JMS[0]) == None
  assert block.find_choice_with_move(bc2, TRS[3]) == 0
  assert block.find_choice_with_move(bc2, TRS[1]) == 0
 
  return
  # ----------------------------------------------------------------------

def test_timing():
  sys.stderr.write("--- testing {min_fabtime} ---\n")
  
  BCS, OPHS, TRS, JMS = block_example.misc_A(mp_fill, mp_jump)
  
  txm_tot_exp = 0
  for kbc in range(len(BCS)):
    bc = BCS[kbc]
    txm_exp = +inf
    for ich in range(block.nchoices(bc)):
      oph = block.choice(bc, ich)
      txm_exp = min(txm_exp, path.fabtime(oph))
    txm_cmp = block.min_fabtime(bc)
    sys.stderr.write("    min_fabtime(bc):  expected %20.16f  returned %20.16f\n" % (txm_exp,txm_cmp))
    assert txm_cmp == txm_exp
    txm_tot_exp += txm_exp
  txm_tot_cmp = block.min_tot_fabtime(BCS)
  sys.stderr.write("tot_min_fabtime(bc):  expected %20.16f  returned %20.16f\n" % (txm_tot_exp,txm_tot_cmp))
  assert txm_tot_cmp == txm_tot_exp
  return
  # ----------------------------------------------------------------------

def test_plot_matter_shadow():
  sys.stderr.write("--- testing {plot_matter_shadow} ---\n")
  tag = "plot_matter_shadow"

  bc =  block_example.raster_rectangle((1,1), 5, 4, False, True, True, mp_fill, mp_jump)
  
  BCS = [bc,]
  
  dp = (2,3)
  B = block.bbox(BCS, True, True)

  c, szx,szy = hacks.make_canvas(hacks.round_box(B,0.5), dp, None, True, True, 1, 1)

  links = True
  block.plot_matter_shadow(c, BCS, dp, links)
  fname = ("tests/out/block_TST_%s" % tag)
  hacks.write_plot(c, fname)
  return
  # ----------------------------------------------------------------------

def test_bbox_barycenter():
  sys.stderr.write("--- testing {bbox_barycenter} ---\n")
  tag = "bbox_barycenter"
  
  mp_thin = move_parms.make(0.5, 100, 100, 0.0)
  mp_wide = move_parms.make(1.0, 100, 100, 0.0)
  
  PHS = [ None, None ]
  for iph in range(2):
    ydir = 1 - 2*iph
    pa = (1,4 + ydir*2)
    pb = (1,4 + ydir*0)
    pc = (5,4 + ydir*0)
    pd = (5,4 + ydir*2)
    mvab = move.make(pa, pb, mp_thin)
    mvbc = move.make(pb, pc, mp_wide)
    mvcd = move.make(pc, pd, mp_wide)
    mvda = move.make(pd, pa, mp_jump)
    ph = path.from_moves([mvab, mvbc, mvcd, mvda ])
    path.set_name(ph, "P%d" % iph, True)
    PHS[iph] = ph
  bc = block.from_paths(PHS)
  
  B = block.bbox([bc,], True, True)
  assert B == ((1, 2), (5,6))
  
  pbar = block.barycenter(bc)
  assert rn.box_include_point(B, pbar) == B
  assert abs(pbar[0] - 46/14) < 1.0e-8
  assert pbar[1] == 4.0

  return
  # ----------------------------------------------------------------------

def test_plot_standard():
  sys.stderr.write("--- testing {plot_standard} ---\n")
  tag = "plot_standard"

  BCS, CTS = block_example.misc_D(mp_fill, True)
  nbc = len(BCS)
  CLRS = hacks.trace_colors(nbc, None) # Colors to use for different blocks.
  nch = 0
  for bc in BCS: nch = max(nch, block.nchoices(bc))

  # Find the bounding box:
  B = block.bbox(BCS, True, True)

  dp = (5,5)
  wd_axes = 0.05*wdf
  rwd = 0.80

  c, szx,szy = hacks.make_canvas(hacks.round_box(B,0.5), dp, None, True, True, 1, nch)
  
  links = True
  contacts = True
  ystep = szy
  
  for ich in range(nch):
    dpi = rn.add(dp, (0, ich*ystep))
    block.plot_standard \
      ( c, BCS, ich, dpi, CLRS = CLRS, rwd=rwd, wd_axes=wd_axes,
        axes=True, dots=True, arrows=True, matter=True, links=links, contacts=contacts
      )
  fname = ("tests/out/block_TST_%s" % tag)
  hacks.write_plot(c, fname)
  return 
  # ----------------------------------------------------------------------

def test_plot_to_files():
  sys.stderr.write("--- testing {plot_to_files} ---\n")
  BCS, CTS = block_example.misc_D(mp_fill, True)
  tag = "plot_to_files"
  fname = ("tests/out/block_TST_%s" % tag)
  CLRS = hacks.trace_colors(len(BCS), None) # Colors to use for different blocks.
  rwd = 0.80
  wd_axes = 0.05*wdf
  links = True
  contacts = True
  block.plot_to_files(fname, BCS, CLRS, rwd, wd_axes, matter=True, links=links, contacts=contacts)
  return
  # ----------------------------------------------------------------------
  
test_bbox_barycenter()
test_basics()
test_show()
test_finding()
test_timing()
test_plot_matter_shadow()
test_plot_standard()
test_plot_to_files()
