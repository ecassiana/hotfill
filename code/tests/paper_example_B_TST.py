# /usr/bin/python3
# Test program for module {paper_example_B}

import paper_example_B
import move
import move_parms
import path
import contact
import input_data
import txt_write
import hacks
import pyx
import rn
import sys
from math import sqrt, sin, cos, floor, ceil, inf, nan, pi

# Some dynamics parameters:
ac = 3000
sp_cont = 20
sp_fill = 40
sp_link = sp_fill
sp_jump = 130
ud_jump = 0.05

# Move parameters:
wd_cont = 0.75; mp_cont = move_parms.make(wd_cont, ac, sp_cont, 0.00)
wd_fill = 1.00; mp_fill = move_parms.make(wd_fill, ac, sp_fill, 0.00)
wd_link = 0.50; mp_link = move_parms.make(wd_link, ac, sp_link, 0.00)
wd_jump = 0.00; mp_jump = move_parms.make(wd_jump, ac, sp_jump, ud_jump)

wd_axes = 0.15*min(wd_fill,wd_cont)

def test_make_simple_examples(which):

  sys.stderr.write("--- testing {make_simple_path,make_simple_moves,make_simple_contacts,make_simple_cover} which = %s ---\n" %  which)
  
  OCRS = []
  OLKS = []
  if which == "path":
    tag = "make_simple_path"
    OPHS = paper_example_B.make_simple_path(mp_fill, mp_jump)
    CTS = []
    plot_sep = False
  elif which == "moves":
    tag = "make_simple_moves"
    OPHS = paper_example_B.make_simple_moves(mp_fill, mp_jump)
    CTS = []
    plot_sep = True
  elif which == "contacts":
    tag = "make_simple_contacts"
    OPHS, CTS = paper_example_B.make_simple_contacts(mp_fill)
    plot_sep = False
  elif which == "cover":
    tag = "make_simple_cover"
    OPHS, CTS = paper_example_B.make_simple_cover(mp_fill,mp_jump)
    plot_sep = False
  else:
    assert False
    
  # Plot:
  o = None
  deco = True
  ydir = (0,1)
  if plot_sep:
    # Plot each path in {OPHS} separately:  
    for iph in range(len(OPHS)):
      oph = OPHS[iph]
      fname = "tests/out/paper_example_B_TST_" + tag + ("_%d" % iph)
      input_data.show_check_and_plot(fname, OCRS, [oph,], OLKS, CTS, o, ydir, wd_axes, deco)
  else:
    # Plot all paths together:
    fname = "tests/out/paper_example_B_TST_" + tag
    input_data.show_check_and_plot(fname, OCRS, OPHS, OLKS, CTS, o, ydir, wd_axes, deco)
  return
  # ----------------------------------------------------------------------

def test_make_turkey():

  sys.stderr.write("--- testing {make_turkey} ---\n")
  
  tag = "make_turkey"
  
  OCRS, OPHS, OLKS, CTS, VGS, EGS = paper_example_B.make_turkey(mp_cont, mp_fill, mp_link, mp_jump)
  o = None

  # Plot:
  deco = True
  ydir = (0,1)
  fname = "tests/out/paper_example_B_TST_" + tag
  sys.stderr.write("wd_axes = %.3f\n" % wd_axes)
  input_data.show_check_and_plot(fname, OCRS, OPHS, OLKS, CTS, o, ydir, wd_axes, deco)
  return
  # ----------------------------------------------------------------------

def test_write_txt():
  sys.stderr.write("--- testing txt output ---\n")
  OCRS,OPHS,OLKS,CTS,VGS,EGS = paper_example_B.make_turkey(mp_cont, mp_fill, mp_link, mp_jump) 

  fname = "tests/out/paper_example_B_TST.txt"
  wr = open(fname, "w")

  Z = 10.0
  angle = 0
  shift = (0,0)
  for oph in OPHS: path.set_group(oph,0)
  txt_write.write(wr, OCRS, OPHS, Z, angle, shift)
  wr.close()
  return
  # ----------------------------------------------------------------------

# Figures:

test_make_simple_examples("cover")
test_make_simple_examples("path")
test_make_simple_examples("moves")
test_make_simple_examples("contacts")
test_make_turkey()
test_write_txt()
