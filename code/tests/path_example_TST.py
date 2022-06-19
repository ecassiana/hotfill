# /usr/bin/python3
# Test program for module {path_example}

import path_example
import move
import move_parms
import path
import job_parms
import hacks
import pyx
import rn
import sys
from math import sqrt, sin, cos, floor, ceil, inf, nan, pi

parms = job_parms.typical_js()
parms['solid_raster_width'] = 1.00
parms['contour_trace_width'] = 0.50

mp_jump = move_parms.make_for_jumps(parms)
mp_cont = move_parms.make_for_contours(parms)
mp_fill = move_parms.make_for_fillings(parms)

wd_fill = move_parms.width(mp_fill)
wd_cont = move_parms.width(mp_cont)

ac,sp,ud = move_parms.dynamics(mp_fill)
mp_link = move_parms.make(0.75*wd_fill, ac,sp,ud)

def show_and_plot(tag, OPHS, deco):
  # Plots the oriented paths in the list {OPHS} to files ending with
  # {tag}. If {deco} is true also prints the matter footprint and draws
  # axes, dots, and arrowheads on the traces.
  
  if OPHS == None: OPHS = []
  assert type(OPHS) is tuple or type(OPHS) is list
  assert len(OPHS) >= 1
  
  path.show_list(sys.stderr, "    ", OPHS, None, True, True)
  for oph in OPHS: path.validate(oph)
  
  dp = None
  B = path.bbox(OPHS)
  B = rn.box_expand(B, (4,4), (4,4))
  wd_axes = 0.15*wd_fill
  scale = None
  frame = False
  grid = True
  c, szx, szy = hacks.make_canvas(B, dp, scale, frame, grid, 1, 1)

  # Plot the links, if any:
  clr_link = pyx.color.rgb( 1.000, 0.000, 0.700 )
  rwd_link = 0.65
  axes_link = False
  dots_link = False
  arrows_link = False
  matter_link = True
  for oph in OPHS:
   for ophr in oph, path.rev(oph):
     OLKS = list(path.get_links(ophr))
     path.plot_standard(c, OLKS, dp, None, [clr_link,], rwd_link, wd_axes, axes_link, dots_link, arrows_link, matter_link)

  # Plot the paths:
  nph = len(OPHS)
  CLRS = hacks.trace_colors(nph, None)
  rwd_path = 0.80
  wd_axes = 0.15*min(wd_fill,wd_cont) # Width of jumps and axis lines.
  axes_path = deco
  dots_path = deco
  arrows_path = deco
  matter_path = True
  path.plot_standard(c, OPHS, dp, None, CLRS, rwd_path, wd_axes, axes_path, dots_path, arrows_path, matter_path)
  
  fname = "tests/out/path_example_TST_" + tag
  hacks.write_plot(c, fname)
  return
  # ----------------------------------------------------------------------

def show_and_plot_separately(tag, OPHS, deco):
  # Plots onto the {pyx} canvas {c} the oriented paths in the list {OPHS}
  # to files ending with {tag}.  Each path is plotted separately, displaced
  # vertically by a suitable amount.  The matter shadow of all paths is
  # ploted under every one. 
  #
  # If {deco} is true also prints the matter
  # footprint and draws axes, dots, and arrowheads on the traces.
  
  if OPHS == None: OPHS = []
  assert type(OPHS) is tuple or type(OPHS) is list
  nph = len(OPHS)
  assert nph >= 1
  
  path.show_list(sys.stderr, "    ", OPHS, None, True, True)
  for oph in OPHS: path.validate(oph)
  
  B = path.bbox(OPHS)
  
  dp = (0, 0)
  
  nph = len(OPHS)
  CLRS = hacks.trace_colors(nph, None)

  c, szx, szy = hacks.make_canvas(hacks.round_box(B,0.5), dp, None, True, True, 1, nph)
  ystep = szy

  for kph in range(nph):
    dpk = rn.add(dp, (0, kph*ystep))
  
    # Plot the matter shadows:
    rwd = 0.80
    wd_axes = 0.05*wd_fill
    axes = False
    dots = False
    arrows = False
    matter = True
    path.plot_standard(c, OPHS, dpk, 0, CLRS, rwd, wd_axes, axes, dots, arrows, matter)

    # Plot path {kph} without matter shadow:
    rwd = 0.80
    wd_axes = 0.05*wd_fill
    axes = deco
    dots = deco
    arrows = True
    matter = False
    path.plot_standard(c, [OPHS[kph],], dpk, None, [CLRS[kph],], rwd, wd_axes, axes, dots, arrows, matter)

  hacks.write_plot(c, "tests/out/path_example_TST_" + tag)
  return
  # ----------------------------------------------------------------------

def test_onion(Rc, Rf, nt):
  sys.stderr.write("--- testing {onion} Rc = %.3f Rf = %.3f nt = %d ---\n" % (Rc,Rf,nt))
  tag = "onion%05.2f_%05.2f_%03d" % (Rc,Rf,nt)
  Rest = max(Rc,Rc) + wd_fill/2
  Rbox = ceil(Rest + wd_fill) + 1
  ctr = (Rbox,Rbox)
  phase = 0.5*pi
  ph, Rin, Rex = path_example.onion(ctr, Rc, mp_cont, Rf, mp_fill, phase, nt, mp_jump)
  OPHS = [ph,]
  deco = True
  show_and_plot(tag, OPHS, deco)
  return
  # ----------------------------------------------------------------------

def test_gearloose():
  sys.stderr.write("--- testing {gearloose} ---\n")
  tag = "gearloose"
  R = 10
  zigzag = True
  ph = path_example.gearloose(R, zigzag, mp_cont, mp_fill, mp_jump)
  OPHS = [ph,]
  deco = False
  show_and_plot(tag, OPHS, deco)
  return
  # ----------------------------------------------------------------------

def test_raster_rectangle(axis, alt):
  sys.stderr.write("--- testing {raster_rectangle} axis = %d alt = %s ---\n" % (axis,str(alt)))
  plo = (1,1)
  nx = 4; ny = 5
  step = wd_fill
  n = ny if axis == 0 else nx
  sz = (nx-1)*step if axis == 0 else (ny-1)*step
  PHS,TRS,LJS0,LJS1 = path_example.raster_rectangle(plo, axis, n, alt, sz, step, mp_fill,mp_jump)
  path.show_list(sys.stderr, "    ", PHS, None, True, True)
  
  # Must plot each path in a separate file because they overlap:
  for kph in range(len(PHS)):
    ph = PHS[kph]
    tag = ("raster_rectangle_axis%d_alt%d_ph%d" % (axis,int(alt),kph))
    deco = True
    show_and_plot(tag, [ph,], deco)
  return
  # ----------------------------------------------------------------------

def test_spiral_rectangle():
  sys.stderr.write("--- testing {spiral_rectangle} ---\n")
  tag = "spiral_rectangle"
  plo = (1,1)    # Starting point of first spiral
  szx = 3.27; szy = 5.31
  PHS = []
  for axis in range(2):
    for dr in range(2):
      sys.stderr.write("  ... axis = %d dr = %d ...\n" % (axis,dr))
      pini = ( plo[0] + dr*(2*szx + 2) , plo[1] + axis*(2*szy + 2) )
      ph = path_example.spiral_rectangle(pini, (1-2*dr)*szx, (1-2*axis)*szy, axis, mp_fill)
      PHS.append(ph)
  deco = True
  show_and_plot(tag, PHS, deco)
  return
  # ----------------------------------------------------------------------

def test_misc_A():
  sys.stderr.write("--- testing {misc_A} ---\n")
  oph = path_example.misc_A(mp_fill, mp_link, mp_jump) 
  deco = True
  show_and_plot("misc_A", (oph,), deco)
  return
  # ----------------------------------------------------------------------

def test_misc_B():
  sys.stderr.write("--- testing {misc_B} ---\n")
  oph = path_example.misc_B(mp_fill, mp_jump) 
  deco = True
  show_and_plot("misc_B", (oph,), deco)
  return
  # ----------------------------------------------------------------------

def test_misc_C():
  sys.stderr.write("--- testing {misc_C} ---\n")
  OPHS = path_example.misc_C(mp_fill, mp_link, mp_jump)
  deco = True
  show_and_plot("misc_C", OPHS, deco)
  return
  # ----------------------------------------------------------------------

def test_misc_D():
  sys.stderr.write("--- testing {misc_D} ---\n")
  OPHS, TRS, JMS = path_example.misc_D(mp_fill, mp_jump)
  deco = True
  show_and_plot("misc_D", OPHS, deco)
  return
  # ----------------------------------------------------------------------

def test_misc_E():
  sys.stderr.write("--- testing {misc_E} ---\n")
  OPHS, TRS, JMS = path_example.misc_E(mp_fill, mp_jump)
  deco = True
  show_and_plot("misc_E", OPHS, deco)
  return
  # ----------------------------------------------------------------------

def test_misc_F(alt):
  sys.stderr.write("--- testing {misc_F} alt = %s ---\n" % str(alt))
  tag = ("misc_F_alt%s" % int(alt))

  ph = path_example.misc_F(alt, mp_fill, mp_jump)
  OPHS = [ ph, ]
  deco = True
  show_and_plot(tag, OPHS, deco)
  return
  # ----------------------------------------------------------------------

def test_misc_G():
  sys.stderr.write("--- testing {misc_G} ---\n")
  tag = ("misc_G")

  PHS, TRS02, TRS1 = path_example.misc_G(mp_cont, mp_fill, mp_jump)
  move.show_list(sys.stderr, "    ", TRS02, None)
  move.show_list(sys.stderr, "    ", TRS1, None)
  deco = True
  show_and_plot_separately(tag, PHS, deco)
  return
  # ----------------------------------------------------------------------

def test_misc_H():
  sys.stderr.write("--- testing {misc_H} ---\n")
  tag = ("misc_H")

  OPHS = path_example.misc_H(mp_fill)
  deco = True
  show_and_plot(tag, OPHS, deco)
  # ----------------------------------------------------------------------

def test_misc_J():
  sys.stderr.write("--- testing {misc_J} ---\n")
  tag = ("misc_J")

  OPHS = path_example.misc_J(mp_cont, mp_fill)
  deco = True
  show_and_plot(tag, OPHS, deco)
  # ----------------------------------------------------------------------

def test_misc_K():
  sys.stderr.write("--- testing {misc_K} ---\n")
  tag = ("misc_K")

  nmv = 10
  step = 2.0*wd_fill
  pert = 1.0*wd_fill
  oph = path_example.misc_K(nmv,step,pert,mp_cont,mp_fill,mp_link,mp_jump)
  deco = False
  show_and_plot(tag, [oph,], deco)
  # ----------------------------------------------------------------------

def test_contours_A():
  sys.stderr.write("--- testing {contours_A} ---\n")
  tag = ("contours_A")

  OCRS, PTSS = path_example.contours_A(mp_cont)
  deco = True
  
  # ??? Should check nesting informattion ???
  
  show_and_plot(tag, OCRS, deco)
  # ----------------------------------------------------------------------

def test_contours_B():
  sys.stderr.write("--- testing {contours_B} ---\n")
  tag = ("contours_B")

  OCRS = path_example.contours_B(mp_cont)
  
  # ??? Should check nesting informattion ???

  deco = True
  show_and_plot(tag, OCRS, deco)
  # ----------------------------------------------------------------------

test_misc_K()

### test_misc_A()
### test_misc_B()
### test_misc_C()
### test_misc_D()
### test_misc_E()
### test_misc_F(False)
### test_misc_F(True)
### test_misc_G()
### test_misc_H()
### test_misc_J()
### 
### test_contours_A()
### test_contours_B()
### 
### for axis in 0, 1:
###   for alt in False, True:
###     test_raster_rectangle(axis,alt)
### 
### test_spiral_rectangle()
### test_gearloose()
### test_onion(5,0.3,36)
### test_onion(5,0.3,35)
### test_onion(5,7,8)
### test_gearloose()
