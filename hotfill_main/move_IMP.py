# Implementation of module {move}
# Last edited on 2021-11-09 16:51:31 by stolfi

import move
import move_parms
import hacks
import rn
import pyx 
from math import nan, inf, sqrt
import sys

class Move_IMP:
  # The field {mv.endpts} is a 2-tuple with the two endpoints of the 
  # axis, in arbitrary order.   An oriented move {(mv,dr)} means 
  # that the motion is from {mv.endpts[dr]} to {mv.endpts[1-dr]}.
  #
  # The field {mv.mp} is a {Move_Parms} object that specifies the width 
  # and timing parameters of the move.
  # The field {mv.fabtime} is the time (in seconds) needed to fabricate it, NOT including the 
  # nozzle up/down or filament suck/refeed if it is a jump.

  def __init__(self, p0, p1, mp, tex):
    self.endpts = (p0, p1)
    self.mp = mp
    self.fabtime = tex
    self.name = None

def make(p0, p1, mp):
  assert hacks.is_point(p0); p0 = tuple(p0) # Make sure it is immutable.
  assert hacks.is_point(p1); p1 = tuple(p1) # Make sure it is immutable.
  assert mp != None and isinstance(mp, move_parms.Move_Parms)
  dpq = rn.dist(p0, p1)
  tex = move_parms.nozzle_travel_time(dpq, None, mp)
  return move.Move(p0, p1, mp, tex)

def parameters(omv):
  mv, dr = unpack(omv)
  return mv.mp

def is_jump(omv):
  if omv == None: return False
  mv, dr = unpack(omv)
  return move_parms.width(mv.mp) == 0

def is_trace(omv):
  if omv == None: return False
  mv, dr = unpack(omv)
  return move_parms.width(mv.mp) > 0

def width(omv):
  mv, dr = unpack(omv)
  return move_parms.width(mv.mp)

def pini(omv):
  mv, dr = unpack(omv)
  return mv.endpts[dr]

def pfin(omv):
  mv, dr = unpack(omv)
  return mv.endpts[1-dr]
  # ----------------------------------------------------------------------

def endpoints(omv):
  mv, dr = unpack(omv)
  pts = (mv.endpts[dr], mv.endpts[1-dr])
  return pts
  # ----------------------------------------------------------------------

def length(omv):
  mv, dr = unpack(omv)
  L = rn.dist(mv.endpts[0], mv.endpts[1])
  return L
  # ----------------------------------------------------------------------

def bbox(OMVS):
  B = None
  for omv in OMVS:
    B = rn.box_include_point(B, pini(omv))
    B = rn.box_include_point(B, pfin(omv))
  return B
  # ----------------------------------------------------------------------

def rev(omv):  
  mv, dr = unpack(omv)
  return (mv, 1-dr)

def spin(omv, dr):
  mv1, dr1 = unpack(omv)
  return (mv1, (dr1 + dr) % 2)

def unpack(omv):
  if isinstance(omv, move.Move):
    return omv, 0
  else:
    # sys.stderr.write("omv =%s\n" % str(omv))
    assert type(omv) is tuple
    assert len(omv) == 2
    mv, dr = omv
    assert isinstance(mv, move.Move)
    assert dr == 0 or dr == 1
    return mv, dr

def displace(omv, ang, disp, mp):
  p = tuple(rn.add(rn.rotate2(pini(omv), ang), disp))
  q = tuple(rn.add(rn.rotate2(pfin(omv), ang), disp))
  return make(p, q, mp)
  
def sort_by_midpoint(omv0, omv1, ydir):
  pm0 = rn.mix(0.5, pini(omv0), 0.5, pfin(omv0))
  pm1 = rn.mix(0.5, pini(omv1), 0.5, pfin(omv1))
  if ydir == None:
    y0 = pm0[1]; y1 = pm1[1]
  else:
    y0 = rn.dot(pm0, ydir)
    y1 = rn.dot(pm1, ydir)
  if y0 < y1: 
    return omv0, omv1
  else:
    return omv1, omv0
  # ----------------------------------------------------------------------

def shared_border(omv0, omv1, mvdir, tol):
  debug = False
  
  if is_jump(omv0) or is_jump(omv1): 
    if debug: sys.stderr.write("one or both moves are jumps\n")
    p0 = None; p1 = None
  else:
    mv0, dr0 = unpack(omv0)
    e0a, e0z = mv0.endpts # Endpoints of axis segment {S0} of {mv0}.
    wd0 = width(mv0)
    
    mv1, dr1 = unpack(omv1)
    e1a, e1z = mv1.endpts # Endpoints of axis segment {S1} of {mv1}.
    wd1 = width(mv1)

    ysep = (wd0 + wd1)/2 # Ideal separation between the two moves.
    
    if mvdir == None:
      # Try to guess a {mvdir} from the directions of the two moves:
      mvdir0, len0 = rn.dir(rn.sub(pfin(omv0), pini(omv0)))
      mvdir1, len1 = rn.dir(rn.sub(pfin(omv1), pini(omv1)))
      # If the moves have opposite orientation, reverse one of {mvdir0,mvdir1}:
      if rn.dot(mvdir0,mvdir1) < 0: mvdir1 = rn.scale(-1, mvdir1)
      if len0 >= tol and len1 >= tol and rn.dot(mvdir0,mvdir1) > 0.90:
        mvdirm, lenm = rn.dir(rn.add(mvdir0, mvdir1))
        assert lenm > 1.0
        mvdir = mvdirm
    
    if mvdir == None:
      # Guess failed, give up:
      p0 = None; p1 = None
    else:
      xdir = mvdir
      ydir = (-xdir[1], +xdir[0])

      # Compute projection of endpoints in direction {xdir}:
      x0a = rn.dot(e0a, xdir); x0z = rn.dot(e0z, xdir)
      x1a = rn.dot(e1a, xdir); x1z = rn.dot(e1z, xdir)

      # Ensure that projections are oriented along {xdir}:
      xlo0 = min(x0a, x0z); xhi0 = max(x0a, x0z)
      xlo1 = min(x1a, x1z); xhi1 = max(x1a, x1z)

      # Compute ovrlap in direction {xdir}:
      xlo = max(xlo0,xlo1); xhi = min(xhi0,xhi1)
      dx = xhi - xlo;  # length of overlap, or negative if no overlap.
      if dx <= 1.0e-6*ysep:
        if debug: 
          sys.stderr.write("move axes do not overlap along {xdir}")
          sys.stderr.write(" mv0 = [ %+10.6f _ %+10.6f ]" % (xlo0,xhi0))
          sys.stderr.write(" mv1 = [ %+10.6f _ %+10.6f ]" % (xlo1,xhi1))
          sys.stderr.write(" ovr = [ %+10.6f _ %+10.6f ]\n" % (xlo,xhi))
        p0 = None; p1 = None
      else:

        # Lift the X-range of overlap onto the segments:
        r0a = (xlo - x0a)/(x0z - x0a)
        c0a = rn.mix(1-r0a, e0a, r0a, e0z)
        r0z = (xhi - x0a)/(x0z - x0a)
        c0z = rn.mix(1-r0z, e0a, r0z, e0z)

        r1a = (xlo - x1a)/(x1z - x1a)
        c1a = rn.mix(1-r1a, e1a, r1a, e1z)
        r1z = (xhi - x1a)/(x1z - x1a)
        c1z = rn.mix(1-r1z, e1a, r1z, e1z)

        # Project the lifted points in the direction {ydir}:
        y0a = rn.dot(c0a, ydir); y0z = rn.dot(c0z, ydir); y0m = (y0a + y0z)/2
        y1a = rn.dot(c1a, ydir); y1z = rn.dot(c1z, ydir); y1m = (y1a + y1z)/2

        # Check for near-parallelism and proper distance perpendicularly to {xdir}:
        dy01 = abs(y0m - y1m) # Distance in {ydir} between midpoints.
        if abs(y0a - y0z) > tol or abs(y1a - y1z) > tol:
          if debug: sys.stderr.write("dy0 = %+.8f dy1 = %+.8f tol = %.8f\n" % (y0z-y0a,y1z-y1a,tol))
          if debug: sys.stderr.write("one of both segments are not parallel to {xdir}\n")
          p0 = None; p1 = None
        elif dy01 > ysep + tol:
          if debug: sys.stderr.write("the segments are too far apart in {ydir}\n")
          p0 = None; p1 = None
        elif dy01 < ysep - tol:
          if debug: sys.stderr.write("the segments overlap too much in {ydir}\n")
          p0 = None; p1 = None
        else:    
          ym = (wd1*y0m + wd0*y1m)/(wd0 + wd1)
          xm = (xlo + xhi)/2
          xr = (xhi - xlo)/2
          p0 = tuple(rn.mix(xm - xr, xdir, ym, ydir))
          p1 = tuple(rn.mix(xm + xr, xdir, ym, ydir))
  return (p0,p1)
  # ----------------------------------------------------------------------

def fabtime(omv):
  mv, dr = unpack(omv)
  return mv.fabtime

def ud_penalty(omv):
  mv, dr = unpack(omv)
  ac, sp, ud = move_parms.dynamics(mv.mp)
  if not is_jump(omv): assert ud == 0
  return ud
 
def transition_penalty(omv0, omv1):
  if omv0 == None or omv1 == None: return 0
  assert pfin(omv0) == pini(omv1), "moves are not connected"
  if is_jump(omv0) and is_trace(omv1):
    return ud_penalty(omv0)
  elif is_trace(omv0) and is_jump(omv1):
    return ud_penalty(omv1)
  else:
    return 0

def cover_time(omv, m):
  mv, dr = unpack(omv) # To type-check.
  ac, sp, ud = move_parms.dynamics(mv.mp)
  # Compute cover time {tomv} from start of move {omv}:
  p, q = endpoints(omv)
  vpm = rn.sub(m, p)
  vpq = rn.sub(q, p)
  dpq = rn.norm(vpq)
  dpm = rn.dot(vpm,vpq)/dpq  # Distance from {p} to point nearest to {m}
  # sys.stderr.write("dpq = %12.8f dpm = %12.8f ratio = %12.8f\n" % (dpq, dpm, dpm/dpq))
  if dpm < 0: dpm = 0
  if dpm > dpq: dpm = dpq
  tc = move_parms.nozzle_travel_time(dpq, dpm, mv.mp)
  return tc

def plot_to_files(fname, OMVS, CLRS, rwd, wd_axes):  
  assert type(OMVS) is list or type(OMVS) is tuple
  
  # Compute the plot's bounding box:
  B = bbox(OMVS)

  c, szx,szy = hacks.make_canvas(hacks.round_box(B,0.5), None, None, True, True, 1, 1)

  axes = True
  dots = True
  arrows = True
  matter = True

  plot_standard(c, OMVS, None, None, CLRS, rwd, wd_axes, axes, dots, arrows, matter)
  hacks.write_plot(c, fname)
  return
  # ----------------------------------------------------------------------

def plot_standard(c, OMVS, dp, layer, CLRS, rwd, wd_axes, axes, dots, arrows, matter):
  
  assert type(OMVS) is list or type(OMVS) is tuple
  nmv = len(OMVS)
  assert nmv > 0
  
  assert wd_axes != None and wd_axes > 0
  
  if CLRS == None: 
    CLRS = [ pyx.color.rgb(0.050, 0.800, 0.000), ] # Default trace color.
  else:
    assert type(CLRS) is list or type(CLRS) is tuple
  nclr = len(CLRS)
  assert nclr == 1 or nclr == nmv
  
  def pick_colors(kmv):
    # Returns the colors for trace sausages and axes of move {OMVS[kmv]}.
    if nclr == 1:
      ctrace = CLRS[0]
    else:
      ctrace = CLRS[kmv]
    caxis =   pyx.color.rgb(0.6*ctrace.r, 0.6*ctrace.g, 0.6*ctrace.b) # Color of trace axis, dots, arrow.
    return ctrace, caxis

  # Colors:
  cmatter = hacks.matter_color(None) # Est. material footprint.
  cjumps = pyx.color.rgb.black                 # Axis lines, dots, and arrowheads of jumps.

  # Dimensions relative to trace nominal widths:
  rwd_matter = 1.13;  # Estimated material.
  
  # Absolute dimenstons (in mm):
  wd_dots = 2.5*wd_axes
  sz_arrows = 6*wd_axes

  # Get list {lys} of layers to plot:
  if layer == None:
    # Plot all four layers.
    lys = (range(4))
  else:
    # Plot only the selected layer.
    assert type (layer) is int
    assert layer >= 0 and layer < 4
    lys = (layer,)

  # Plot the layers:
  for ly in lys:
    # sys.stderr.write("{move.plot_standard}: layer %d\n" % ly)
    for kmv in range(nmv):
      omv = OMVS[kmv]
      mv, dr = unpack(omv) # For the typechecking.
      jmp = is_jump(omv)
      wd = width(omv)
      if ly == 0 and not jmp and matter:
        # plots the estimate of actual material:
        wd_matter = rwd_matter*wd
        plot_layer \
          (c, omv, dp, clr=cmatter, wd=wd_matter, dashed=False, wd_dots=0, sz_arrow=None)
      elif ly == 1 and not jmp:
        # plots the nominal trace material:
        wd_trace = rwd*wd
        ctrace, caxis = pick_colors(kmv)
        plot_layer \
          (c, omv, dp, clr=ctrace, wd=wd_trace, dashed=False, wd_dots=0, sz_arrow=None)
      elif ly == 2 and not jmp and (axes or dots or arrows):
        # Trace axis and/or dots and/or arrowhead:
        ctrace, caxis = pick_colors(kmv)
        t_wd_axis = wd_axes if axes else 0
        t_wd_dots = wd_dots if dots else 0
        t_sz_arrow = sz_arrows if arrows else 0
        plot_layer(c, omv, dp, clr=caxis, wd=t_wd_axis, dashed=False, wd_dots=t_wd_dots, sz_arrow=t_sz_arrow)
      elif ly == 3 and jmp:
        # Jump axis, dots, and arrowhead:
        j_wd_axis = wd_axes
        j_wd_dots = wd_dots
        j_sz_arrow = sz_arrows
        plot_layer(c, omv, dp, clr=cjumps, wd=j_wd_axis, dashed=True, wd_dots=j_wd_dots, sz_arrow=j_sz_arrow)
  return
  # ----------------------------------------------------------------------

def plot_layer(c, omv, dp, clr, wd, dashed, wd_dots, sz_arrow):
  
  if clr == None: return
  
  # Simplifications:
  if wd == None: wd = 0 
  if wd_dots == None: wd_dots = 0 
  if sz_arrow == None: sz_arrow = 0
  assert wd >= 0 and wd_dots >= 0 and sz_arrow >= 0 
  
  # Get the move endpoints:
  p, q = endpoints(omv)
  vpq = rn.sub(q,p)
  dpq = rn.norm(vpq)
  
  # Omit the arrowhead if not enough space for it:
  if dpq <= 2*sz_arrow: sz_arrow = 0
  
  # Nothing to plot if nothing is requested:
  if wd == 0 and wd_dots == 0 and sz_arrow == 0: return
  
  arrowpos = 0.5  # Position of arrow along axis.

  # Perturbation to force painting of zero-length lines.
  eps = 0.0001*max(wd,wd_dots,sz_arrow)  # For plotting the dots.
  assert eps > 0
  
  if wd == 0 and sz_arrow > 0:
    # We want an invisible line but still with the arrowhead.
    # Unfortunately setting linewidth(0) still draws a thin line.
    # So we cook things up, by moving {p} and {q} right _next 
    # to the arrowhead.
    m = rn.mix(1-arrowpos, p, arrowpos, q)  # Posititon of arrow.
    shaft = 0.9*sz_arrow  # Length of reduced arrow shaft.
    pa = arrowpos*shaft/dpq
    qa = (1-arrowpos)*shaft/dpq
    paxis = rn.mix(1, m, -pa, vpq)
    qaxis = rn.mix(1, m, +qa, vpq)
  elif dpq == 0:
    # Perturb {p,q} to ensure that the zero-length line is drawn as a dot:
    paxis = rn.sub(p,(eps,eps))
    qaxis = rn.add(q,(eps,eps))
  else:
    paxis = p; qaxis = q
  
  # Define styles {sty_axis} for the axis, {sty_dots} for the dots:
  sty_comm = [ pyx.style.linecap.round, pyx.style.linejoin.round, clr, ] # Common style.
  
  sty_axis = sty_comm + [ pyx.style.linewidth(wd) ]
  if dp != None: sty_axis.append(pyx.trafo.translate(dp[0], dp[1]))

  if sz_arrow > 0: 
    # Define the arrow style {sty_deco}:
    wdarrow = sz_arrow/8 # Linewidth for stroking the arrowhead (guess).
    sty_deco = sty_comm + [ pyx.style.linewidth(wdarrow) ]
    sty_deco = sty_deco + [ pyx.deco.stroked([pyx.style.linejoin.round]), pyx.deco.filled([]) ]
    
    # Add an arrowhead in style {sty_deco} to {sty_axis}:
    sty_axis = sty_axis + [ 
      pyx.deco.earrow(sty_deco, size=sz_arrow, constriction=None, pos=arrowpos, angle=35)
    ]
  
  if wd_dots > 0:
    # Plot dots:
    sty_dots = sty_comm + [ pyx.style.linewidth(wd_dots) ]
    if dp != None: sty_dots.append(pyx.trafo.translate(dp[0], dp[1]))
    c.stroke(pyx.path.line(p[0]+eps, p[1]+eps, p[0]-eps, p[1]-eps), sty_dots)
    c.stroke(pyx.path.line(q[0]+eps, q[1]+eps, q[0]-eps, q[1]-eps), sty_dots)

  if wd > 0 and dashed: 
    # Define dash pattern, or turn off dashing if too short:
    rdashpat = hacks.adjust_dash_pattern(dpq/wd, 1.75, 2.00)
    if rdashpat != None:
      # Make the line style dashed:
      sty_axis = sty_axis + [ pyx.style.linestyle(pyx.style.linecap.round, pyx.style.dash(rdashpat)) ]

  if wd > 0 or sz_arrow > 0:
    # Stroke the axis, with fixed endpoints:
    c.stroke(pyx.path.line(paxis[0], paxis[1], qaxis[0], qaxis[1]), sty_axis)

# DEBUGGING AND TESTING

def has_name(omv):
  mv, dr = unpack(omv)
  return mv.name != None
  # ----------------------------------------------------------------------

def get_name(omv):
  mv, dr = unpack(omv)
  name = mv.name
  if name == None: name = "J?" if is_jump(mv) else "T?"
  if dr == 1: name = "~" + name
  return name
  # ----------------------------------------------------------------------

def set_name(omv, name):
  assert type(name) is str
  mv, dr = unpack(omv)
  mv.name = name
  return
  # ----------------------------------------------------------------------

def tag_names(OMVS, tag):
  if tag != None and tag != "":
    assert type(tag) is str
    for omv in OMVS:
      mv, dr = unpack(omv)
      mv.name = tag + get_name(mv)
  return
  # ----------------------------------------------------------------------

def show(wr, pref, omv, suff, wna):
  mv, dr = unpack(omv)
  if pref != None: wr.write(pref)
  wr.write(" " if dr == 0 else "~")
  wr.write("%-*s" % (wna,get_name(mv)))
  p = pini(omv)
  q = pfin(omv)
  wr.write(" (%6.1f,%6.1f)" % (p[0],p[1]))
  wr.write(" (%6.1f,%6.1f)" % (q[0],q[1]))
  wd = width(mv)
  if wd == 0:
    wr.write(" jmp")
  else:
    wr.write(" %3.1f" % wd)
  if p == q:
    wr.write(" trivial")
  elif p[0] == q[0]:
    wr.write(" vertical")
  elif p[1] == q[1]:
    wr.write(" horizontal")
  if suff != None: wr.write(suff)
  return
  # ----------------------------------------------------------------------

def show_list(wr, pref, OMVS, suff):
  nmv = len(OMVS)
  if nmv == 0: return
  wix = len(str(nmv-1)) # Digits for index.
  wna = 3 # Width of "name" column; min 3 because of the header.
  for omv in OMVS:
    mv, dr = unpack(omv)
    wna = max(wna, len(get_name(mv)))
  
  wr.write("\n")

  # Write header:
  wr.write("%*s%*s %-*s %*s %*s  wd obs\n" % (len(pref),'',wix,"k",wna+1,"name",15,"pini",15,"pfin"))
  wr.write("%*s%s %s %s %s --- %s\n" % (len(pref),'',"-"*wix,"-"*(wna+1),"-"*15,"-"*15,"-"*15))
  
  # Write moves:
  for kmv in range(len(OMVS)):
    omv = OMVS[kmv]
    if pref != None: wr.write(pref)
    wr.write("%*d " % (wix, kmv))
    show(wr, None, omv, suff, wna)
    wr.write("\n")

  wr.write("\n")
  return
  # ----------------------------------------------------------------------
  
