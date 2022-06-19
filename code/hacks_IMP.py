# Implementation of module {hacks}.

import color
import palette
import rn
import pyx
from math import sqrt, exp, log, floor, ceil, sin, cos, inf, nan, pi
from PIL import Image 
import sys
import os

def real_quadratic_roots(A, B, C):
  assert A != 0
  if A < 0:
    # Reverse signs to make {A} positive:
    A = -A; B = -B; C = -C  
  # Rescale by replacing {t = M s} so that {A} is 1:
  M = 1/sqrt(A)
  A = 1; B = B*M
  # sys.stderr.write("A = %10.7f B = %10.7f C = %10.7f M = %12.7f \n" % (A, B, C, M))
  # Now solve {s^2 + B*s + C = 0} 
  Delta = B*B - 4*C
  if Delta <= 0:
    # Zero or one root
    return None, None
  elif Delta == inf:
    # Pretend the roots are infinite but symmetric:
    return -inf, +inf
  else:
    sD = sqrt(Delta)
    s0 = (-B - sD)/2; t0 = M*s0
    s1 = (-B + sD)/2; t1 = M*s1
    assert -inf < t0 and t0 < t1 and t1 < +inf
    return t0, t1
  # ----------------------------------------------------------------------

# GEOMETRY

def is_point(p):
  if not (type(p) is tuple or type(p) is list): return False
  if len(p) != 2: return False
  for c in p:
    if (type(c) is not float) and (type(c) is not int): return False
  return True
  # ----------------------------------------------------------------------

def same_point(p, q, tol):
  return abs(p[0] - q[0]) <= tol and abs(p[1] - q[1]) <= tol
  # ----------------------------------------------------------------------

def poly_cleanup(PTS, dmin):
  if len(PTS) == 0: return PTS
  closed = same_point(PTS[0], PTS[-1], dmin) # Is polygon effectively closed?
  PTSnew = list.copy(PTS)
  n = 0 # Points {PTSnew[0..n-1]} are the non-deleted points.
  for p in PTS:
    while n > 0 and same_point(PTSnew[n-1], p, dmin):
      n = n - 1
    PTSnew[n] = p; n = n + 1
  # Check first and last points:
  while n >= 2 and same_point(PTSnew[0], PTSnew[n-1], dmin):
    n = n - 1
  PTSnew = PTSnew[0:n]
  if closed and PTSnew[0] != PTSnew[n-1]:  PTSnew.append(PTSnew[0])
  return PTSnew
  # ----------------------------------------------------------------------

def poly_orientation(PTS):
  A = rn.poly_area(PTS)
  if abs(A) < 1.0e-6: return 0
  return +1 if A > 0 else -1
  # ----------------------------------------------------------------------

def filet_center(p, q, rt, ro):
  rto = rt + ro
  upq, dpq = rn.dir(rn.sub(q, p))     # Unit vector directed from {p} to {q}
  vpq = (-upq[1], upq[0])  # Unit vector directed 90 degrees to the left of {p-->q}
  o = rn.mix(0.5, p, 0.5, q)
  h = sqrt(max(0, rto*rto - dpq*dpq/4)) 
  assert h > 0, "can't compute filet"
  t = rn.mix(1, o, -h, vpq)
  return t
  # ----------------------------------------------------------------------
      
def circle_x(y, yctr, rad):
  dy = yctr - y
  x = sqrt(max(0, rad*rad - dy*dy))
  return x
  # ----------------------------------------------------------------------

def cut_poly_with_band(CS, ydir, y0, y1):
  assert rn.dist(CS[0],CS[-1]) < 1.0e-8  # Contour component must be closed.
  npt = len(CS) # Number of points in contour component.    
  
  LS = [] 

  # Find first index {kpt} such that {CS[kpt]} is outside or on boundary of strip
  # and {CS[kpt+1]} is not on the same scan-line (indices modulo {npt}):
  k_start = None
  for kpt in range(npt):
    p = CS[kpt]
    q = CS[(kpt+1)%npt]
    yp = rn.dot(p, ydir)
    yq = rn.dot(q, ydir)
    if yp != yq:
      if yp <= y0 or yp > y1:
        k_start = kpt; break
  if k_start == None:
    # Component is entirely inside strip (?!).
    return LS
    
  PS = None # Points in current link path.
  eloc = None # Line through which the strip was entered (0, 1, or {None}).
  for dk in range(npt):
    kpt = k_start + dk
    p = CS[kpt % npt]
    q = CS[(kpt+1) % npt]
    ps, ploc, qs, qloc = clip_seg_to_strip(p, q, ydir, y0, y1)
    assert ploc == None or qloc == None or ploc != qloc
    if ps == None:
      assert qs == None
      assert PS == None # We should have closed it by hitting {y0} or {y1}.
    else:
      assert qs != None
      if eloc == None:
        # First seg of a new link:
        eloc = ploc
        PS = [ ps ]
      assert eloc != None
      assert PS != None
      if qs != ps:
        PS.append(qs)
        if qloc == 0 or qloc == 1:
          # About to exit strip:
          if qloc != eloc:
            # Exiting from opposite side of strip -- completed a link path:
            LS.append(PS)
          # Either way, we are done with the currenr link path:
          PS = None; eloc = None
  # Since we started outside or on the boundary of the strip, we must be quiescent:
  assert eloc == None
  assert PS == None
  return LS
  # ----------------------------------------------------------------------

def clip_seg_to_strip(p,q, ydir, y0,y1):
  assert y0 < y1
  
  yp = rn.dot(p, ydir)
  yq = rn.dot(q, ydir)

  # Sort the points by Y:
  swap = False
  if yp > yq:
    q,p = p,q; yq,yp = yp,yq; swap = True
    
  assert yp <= yq
  if yp >= y1 or yq <= y0: 
    # Segment outside strip or on its border.
    return None, None, None, None

  if yp >= y0:
    # Starts inside strip:
    ps = p; ploc = None
  else:
    # Enter strip at {y0}
    assert yq > yp
    rp = (y0 - yp)/(yq - yp)
    ps = rn.mix(1-rp, p, rp, q)
    ploc = 0

  if yq <= y1:
    # Ends inside strip:
    qs = q; qloc = None
  else:
    # Exits strip at {y1}
    assert yq > yp
    rq = (y1 - yp)/(yq - yp)
    qs = rn.mix(1-rq, p, rq, q)
    qloc = 1

  if swap:
    return qs, qloc, ps, ploc
  else:
    return ps, ploc, qs, qloc
  # ----------------------------------------------------------------------

def find_point_on_poly(CS,idr,g):

  npt = len(CS)
  assert npt >= 2, "polygonal line must have at least 2 vertices"
  assert g >= 0, "trim distance cannot be negative"

  i_start = None
  q_start = None

  i_prev = (0, npt-1)[idr]
  d_prev = 0 # Distance from endpoint {idr} to point {i_prev} along path.
  for k in range(npt-1):
    i = (k+1, npt-k-2)[idr]
    ei = rn.dist(CS[i], CS[i_prev])
    di = d_prev + ei
    if di > g:
      r = (di - g)/ei
      i_start = i
      q_start = rn.mix(r,CS[i_prev], 1-r,CS[i])
      break
    i_prev = i
    d_prev = di
  assert i_start != None, "trim distance must be less than polygonal length"
  assert i_start >= 0 and i_start < npt # Paranoia.
  
  return i_start, q_start
  # ----------------------------------------------------------------------

def trim_poly(CS, g0, g1):

  npt = len(CS)
  assert npt >= 3, "polygonal must have at least 3 vertices"
  
  i_start = [None,None]
  q_start = [None,None]
  for idr in 0, 1:
    g = (g0, g1)[idr]
    i_start[idr], q_start[idr] = find_point_on_poly(CS,idr,g)

  assert i_start[0] < i_start[1], "polygonal is too short" 

  CSnew = [ q_start[0], ]
  for i in range(i_start[0], i_start[1]+1): CSnew.append(CS[i])
  CSnew.append(q_start[1])

  return CSnew
  # ----------------------------------------------------------------------

def fbomb(die):
  # Returns {False} if die is {False}, else fails.
  assert not die
  return False
  # ----------------------------------------------------------------------

# PLOTTING

def make_canvas(B, dp, scale, frame, grid, nx, ny):

  if dp == None: dp = (0,0)
  szx, szy = rn.box_size(B)

  sz_max = max(szx, szy)
  gstep_big, gstep_sma = choose_grid_steps(sz_max)
  sys.stderr.write("grid steps = %8.4f %8.4f\n" % (gstep_big, gstep_sma))

  if scale == None:
    # Compute the {pyx} scale of the canvas for reasonable images size:
    scale = 15.0/max(nx*szx, ny*szy)
  sys.stderr.write("scale = %8.4f\n" % scale)
  pyx.unit.set(uscale=scale, wscale=scale, vscale=scale, xscale=scale)

  # Set the text rendering engine:
  pyx.text.set(pyx.text.LatexEngine)
  # pyx.text.set(pyx.text.TexEngine)
  # pyx.text.set(pyx.text.UnicodeEngine)

  # Load some packages:
  pyx.text.preamble( ""
    + r"\usepackage{lmodern}" + "\n"
    + r"\usepackage{anyfontsize}" + "\n"
    + r"\usepackage{scalefnt}" + "\n"
    + r"\usepackage{color}" + "\n"
    + r"\usepackage[scaled=1.15]{BOONDOX-calo}" + "\n"
    + r"\newcommand{\su}[1]{\hbox{\scalefont{0.45}\textsf{#1}}}" + "\n"
  )
  #   + r"\usepackage{rotating}" + "\n"
  #   + r"\usepackage{graphicx}" + "\n" 

  # Create the canvas:
  c = pyx.canvas.canvas()
  
  wd_frame = sz_max/500
  inset_frame = 3*wd_frame
  
  wd_grid = sz_max/500
  inset_grid = 5*wd_grid
  dash_grid = ( 3*wd_grid, 2*wd_grid ) # Should suppress dash adjustment and sync dashes to {gstep_big}.
  
  white = pyx.color.rgb.white

  for ix in range(nx):
    for iy in range(ny):
      dpi = rn.add(dp, (ix*szx, iy*szy))

      # White frame to force image bounding box:
      plot_frame(c, B, dpi, white, wd_frame, None, 0.5*wd_frame)
      
      if grid:
        plot_grid(c, B, dpi, None, 0.5*wd_grid, dash_grid, inset_grid, gstep_sma, gstep_sma)
        plot_grid(c, B, dpi, None, 1.0*wd_grid, None,      inset_grid, gstep_big, gstep_big)

      if frame:
        plot_frame(c, B, dpi, None, wd_frame, None, inset_frame)

  # wd_axes = 0.15*min(wd_fill,wd_cont) # Width of jumps and axis lines.

  return c, szx, szy
  # ----------------------------------------------------------------------

def round_box(B, mrg):
  plo = ( floor(B[0][0] - mrg), floor(B[0][1] - mrg)  )
  phi = ( ceil(B[1][0] + mrg), ceil(B[1][1] + mrg)  )
  return (plo, phi)
  # ----------------------------------------------------------------------

def plot_line(c, p, q, dp, clr, wd, dashpat):
  if wd == None or wd <= 0: return
  if clr == None: clr = pyx.color.rgb.black
  
  # Fudge the endpoints if practically equal:
  distpq = rn.dist(p, q)
  distmin = 1.0e-5*wd
  peps = 0 if distpq > distmin else 2*distmin

  if dp != None: p = rn.add(p, dp); q = rn.add(q, dp)

  sty_jround = pyx.style.linejoin.round
  sty_cround = pyx.style.linecap.round
  sty = [ sty_jround, sty_cround, clr ]
  if distpq > wd and dashpat != None: 
    # Convert dash pattern to relative, then adjust it or turn off dashing if line is too short:
    assert type(dashpat) is tuple or type(dashpat) is list
    assert len(dashpat) == 2
    dashpat = adjust_dash_pattern(distpq, dashpat[0], dashpat[1])
  else:
    dashpat = None
  if dashpat != None:
    dashfac = 0.568/wd # Magic scaling factor for dash/gap sizes.
    rdashpat = [dashfac*dashpat[0], dashfac*dashpat[1]]
    # Make the line style dashed:
    sty_dash = pyx.style.dash(rdashpat, 0)
    # sys.stderr.write("in: rdashpat = %s\n" % str(rdashpat))
    # sys.stderr.write("dash = %s\n" % str(sty_dash))
    sty = sty + [ pyx.style.linestyle(sty_cround, sty_dash) ] 
    # sty = sty + [ pyx.style.dash(rdashpat) ]
  sty_wd = pyx.style.linewidth(wd)
  sty += [ sty_wd ]
  c.stroke(pyx.path.line(p[0]-peps, p[1]-peps, q[0]+peps, q[1]+peps), sty)
  # ----------------------------------------------------------------------

def plot_text(c, tx, p, dp, fsize, texclr):
  fspc = fsize  # For now.
  if dp != None: p = rn.add(p, dp)
  if texclr == None:
    tx_tex = (r'\fontsize{%.0f}{%.0f}\selectfont %s' % (fsize,fspc,tx))
  else:
    tx_tex = (r'\fontsize{%.0f}{%.0f}\selectfont \textcolor{%s}{%s}' % (fsize,fspc,texclr,tx))
  c.text(p[0], p[1], tx_tex)
  # ----------------------------------------------------------------------

def plot_box(c, B, dp, clr):
  if clr == None: clr = pyx.color.rgb.black
  if dp == None: dp = (0,0)
  xp = B[0][0] + dp[0]; xsz = B[1][0] - B[0][0]
  yp = B[0][1] + dp[1]; ysz = B[1][1] - B[0][1]
  sty = [ clr ]
  c.fill(pyx.path.rect(xp, yp, xsz, ysz), sty)
  # ----------------------------------------------------------------------

def plot_frame(c, B, dp, clr, wd, dashpat, inset):
  if inset == None: inset = 0
  if wd == None or wd <= 0: return
  
  # Frame coordinate ranges including {dp} and {inset}:
  xmin = B[0][0] + inset; xmax = B[1][0] - inset
  ymin = B[0][1] + inset; ymax = B[1][1] - inset
  
  plot_line(c, (xmin, ymin), (xmax, ymin), dp, clr, wd, dashpat)
  plot_line(c, (xmax, ymin), (xmax, ymax), dp, clr, wd, dashpat)
  plot_line(c, (xmax, ymax), (xmin, ymax), dp, clr, wd, dashpat)
  plot_line(c, (xmin, ymax), (xmin, ymin), dp, clr, wd, dashpat)
  # ----------------------------------------------------------------------

def plot_grid(c, B, dp, clr, wd, dashpat, inset, xstep, ystep):
  if clr == None: clr = pyx.color.grey(0.8)
  if inset == None: inset = 0
  if wd == None or wd <= 0: return

  if xstep == None: xstep = 1
  if ystep == None: ystep = 1

  # Grid coordinate ranges including {inset} but excluding {dp}:
  xmin = B[0][0] + inset; xmax = B[1][0] - inset
  ymin = B[0][1] + inset; ymax = B[1][1] - inset
 
  # Grid line index ranges, with some extra:
  kxmin = int(floor(xmin/xstep)) - 1; kxmax = int(ceil(xmax/xstep)) + 1
  kymin = int(floor(ymin/ystep)) - 1; kymax = int(ceil(ymax/ystep)) + 1

  eps = 1.0e-6*wd
  
  for dkx in range(kxmax-kxmin+1):
    x = (kxmin+dkx)*xstep
    if x >= xmin+eps and x < xmax-eps:
      plot_line(c, (x, ymin), (x, ymax), dp, clr, wd, dashpat)
  for dky in range(kymax-kymin+1):
    y = (kymin+dky)*ystep
    if y >= ymin+eps and y < ymax-eps:
      plot_line(c, (xmin, y), (xmax, y), dp, clr, wd, dashpat)
  # ----------------------------------------------------------------------

def write_plot(c, name, scale_n):
  sys.stderr.write("writing %s.eps ...\n" % name)
  c.writeEPSfile(name + ".eps")
  sys.stderr.write("writing %s.png ...\n" % name)
  convert_eps_to_png(name, scale_n)
  #sys.stderr.write("writing %s.jpg ...\n" % name)
  #convert_eps_to_jpg(name)
  os.remove(name + ".eps")
  sys.stderr.write("done.\n")

def convert_eps_to_png(name, scale_n):
   im = Image.open(name + '.eps')
   im.load(scale=scale_n)
   fig = im.convert('RGBA')
   fig.save(name + '.png', lossless = True)
   im.close()
  # ----------------------------------------------------------------------

def convert_eps_to_jpg(name):
   im = Image.open(name + '.eps')
   im.load(scale=4)
   fig = im.convert('RGB')
   fig.save(name + '.jpg', quality = 95)
   im.close()
  # ----------------------------------------------------------------------

def adjust_dash_pattern(rdist, rlen, rgap):
  assert rlen > 0
  if rgap <= 0: return None # No gaps -- might as well use solid.
  if rdist <= rlen: return None # Not enough space for a single dash.
  
  # Compute the approx number {ngap} of gaps:
  xgap = (rdist - rlen)/(rgap + rlen) # Fractional number of gaps with ideal pattern.
  assert xgap >= 0
  ngap = floor(xgap + 0.5)
  if ngap == 0: ngap = 1
  
  # Check whether {ngap+1} or {ngap-1} may be better:
  nbest = None
  ebest = +inf
  for dn in 0, -1, +1:
    nk = ngap + dn
    if nk >= 0:
      rdk = nk*(rlen+rgap) + rlen
      edk = rdk/rdist if rdk >= rdist else rdist/rdk
      if edk < ebest: ebest = edk; nbest = nk
  assert nbest != None
  
  # Compute the adjustment factor {fac}:
  fac = rdist/(nbest*(rlen+rgap) + rlen)
  minfac = 0.75; maxfac = 1/minfac
  if fac < minfac or fac > maxfac:
    # Would have to shrink or squeeze too much. Only if {rdist} is small, so:
    sys.stderr.write("!! warning: (adjust_dash_pattern} failed fac = %.3f\n" % fac)
    return None
  else:
  # Okeey:
    rdashpat = ( fac*rlen, fac*rgap )
    return rdashpat
  # ----------------------------------------------------------------------

def choose_grid_steps(sz):
  # Ideally the image should be at most 100 times the minor step size.
  nmax = 20
  # Let's set {gstep_sma} to be the largest power of 10 that is less than {sz/nmax}:
  k = floor(log(1.0001*sz/nmax)/log(10))
  gstep_sma = 10.0**k; gstep_big = 10*gstep_sma
  # If {gstep_sma} is too small, increase it:
  if 2*gstep_sma < 1.0001*sz/nmax: gstep_sma = 2*gstep_sma
  return gstep_big, gstep_sma
  # ----------------------------------------------------------------------
  
def trace_colors(nc, Y):
  if Y == None: Y = 0.650
  Ymin = max(0.000, Y - 0.020); Ymax = min(1.000, Y + 0.020)
  Tmin = 0.700; Tmax = 0.950
  return make_colors(nc, Ymin, Ymax, Tmin, Tmax)
  # ----------------------------------------------------------------------

def link_colors(nc, Y):
  if Y == None: Y = 0.700
  Ymin = max(0.000, Y - 0.100); Ymax = min(1.000, Y + 0.100)
  Tmin = 0.900; Tmax = 1.000
  return make_colors(nc, Ymin, Ymax, Tmin, Tmax)
  # ----------------------------------------------------------------------
  
def matter_color(Y):
  if Y == None: Y = 0.870
  YHT = ( Y, 0.300, 0.300 )
  RGB = color.RGB_from_YUV(color.YUV_from_YHS(color.YHS_from_YHT(YHT)))
  clr = pyx.color.rgb( RGB[0], RGB[1], RGB[2] )
  return clr
  # ----------------------------------------------------------------------

def make_colors(nc, Ymin, Ymax, Tmin, Tmax):
  # Used by {trace_colors,link_colors}
  colors = [None]*nc
  if nc == 1:
    colors[0] = pyx.color.rgb(0.200, 0.600, 0.000)
  else:
    RGBS = palette.table(nc, Ymin, Ymax, Tmin, Tmax)
    colors = [ pyx.color.rgb( rgb[0], rgb[1], rgb[2] ) for rgb in RGBS ]
  return colors
  # ----------------------------------------------------------------------
  


