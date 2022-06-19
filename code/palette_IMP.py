# Implementation of module {palette}.

import palette
import color
import rn
import sys
from math import sqrt, hypot, log, exp, sin, cos, floor, ceil, gcd, inf, nan, pi

def smooth_single(x, x0,RGB0, x1,RGB1, blog):
  dtot = x1 - x0
  assert dtot != 0, "end values mut be different"
  r = (x - x0)/dtot
  if r <= 0:
    RGB = RGB0
  elif r >= 1:
    RGB = RGB1
  else:
    if blog != 0:
      # Use log scale:
      runit = abs(blog/dtot)
      if r <= runit:
        r = 0
      else:
        r = log(runit/r)/log(runit)
        assert 0 <= r and r <= 1
    # Interpolate with ratio {r}:
    RGB = color.interp_vis_RGB(r, RGB0, RGB1)
  return RGB
  # ----------------------------------------------------------------------

def smooth_double(x, xmin, RGBmin, xzer, eps, RGBzer, xmax, RGBmax, ulog):
  assert xmin <= xzer-eps
  assert xmax >= xzer+eps
  if abs(x-xzer) < eps:
    RGB = RGBzer
  elif x <= xmin:
    RGB = RGBmin
  elif x >= xmax:
    RGB = RGBmax
  else:
    # Must interpolate. 
    assert x != xzer
    # Define the endpoint away from {xzer}, and sign {sgn} of {x-xzer}
    if x < xzer:
      x1 = xmin; RGB1 = RGBmin; sgn = -1
    elif x > xzer:
      x1 = xmax; RGB1 = RGBmax; sgn = +1
    # Define the endpoint near {xzer}
    RGB0 = RGBzer;
    if ulog:
      x0 = xzer; blog = eps
    else:
      x0 = xzer + sgn*eps; blog = 0
    RGB = smooth_single(x, x0,RGB0, x1,RGB1, blog)
  return RGB
  # ----------------------------------------------------------------------

def table(n, Ymin, Ymax, Tmin, Tmax):
  assert 0 <= Ymin and Ymin <= Ymax and Ymax <= 1.000, "invalid {Ymin,Ymax}"
  assert 0 <= Tmin and Tmin <= Tmax and Tmax <= 1.000, "invalid {Tmin,Tmax}"

  CLRS = []
  phi = (sqrt(5)-1)/2

  def unit_to_range(s, Vbias, Vmin, Vmax):
    # Maps a value {sV} from {[-1 + +1]} to {[Vmin _ Vmax]} in log scale:
    r = (1 + s)/2
    V = exp((1-r)*log(Vmin+Vbias) + r*log(Vmax+Vbias)) - Vbias
    V = max(Vmin, min(Vmax, V))
    return V
    # ....................................................................
  
  # Compute a table whith linearly varying hues:
  for k in range(n):

    # Compute the hue {H}:
    H = k / n
    
    # The saturation {T} and brightness {Y}
    # start as a point {sT,sY} on the unit circle:
    ele = k*phi*2*pi      # Elevation, latitude: angle (radians) on meridian circles.
    sT = cos(ele)
    sY = sin(ele)

    # Map Y from {[-1 _ +1]} to {[Ymin _ Ymax]} in log scale:
    Ybias = 0.010
    Yi = unit_to_range(sY, Ybias, Ymin, Ymax)
    Y = min(1, max(0, Yi)) # Just in case.
    assert Ymin <= Y and Y <= Ymax

    # Compute the {T} of the ideal color, in log scale:
    Tbias = 0.010
    Ti = unit_to_range(sT, Tbias, Tmin, Tmax)
    T = min(1, max(0, Ti)) # Just in case.
    assert Tmin <= T and T <= Tmax

    # Convert to {RGB} coordinates:
    YHT = ( Y, H, T )
    YHS = color.YHS_from_YHT(YHT)
    YUV = color.YUV_from_YHS(YHS)
    RGB = color.RGB_from_YUV(YUV)
    RGB = tuple( min(1, max(0, RGB[k])) for k in range(3) ) # Just in case.

    # Paranoia:
    Ycheck = color.Y_from_RGB(RGB)
    assert Ymin - 1.0e-3 <= Ycheck <= Ymax + 1.0e-3

    CLRS.append(RGB)
    
  # Now scramble the order:
  m0 = int(floor(phi*n + 0.5)) % n
  m = m0;
  mok = None
  while abs(m - m0) < n:
    if m >= 2 and m < n and gcd(m,n) == 1: mok = m; break
    if m < m0: 
      m = 2*m0 - m
    else:
      m = 2*m0 - m - 1
  if mok != None:
    sys.stderr.write("scrambing %d colors with step %d\n" % (n,mok))
    assert gcd(mok,n) == 1
    # Rotate with step {mok}:
    CLRSnew = [ CLRS[(k*mok) % n] for k in range(n) ]
    CLRS = CLRSnew

  return CLRS
  # ----------------------------------------------------------------------
 
