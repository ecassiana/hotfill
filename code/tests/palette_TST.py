#! /usr/bin/python3
# Test program for module {palette}
# Last edited on 2021-10-22 22:30:05 by stolfi

import palette
import color
import hacks
import rn
import pyx
import sys
from math import sqrt, hypot, floor, ceil, sin, cos, inf, nan, pi

def plot_sample(c, k, RGBk, dp, szx_sm, szy_sm, ns):
  #  Plot a sample rectangle with index {k} (in {0..ns-1}) and color {RGBk}.
  # Assumes that the row of samples starts at
  # {dp} and samples have size {szx_sm,szy_sm} with a margin 1 all around the sample row. 
  xlo = dp[0] + 1 + k*szx_sm; xhi = xlo + szx_sm;
  ylo = dp[1] + 1; yhi = ylo + szy_sm
  pyck = pyx.color.rgb( RGBk[0], RGBk[1], RGBk[2] )
  hacks.plot_box(c, ((xlo,ylo), (xhi,yhi)), None, pyck)
  return
  # ----------------------------------------------------------------------

def plot_tic(c, r, dp, szx_sm, szy_sm, ns):
  #  Plots a tic mark at position {r} between the centers
  # of samples {k=2} and {k=ns-3}. 
  kf = r*(ns-5) + 2 # Fractional palette index corresp. {r}
  xlo = 1 + kf*szx_sm; xhi = xlo + szx_sm; xmd = (xlo + xhi)/2
  ylo = 1; yhi = ylo + szy_sm; ymd = (ylo + yhi)/2
  yr0 = ymd - 0.2; yr1 = ymd + 0.2 
  wd_tic = 0.05
  dashed = None
  hacks.plot_line(c, (xmd, yr0), (xmd, yr1), dp, pyx.color.rgb.black, wd_tic, dashed)
  return
  # ----------------------------------------------------------------------

def test_smooth_single(c, ns, dp, szx_sm, szy_sm):
  # Plots two rows of {ns} samples of {palette.smooth_single} on {c}
  # starting at position {dp}. Each sample has size {szx_sm,szy_sm} with
  # 1 unit of space below and above the two rows.
  
  sys.stderr.write("--- testing{smooth_single} ---\n")

  t0 = 17; RGB0 = ( 0.000, 0.500, 0.000 )
  t1 = -3; RGB1 = ( 1.000, 0.000, 1.000 )
  blog = 2
  
  # Testing {smooth_single}:
  for ulog in False, True:
    for k in range(ns):
      rk = (k-2)/(ns-5)
      tk = (1-rk)*t0 + rk*t1
      RGBk = palette.smooth_single(tk, t0,RGB0, t1,RGB1, blog*int(ulog))
      sys.stderr.write("k = %d RGBk = ( %5.3f %5.3f %5.3f )\n" % ((k,) + RGBk)) 
      plot_sample(c, k, RGBk, dp, szx_sm, szy_sm, ns)
    plot_tic(c, 0.0, dp, szx_sm, szy_sm, ns)
    plot_tic(c, 1.0, dp, szx_sm, szy_sm, ns)
    if ulog: plot_tic(c, abs(blog/(t1 - t0)), dp, szx_sm, szy_sm, ns)
    dp = rn.add(dp, (0, szy_sm))

  return
  # ----------------------------------------------------------------------
  
def test_smooth_double(c, ns, dp, szx_sm, szy_sm):
  # Plots two rows of {ns} samples of {palette.smooth_double} on {c}
  # starting at position {dp}. Each sample has height {szx_sm, szy_sm}
  # with 1 unit of space below and above the two rows.
  
  sys.stderr.write("--- testing{smooth_double} ---\n")
  
  tmin =  -9; RGBmin = ( 0.000, 0.167, 1.000 )            # Y = 0.200
  tzer =  +4; RGBzer = ( 0.350, 0.350, 0.350 ); eps = 0.5 # Y = 0.350
  tmax = +10; RGBmax = ( 1.000, 0.500, 0.000 )            # Y = 0.600

  # Testing {smooth_double}:
  for ulog in False, True:
    for k in range(ns):
      rk = (k-2)/(ns-5)
      tk = (1-rk)*tmin + rk*tmax
      RGBk = palette.smooth_double(tk, tmin, RGBmin, tzer, eps, RGBzer, tmax, RGBmax, ulog)
      sys.stderr.write("k = %d RGBk = ( %5.3f %5.3f %5.3f )\n" % ((k,) + RGBk)) 
      plot_sample(c, k, RGBk, dp, szx_sm, szy_sm, ns)
    plot_tic(c, 0.0, dp, szx_sm, szy_sm, ns)
    plot_tic(c, 1.0, dp, szx_sm, szy_sm, ns)
    plot_tic(c, (tzer-eps-tmin)/(tmax - tmin), dp, szx_sm, szy_sm, ns)
    plot_tic(c, (tzer-tmin)/(tmax - tmin), dp, szx_sm, szy_sm, ns)
    plot_tic(c, (tzer+eps-tmin)/(tmax - tmin), dp, szx_sm, szy_sm, ns)
    dp = rn.add(dp, (0, szy_sm))

  return
  # ----------------------------------------------------------------------
  
def test_table(c, ns, dp, szx_sm, szy_sm, Ymin, Ymax, Tmin, Tmax):
  # Plots {ns} samples of {palette.table} on {c} starting at position {dp}.
  # Each sample has height {szx_sm, szy_sm} with 1 unit of space below and above.
  
  sys.stderr.write("--- testing{table} ---\n")
  
  sys.stderr.write("Y = [ %5.3f _ %5.3f ] S = [ %5.3f _ %5.3f ]\n" % (Ymin,Ymax,Tmin,Tmax))
  RGBS = palette.table(ns, Ymin, Ymax, Tmin, Tmax)
  for k in range(ns):
    RGBk = RGBS[k]
    sys.stderr.write("k = %d RGBk = ( %5.3f %5.3f %5.3f )" % ((k,) + RGBk)) 
    plot_sample(c, k, RGBk, dp, szx_sm, szy_sm, ns)
    YUVk = color.YUV_from_RGB(RGBk)
    sys.stderr.write("  YUVk = ( %5.3f %+6.3f %+6.3f )" % YUVk) 
    YHSk = color.YHS_from_YUV(YUVk)
    sys.stderr.write("  YHSk = ( %5.3f %5.3f %5.3f )\n" % YHSk) 
    assert Ymin-1.0e-5 <= YUVk[0] and YUVk[0] <= Ymax+1.0e-5
  plot_tic(c, 0.0, dp, szx_sm, szy_sm, ns)
    
  return 
  # ----------------------------------------------------------------------

def test_all():

  ns = 25 # Number of color samples per row
  
  tests_table = (
      (0.100, 0.100, 1.000, 1.000),
      (0.100, 0.100, 0.500, 0.700),
      (0.100, 0.500, 0.500, 1.000),
      (0.100, 0.900, 0.300, 0.500),
      
      (0.300, 0.300, 1.000, 1.000), 
      (0.300, 0.300, 0.500, 0.700), 
      (0.300, 0.600, 0.500, 1.000), 
      (0.300, 0.600, 0.300, 0.500),

      (0.500, 0.500, 1.000, 1.000), 
      (0.500, 0.500, 0.500, 0.700), 
      (0.500, 0.800, 0.500, 1.000), 
      (0.500, 0.800, 0.300, 0.500), 

      (0.700, 0.700, 1.000, 1.000), 
      (0.700, 0.700, 0.500, 0.700), 
      (0.700, 1.000, 0.700, 1.000), 
      (0.700, 1.000, 0.300, 0.300), 
    )

  nt = len(tests_table) # Number of {palette.table} tests.
  
  szx_sm = 4; szy_sm = 4
  ysep = 0.5
  dp = (0, 0)

  xtot = 1 + ns*szx_sm + 1
  ytot = 1 + 2*(2*szy_sm + ysep) + nt*(szy_sm + ysep) - ysep + 1
  B = ((0,0), (xtot, ytot))
  c,szx,szy = hacks.make_canvas(hacks.round_box(B,0.5), dp, None, False, False, 1, 1)

  dpk = dp
  test_smooth_single(c, ns, dpk, szx_sm, szy_sm); dpk = rn.add(dpk, (0, 2*szy_sm + ysep)) 
  test_smooth_double(c, ns, dpk, szx_sm, szy_sm); dpk = rn.add(dpk, (0, 2*szy_sm + ysep)) 
    
  for Ymin, Ymax, Tmin, Tmax in tests_table:
    test_table(c, ns, dpk, szx_sm, szy_sm, Ymin, Ymax, Tmin, Tmax); dpk = rn.add(dpk, (0, szy_sm + ysep)) 

  hacks.write_plot(c, "tests/out/palette_TST_")
  return
  # ----------------------------------------------------------------------
  
test_all()
