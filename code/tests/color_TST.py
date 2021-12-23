#! /usr/bin/python3
# Test program for the {color} module.
# Last edited on 2021-05-31 07:14:39 by jstolfi

import color
import rn
import sys
from math import sqrt, hypot, floor, ceil, sin, cos, inf, nan, pi

# Pure colors:
pures = \
  ( ( 0.000, 0.000, 0.000 ),
    ( 1.000, 0.000, 0.000 ),
    ( 0.000, 1.000, 0.000 ),
    ( 0.000, 0.000, 1.000 ),
    ( 1.000, 1.000, 0.000 ),
    ( 1.000, 0.000, 1.000 ),
    ( 0.000, 1.000, 1.000 ),
    ( 1.000, 1.000, 1.000 ),
  )

grays = \
  (
    ( 0.001, 0.001, 0.001 ),
    ( 0.010, 0.010, 0.010 ),
    ( 0.500, 0.500, 0.500 ),
    ( 0.990, 0.990, 0.990 ),
    ( 0.999, 0.999, 0.999 ),
  )

def closenough(c1, c2):
  # Returns true if and only if color {c1} is close enough to color {c2}.
  d12 = rn.dist(c1,c2)
  return d12 <= 1.0e-5, "colors do not match"

def test_RGB_YUV_YHS_YHT():
  sys.stderr.write("--- testing RGB/YUV/YHS/YHT conversion ---\n")
  
  for RGB in pures + grays: 
    
    YUV = color.YUV_from_RGB(RGB)
    YHS = color.YHS_from_YUV(YUV)
    YHT = color.YHT_from_YHS(YHS)

    Y = color.Y_from_RGB(RGB)
    assert Y == YUV[0]
    if min(RGB) == max(RGB): assert abs(Y - RGB[0]) < 1.0e-8
    
    if Y <= 0:
      assert YUV[1] == 0
      assert YUV[2] == 0
      assert YHS[2] == 0
      assert YHT[2] == 0
    elif Y >= 1:
      assert YHT[2] == 0
    elif Y > 1.0e-6 and Y < 1 - 1.0e-6:
      assert abs(YUV[0] - Y) <= 1.0e-6
      assert abs(YHS[0] - Y) <= 1.0e-6
      assert abs(YHT[0] - Y) <= 1.0e-6
      if max(RGB) == 1 or min(RGB) == 0:
        # A color on the boundary of the RGB cube:
        # sys.stderr.write("RGB = ( %.6f %.6f %.6f )\n" % RGB)
        # sys.stderr.write("YHT = ( %.6f %.6f %.6f )\n" % YHT)
        assert abs(YHT[2] - 1) < 1.0e-5
    
    YHS2 = color.YHS_from_YHT(YHT)
    assert closenough(YHS, YHS2)
    
    YUV2 = color.YUV_from_YHS(YHS)
    assert closenough(YUV, YUV2)
    
    RGB2 = color.RGB_from_YUV(YUV)
    assert closenough(RGB, RGB2)
      
  return
  # ----------------------------------------------------------------------

def test_clip():
  sys.stderr.write("--- testing RGB clipping conversion ---\n")
  
  # Saturation scaling:
  nH = 11
  for Y in ( 0.00, 0.25, 0.99, 1.00 ):
    gry = (Y, Y, Y)
    for kH in range(nH):
      H = kH/nH
      ang = 2*pi*kH/nH
      if Y == 0:
        Slist = ( 0.000, )
      else:
        Smax = color.max_S_from_YH(Y, H)
        # sys.stderr.write("Smax = %s\n" % str(Smax))
        Slist = ( 0.00, 0.20, 0.80, 2.00, 4.00, Smax, )
      for S in Slist:
        YHS0 = (Y, H, S)
        # sys.stderr.write("YHS0 = ( %.6f %.6f %.6f )\n" % YHS0)

        YUV0 = color.YUV_from_YHS(YHS0)
        # sys.stderr.write("YUV0 = ( %.6f %.6f %.6f )\n" % YUV0)

        RGB0 = color.RGB_from_YUV(YUV0)
        # sys.stderr.write("RGB0 = ( %.6f %.6f %.6f )\n" % RGB0)
        
        YHS1 = color.YHS_from_YUV(YUV0)
        # sys.stderr.write("YHS1 = ( %.6f %.6f %.6f )\n" % YHS1)
        if Y > 0 and S > 0: assert rn.dist(YHS1, YHS0) < 1.0e-6

        U = Y*S*cos(ang)
        V = Y*S*sin(ang)
        YUV1 = (Y, U, V)
        # sys.stderr.write("YUV1 = ( %.6f %.6f %.6f )\n" % YUV1)
        if Y > 0: assert rn.dist(YUV1, YUV0) < 1.0e-6
        
        cmax = max(RGB0); cmin = min(RGB0)
        ovf0 = max(0, max(cmax - 1, -cmin)) # Amount of overflow from unit cube.
        f = color.sat_clip_factor(RGB0)
        # sys.stderr.write("cmin = %.8f cmax = %.8f ovf0 = %.8f f = %.8f\n" % (cmin, cmax, ovf0, f))
        assert 0 <= f and f <= 1
        if ovf0 <= 0: assert f == 1
        RGB2 = rn.mix(1, gry, f, rn.sub(RGB0, gry))
        ovf1 = max(0, max(max(RGB2) - 1, -min(RGB2)))
        # sys.stderr.write("ovf1 = %.6f\n" % ovf1)
        if f == 1: assert ovf1 <= 0 and rn.dist(RGB0, RGB2) < 1.0e-6

        if Y > 0 and S == Smax:
          # Color must be on surface of the {RGB} unit cube:
          assert abs(cmax - 1) < 1.0e-5 or abs(cmin) < 1.0e-5

        RGB3 = color.clip_RGB(RGB0)
        Y2 = color.Y_from_RGB(RGB3)
        if Y >= 0 and Y <= 1:
          assert abs(Y2 - Y) < 1.0e-5
          ZUV0 = rn.sub(RGB0, gry)
          ZUV1 = rn.scale(f, ZUV0)
          ZUV3 = rn.sub(RGB3, gry)
          # sys.stderr.write("ZUV0 = ( %.6f %.6f %.6f )\n" % ZUV0)
          # sys.stderr.write("ZUV1 = ( %.6f %.6f %.6f )\n" % ZUV1)
          # sys.stderr.write("ZUV3 = ( %.6f %.6f %.6f )\n" % ZUV3)
          assert rn.dist(ZUV1, ZUV3) < 1.0e-5
        
    # Test {purish_RGB_colors}:
    for Y in \
      0.000, 0.001, \
      0.114, 0.115, \
      0.298, 0.299, \
      0.413, 0.414, \
      0.586, 0.587, \
      0.701, 0.702, \
      0.885, 0.886, \
      0.999, 1.000, :
      L = color.purish_RGB_colors(Y)
      Hant = None; dh = 0
      for RGB in L:
        Yr = color.Y_from_RGB(RGB)
        assert (Yr - Y) < 1.0e-6
        C = list(RGB); C.sort()
        if C[0] == 0 and C[1] == 0:
          r = C[2]
        elif C[0] == 0 and C[2] == 1:
          r = C[1]
        elif C[1] == 1 and C[2] == 1:
          r = C[0]
        else:
          assert False, ("not purish: " + str(RGB) + " " + str(C))
        assert 0 <= r and r <= 1
        YHS = color.YHS_from_YUV(color.YUV_from_RGB(RGB))
        H = YHS[1]
        if Hant == None:
          Hant = H
        else:
          if H < Hant: dH = 1 # Made a full turn?
          assert H + dH > Hant, "hues out of order"
          Hant = H + dH
  return      
  # ----------------------------------------------------------------------
  
def test_interp():
  sys.stderr.write("--- testing visual interpolation ---\n")

  def do_test_interp(YHS0, YHS1, r):
    YHSi = color.interp_vis_YHS(r, YHS0, YHS1)
    sys.stderr.write("YHSi = ( %.8f %.8f %.8f)\n" % YHSi)

    YUVi = color.interp_vis_YUV(r, YUV0, YUV1)
    sys.stderr.write("YUVi = ( %.8f %.8f %.8f)\n" % YUVi)

    RGBi = color.interp_vis_RGB(r, RGB0, RGB1)
    sys.stderr.write("RGBi = ( %.8f %.8f %.8f)\n" % RGBi)

    Yi, Hi, Si = YHSi
    assert Yi >= 0 and Yi <= 1
    if Si > 1.0e-6 and Y0 > 0 and Y0 < 1 and Y1 > 0 and Y1 < 1:
      
      # Check saturations relative to the RGB cube:
      Sm0 = color.max_S_from_YH(Y0, H0); Sr0 = S0/Sm0
      sys.stderr.write("S0 = %.8f Sm0 = %.8f Sr0 = %.8f\n" % (S0,Sm0,Sr0))
      
      Sm1 = color.max_S_from_YH(Y1, H1); Sr1 = S1/Sm1
      sys.stderr.write("S1 = %.8f Sm1 = %.8f Sr1 = %.8f\n" % (S1,Sm1,Sr1))
      
      Smi = color.max_S_from_YH(Yi, Hi); Sri = Si/Smi
      sys.stderr.write("Si = %.8f Smi = %.8f Sri = %.8f\n" % (Si,Smi,Sri))
      
      assert min(Sr0,Sr1)-1.0e-6 <= Sri and Sri <= max(Sr0,Sr1)+ 1.0e-6

    return RGBi
    # ----------------------------------------------------------------------
  
  for RGB0 in pures + grays: 
    for RGB1 in pures + grays: 
      sys.stderr.write("\n- - - - - - - - - - - - - - - - - - - - \n")
      sys.stderr.write("RGB0 = ( %.8f %.8f %.8f)\n" % RGB0)
      sys.stderr.write("RGB1 = ( %.8f %.8f %.8f)\n" % RGB1)

      YUV0 = color.YUV_from_RGB(RGB0)
      YUV1 = color.YUV_from_RGB(RGB1)
      sys.stderr.write("YUV0 = ( %.8f %.8f %.8f)\n" % YUV0)
      sys.stderr.write("YUV1 = ( %.8f %.8f %.8f)\n" % YUV1)

      YHS0 = color.YHS_from_YUV(YUV0); Y0, H0, S0 = YHS0
      YHS1 = color.YHS_from_YUV(YUV1); Y1, H1, S1 = YHS1
      sys.stderr.write("YHS0 = ( %.8f %.8f %.8f)\n" % YHS0)
      sys.stderr.write("YHS1 = ( %.8f %.8f %.8f)\n" % YHS1)
      for r in 0.00, 0.25, 0.99, 1.00:
        RGBi = do_test_interp(YHS0, YHS1, r)
        if RGB0 == RGB1: assert rn.dist(RGBi, RGB0) < 1.0e-6
        sys.stderr.write("\n")


  sys.stderr.write("??? NEEDS MORE TESTS ???\n")
  
  return      
  # ----------------------------------------------------------------------

test_RGB_YUV_YHS_YHT()
test_clip()
test_interp()
