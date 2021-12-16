# Implementation of module {color}.
# Last edited on 2021-05-30 20:49:16 by jstolfi

# import color
import rn
import sys
from math import sqrt, hypot, floor, ceil, sin, cos, atan2, log, exp, inf, nan, pi

def Y_from_RGB(RGB):
  R = RGB[0]
  G = RGB[1]
  B = RGB[2]
  Y = +0.298911*R +0.586611*G +0.114478*B
  return Y
  # ----------------------------------------------------------------------

def YUV_from_RGB(RGB):
  R = RGB[0]
  G = RGB[1]
  B = RGB[2]
  Y = +0.298911*R +0.586611*G +0.114478*B
  U = -0.147000*R -0.289000*G +0.436000*B
  V = +0.615000*R -0.515000*G -0.100000*B
  return (Y, U, V,)
  # ----------------------------------------------------------------------
  
def RGB_from_YUV(YUV):
  Y = YUV[0];
  U = YUV[1];
  V = YUV[2];
  R = +1.000000*Y -0.001164*U +1.139704*V;
  G = +1.000000*Y -0.395735*U -0.580624*V;
  B = +1.000000*Y +2.030875*U -0.000606*V;
  return (R, G, B,)
  # ----------------------------------------------------------------------

def YHS_from_YUV(YUV):
  Y = YUV[0]
  if Y <= 0:
    YHS = ( 0, 0, 0 )
  else:
    S = hypot(YUV[1], YUV[2])/Y
    H = 0 if S == 0 else atan2(YUV[2],YUV[1])/(2*pi)
    if H >= 1.0: H = H - 1
    if H < 0.0: H = H + 1
    YHS = (Y, H, S)
  return YHS
  # ----------------------------------------------------------------------
 
def YUV_from_YHS(YHS):
  Y = YHS[0]
  S = YHS[2]
  if Y <= 0:
    YUV = ( 0, 0, 0 )
  else:
    if S == 0: 
      YUV = (Y, 0, 0)
    else:
      assert S > 0
      H = YHS[1]
      U = Y*S*cos(H*2*pi)
      V = Y*S*sin(H*2*pi)
      YUV = (Y, U, V)
  return YUV
  # ----------------------------------------------------------------------
  
def YHT_from_YHS(YHS):
  Y, H, S = YHS
  if Y <= 0 or Y >= 1:
    T = 0
  else:
    S_max = max_S_from_YH(Y, H)
    assert S_max > 0
    T = S/S_max
  YHT = (Y, H, T)
  return YHT
    
def YHS_from_YHT(YHT):
  Y, H, T = YHT
  if Y <= 0 or Y >= 1:
    S = 0
  else:
    S_max = max_S_from_YH(Y, H)
    assert S_max > 0
    S = T*S_max
  YHS = (Y, H, S)
  return YHS

def max_S_from_YH(Y, H):
  if Y <= 0 or Y >= 1:
    S_max = 0
  else:
    S_ref= 10 # Guaranteed to be outside the RBG unit cube.
    RGB_ref = RGB_from_YUV(YUV_from_YHS((Y,H,S_ref)))
    f = sat_clip_factor(RGB_ref)
    S_max = f*S_ref
  return S_max
  # ----------------------------------------------------------------------

def sat_clip_factor(RGB):
  Y = Y_from_RGB(RGB)
  if Y < 0 or Y > 1: return 0
  if Y == 0 or Y == 1: return 1
  f = 1.0
  for i in range(3):
    di = RGB[i] - Y
    if di < 0:
      f = min(f, (0 - Y)/di)
    elif di > 0:
      f = min(f, (1 - Y)/di)
    else:
      pass
    ci = Y + f*di
    assert ci > -1.0e-8 and ci < 1 + 1.0e-8
  assert f >= 0 and f <= 1
  return f
  # ----------------------------------------------------------------------

def clip_RGB(RGB):
  if max(RGB) <= 1 and min(RGB) >= 0: 
    RGBc = RGB
  else:
    Y = Y_from_RGB(RGB)
    if Y <= 0:
      RGBc = ( 0, 0, 0 )
    elif Y >= 1:
      RGBc = ( 1, 1, 1 )
    else:
      f = sat_clip_factor(RGB)
      RGBc = tuple( Y + f*(RGB[i]-Y) for i in range(3) )
  return RGBc
  # ----------------------------------------------------------------------

def purish_RGB_colors(Y):
  if Y <= 0:
    L = [(0,0,0)]
  elif Y >= 1:
    L = [(1,1,1)]
  else:
    L = []
    for RGB0, RGB1 in \
      ( (0,0,0), (0,0,1) ), \
      ( (0,0,1), (1,0,1) ), \
      ( (1,0,1), (1,1,1) ), \
      ( (1,0,0), (1,0,1) ), \
      ( (0,0,0), (1,0,0) ), \
      ( (1,0,0), (1,1,0) ), \
      ( (1,1,0), (1,1,1) ), \
      ( (0,1,0), (1,1,0) ), \
      ( (0,0,0), (0,1,0) ), \
      ( (0,1,0), (0,1,1) ), \
      ( (0,0,1), (0,1,1) ), \
      ( (0,1,1), (1,1,1) ), :
      Y0 = Y_from_RGB(RGB0)
      Y1 = Y_from_RGB(RGB1)
      assert Y0 < Y1
      if Y0 <= Y and Y < Y1:
        r = (Y - Y0)/(Y1 - Y0)
        RGBr = rn.mix(1-r, RGB0, r, RGB1)
        L.append(RGBr)
  return L
  # ----------------------------------------------------------------------

def interp_vis_RGB(r, RGB0, RGB1):
  YUV0 = YUV_from_RGB(RGB0)
  YUV1 = YUV_from_RGB(RGB1)
  YUV = interp_vis_YUV(r, YUV0, YUV1)
  RGB = RGB_from_YUV(YUV)
  RGB = clip_RGB(RGB)
  return RGB
  # ----------------------------------------------------------------------
  
def interp_vis_YUV(r, YUV0, YUV1):
  YHS0 = YHS_from_YUV(YUV0)
  YHS1 = YHS_from_YUV(YUV1)
  YHS = interp_vis_YHS(r, YHS0, YHS1)
  YUV = YUV_from_YHS(YHS)
  return YUV

def interp_vis_YHS(r, YHS0, YHS1):
  YHT0 = YHT_from_YHS(YHS0)
  YHT1 = YHT_from_YHS(YHS1)
  YHT = interp_vis_YHT(r, YHT0, YHT1)
  YHS = YHS_from_YHT(YHT)
  return YHS

def interp_vis_YHT(r, YHT0, YHT1):
  assert r >= 0 and r <= 1, "invalid {r}"

  Y0, H0, T0 = YHT0; assert Y0 >= 0 and T0 >= 0
  Y1, H1, T1 = YHT1; assert Y1 >= 0 and T1 >= 0

  # If either color is gray, use the other color's hue:
  if T0 == 0: H0 = H1
  if T1 == 0: H1 = H0

  # Interpolate the brightness in biased log scale:
  Y0 = YHT0[0]; Y1 = YHT1[0] 
  assert Y0 > -1.0e-6 and Y1 > -1.0e-6, "negative {Y}" 
  Ybias = 0.010
  Y = exp((1-r)*log(Y0 + Ybias) + r*log(Y1 + Ybias)) - Ybias
  Y = min(1, max(0, Y))
  
  if Y <= 0:
    YHT = ( 0, 0, 0 )
  elif Y >= 1:
    YHT = ( 1, 0, 0 )
  else:
    # Interpoate angular {UV} hues in linear scale:
    # Clip hue difference to {[-0.5, +0.5]}
    dH = H1 - H0
    if dH > +0.5: dH = dH - 1
    if dH < -0.5: dH = dH + 1
    assert -0.5 <= dH and dH <= +0.5
    H = H0 + r*dH

    # Interpolate relative saturations in biased log scale:
    Tbias = 0.010
    T = exp((1-r)*log(T0 + Tbias) + r*log(T1 + Tbias)) - Tbias
    T = max(0, T)

    # Convert {H,T} back to {U,V}:
    YHT = (Y, H, T)
    
  return YHT
  # ----------------------------------------------------------------------
  

  
