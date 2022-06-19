# Tools to handle colors as triples of floats.

import color_IMP

# For this module, a color is represented mainly by its /RGB
# coordinates/ a triple of floats {(R,G,B)} that are the relative
# amounts of the NTSC red, green, and blue primary lights that are
# perceived as that color; with each factor scaled so that the RGB
# triple {(1,1,1)} is perceived as "white".
#
# Note that in many displays and image files each component of a
# specified "RGB" triple undergoes a non-linear encoding or decoding
# ("gamma correction"). This module assumes that RGB triples are in
# /linear/ color space: where any "gamma" encoding has been undone, so
# that each component of the triple is proportional to the luminous
# power (energy per second per solid angle) of the corresponding primary
# light. Improper correction of "gamma" encoding may change the hue and
# saturation of the perceived color, as well as its brightness.
  
# RGB-YUV CONVERSION

def Y_from_RGB(RGB):
  # Returns the {Y} (/brightness/) value of the RGB color coordinate
  # triple {RGB=(R,G,B)}.
  #
  # The brightness {Y} is a linear combination of the {R}, {G}, and {B}
  # coordinates that is a monotonic funtion of the apparently
  # "brightness", or "lightness" as perceived by the typical human
  # eye. The {Y} value is positive for all pysically realizable colors,
  # is 0 for {(0,0,0)}, and is 1 for the RGB "white" {(1,1,1)}. More
  # generally, the RGB triple {(Y,Y,Y)} represents a neutral gray of
  # brightness {Y}.
  #
  # The formula for {Y} accounts for the fact that the pure colors
  # {(1,0,0)}, {(0,1,0)}, and {(0,0,1)} have vastly different perceived
  # brightness. 
  #
  # This procedure assumes that the {R}, {G}, and {B} coordinates are in
  # linear scale. The visual perception of brightness is non-linear,
  # however; it is more like {log(Y+bias)} for some positive {bias}.
  return color_IMP.Y_from_RGB(RGB)

def YUV_from_RGB(RGB):
  # Returns the YUV coordinate triple {YUV=(Y,U,V)} of a color, given its
  # RGB color coordinate triple {RGB=(R,G,B)}.
  #
  # The YUV coordinates are a linear mapping of the RGB coordinates. The
  # {Y} axis is the color's brightness, as described in {Y_from_RGB}.
  # The {U} and {V} are /chroma/ components, respecively "blueness"
  # (negative: "yellowness") and "redness" (negative: "greenness").
  # Changing the {U} and/or {V} coordinates does not change the color's
  # brightness {Y}. The {U} and {V} axes are nearly orthogonal in RGB
  # space, and are scaled roughly in proportion of visual impact.
  #
  # The entire unit RGB cube, mapped to{YUV} space, is contained in the
  # upended cone about the {Y} axis, apex at the origin, height 1, and base
  # radius 3.908.
  return color_IMP.YUV_from_RGB(RGB)
 
def RGB_from_YUV(YUV):
  # Returns the RGB color coordinate triple {RGB=(R,G,B)} of a color, given its
  # YUV coordinate triple {YUV=(Y,U,V)}.
  return color_IMP.RGB_from_YUV(YUV)
  
# YHS-YUV CONVERSION

def YHS_from_YUV(YUV):
  # Returns the YHS coordinate triple {YHS=(Y,H,S)} of a color, given its YUV coordinate
  # triple {YUV=(Y,U,V)}.
  #
  # The hue coordinate {H} is the angle of the vector {(U,V)} with the {U} axis
  # towards the {V} axis, expressed as a fraction of a turn.
  #
  # The saturation {S} is the modulus of the vector {(U,V)} divided by 
  # the brightness {Y}.In this space, {Y} and {S} are always non-negative, and {S} is zero if
  # {Y} is zero.
  #
  # If the input {Y} is negative or zero, the result is {{(0,0,0)}.
  # Otherwise, if the input {U} and {V} are zero, the hue {H} is zero by
  # definition. The procedure returns {H} normalized in the range {[0_1)}.
  return color_IMP.YHS_from_YUV(YUV)

def YUV_from_YHS(YHS):
  # Returns the YUV color coordinate triple {YUV=(Y,U,V)} of a color,
  # given the YHS coordinates {YHS=(Y,H,S)}. See {YHS_from_YUV(YUV)}.
  #
  # The input {H} can be any float; only the fractional part is
  # considered. If the input {Y} is zero or negative, {S} and {H} are
  # ignored and the result is {(0,0,0)}.
  return color_IMP.YUV_from_YHS(YHS)

# YUV-YHT CONVERSION

def YHT_from_YHS(YHS):
  # Returns the YHT coordinate triple {YHT=(Y,H,T)} of a color, given its YHS
  # coordinate triple (YHS=(Y,H,S)}.
  #
  # The {T} coordinate is the /relative saturation/, a rescaling of the
  # {S} coordinate so that {T=1} means the maximum saturation {S} that
  # will give a point inside the RGB unit cube in RGB space. Thus the RGB
  # unit cube is mapped to the cylinder of unit radius and height in
  # cylindrical coordinates. 
  #
  # Whereas the mappings between {RGB},{YUV} and {YHS} are smooth
  # functions, except at edge cases, the {T} coordinate is a contunuous
  # but non-differentiable function of those other coordinates.
  return color_IMP.YHT_from_YHS(YHS)

def YHS_from_YHT(YHT):
  # Returns the YHS coordinate triple (YHS=(Y,H,S)} of a color, given
  # its YHT coordinate triple {YHT=(Y,H,T)}. See {YHT_from_YHS}.
  return color_IMP.YHS_from_YHT(YHT)

# RGB CUBE CLIPPING
  
def max_S_from_YH(y, h):
  # Returns the maximum {S} coordinate in the {YHS} system for any 
  # RGB color in the RGB unit cube that has the given {Y} and {H}
  # coordinates.
  #
  # If the input {Y} is negative or zero, or 1 or more, returns 0.
  return color_IMP.max_S_from_YH(y, h)
  
def sat_clip_factor(RGB):
  # Finds the factor {f} in {[0 _ 1]} by which the chroma component of
  # {RGB} must be reduced, without changing its brightness {Y} or hue
  # {H}, so that it will lie inside the unit RGB cube.
  #
  # Returns 0 if the input brightness {Y} is outside the range {[0 _ 1]}.
  return color_IMP.sat_clip_factor(RGB)
  
def clip_RGB(RGB):
  # If the given RGB color coordinate triple {RGB} is in the unit cube
  # {[0_1]^3}, returns it unchanged. Otherwise reduces the saturation,
  # preserving its hue and brigtness, so as to bring it into the cube.
  #
  # If the input {Y} is greater than 1, the result is {(1,1,1)}.
  # If it is less than 0, the result is {(0,0,0)}.
  return color_IMP.clip_RGB(RGB)

def purish_RGB_colors(Y):
  # Returns a list of the colors with brightness {Y} that lie on the edges of
  # the unit RGB cube.  
  #
  # If {Y} is zero or negative, the result is {[(0,0,0)]}.
  # If {Y} is 1 or greater, the result is {[(1,1,1)]}. Otherwise
  # the result will be a list of between 3 and 6 colors, in order
  # of increasing hue.
  return color_IMP.purish_RGB_colors(Y)

# COLOR INTERPOLATION

def interp_vis_RGB(r, RGB0, RGB1):
  # Interpolates {r} of the way from colors described by the RGB triples 
  # {RGB0} to {RGB1}, in a hopefully visually uniform way. 
  # The colors are converted to {YHT} space and smoothly interpolated there.
  return color_IMP.interp_vis_RGB(r, RGB0, RGB1)
  
def interp_vis_YUV(r, YUV0, YUV1):
  # Interpolates {r} of the way from colors described by the YUV triples 
  # {YUV0} to {YUV1}, hopefully in a visually uniform way.
  # The colors are converted to {YHT} space and smoothly interpolated there.
  return color_IMP.interp_vis_YUV(r, YUV0, YUV1)
  
def interp_vis_YHS(r, YHS0, YHS1):
  # Interpolates {r} of the way from colors described by the YHS triples 
  # {YHS0} to {YHS1}, hopefully in a visually uniform way.
  # The colors are converted to {YHT} space and smoothly interpolated there.
  return color_IMP.interp_vis_YHS(r, YHS0, YHS1)
  
def interp_vis_YHT(r, YHT0, YHT1):
  # Interpolates {r} of the way from colors described by the YHT triples 
  # {YHT0} to {YHT1}, hopefully in a visually uniform way.
  #
  # The RELATIVE saturation {T} of the result will interpolate between 
  # those of the two colors.
  return color_IMP.interp_vis_YHT(r, YHT0, YHT1)
