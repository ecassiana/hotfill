# Tools to generate color palettes and color sets for plots etc..
# Last edited on 2021-10-01 21:58:45 by stolfi

import color
import palette_IMP

def smooth_single(x, x0,RGB0, x1,RGB1, blog):
  # Returns a color (as an RGB triple) that is a smooth function of the
  # float argument {x}, ranging from colors {RGB0} to {RGB1} as {x} ranges
  # from {x0} to {x1}.
  #
  # Specifically, if the float {blog} is zero, values of {x} between {x0}
  # and {x1} are interpolated affinely between {RGB0} and {RGB1}. Values
  # of {x} outside that range are mapped to {RGB0} or {RGB1}, depending on
  # which end value is closer to {x}. The end values {x0} and {x1} must be
  # distinct.
  #
  # If {blog} is nonzero, uses a log scale from {x0} to {x1}, with unit
  # value {blog}. In this case, {abs(x1-x0)} must be greater than
  # {abs(blog)}, and values of {x} such that {abs(x-x0) < abs(blog)} are
  # mapped to {RGB0}.
  return palette_IMP.smooth_single(x, x0,RGB0, x1,RGB1, blog)

def smooth_double(x, xmin, RGBmin, xzer, eps, RGBzer, xmax, RGBmax, ulog):
  # Returns a color (as an RGB triple) that is a smooth 
  # function of the float argument {x}, using two color 
  # scales for two sub-ranges.
  #
  # Requires {xmin < xzer-eps} and {xmax > xzer +eps}. 
  #
  # If the boolean {ulog} is false, values of {x} between {xmin} and
  # {xzer-eps} are affinely interpolated between {RGBmin} and {RGBzer}, while
  # values of {x} between {xzer+eps} and {xmax} are interpolated between
  # {RGBzer} and {RGBmax}.
  #
  # If {ulog} is true, the interpolaton uses {log(|x-xzer|)} instead of {x}.
  #
  # In any case, values of {x} between smaller than {xmin} are mapped to {RGBmin},
  # values larger than {xmax} are mapped to [RGBmax} and values betwen {xzer-eps} 
  # and {xzer+eps} are mapped to {RGBzer}.
  return palette_IMP.smooth_double(x, xmin, RGBmin, xzer, eps, RGBzer, xmax, RGBmax, ulog)

def table(n, Ymin, Ymax, Tmin, Tmax):
  # Returns a table of {n} colors chosen to be 
  # as distinct as possible from each other.  
  #
  # The brightness of the colors will vary between {Ymin} and {Ymax},
  # and their relative saturation (in the {YHT} color system) between
  # {Tmin} and {Tmax}. The {Y} and {T} values must be in {[0.005 _
  # 0.995]}.
  return palette_IMP.table(n, Ymin, Ymax, Tmin, Tmax)
