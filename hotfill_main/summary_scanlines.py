# Implementation of the module {hotfill}.
# Last edited on 2021-11-25 01:08:56 by stolfi

import hotfill
import bandpath
import raster
import path
import contact
import move
import move_parms
import hacks
import rn
import sys
import time
from math import sqrt, sin, cos, log, exp, floor, ceil, inf, nan, pi
import txt_read
import os

# --------------------------------------------------------------- #

def solve(OPHS):
  xdir = (1,0)
  ydir = (0,1)
  ystep,yphase = raster.get_spacing_and_phase(OPHS, xdir, ydir)
  OPHS = raster.sort_by_scanline(OPHS, xdir, ydir, ystep, yphase)
  SCS = raster.separate_by_scanline(OPHS, xdir, ydir, ystep, yphase)

  x_min = None
  x_max = None
  y_min = None
  y_max = None

  for ph in OPHS:
    p, q = path.endpoints(ph)

    x_min_aux = min(p[0], q[0])
    x_max_aux = max(p[0], q[0])
    y_min_aux = min(p[1], q[1])
    y_max_aux = max(p[1], q[1])

    if x_min == None or x_min > x_min_aux:
      x_min = x_min_aux
    if y_min == None or y_min > y_min_aux:
      y_min = y_min_aux
    if x_max == None or x_max < x_max_aux:
      x_max = x_max_aux
    if y_max == None or y_max < y_max_aux:
      y_max = y_max_aux

  return len(SCS), x_max-x_min, y_max-y_min

# --------------------------------------------------------------- #

def read_txt(fname_in):
  ##
  # Some arbitrary dynamics parameters:
  acc = 3000      # Acceleration/deceleration for moves and jumps (mm/s^2).
  csp_trace = 40  # Cruise speed for traces (mm/s).

  angle = 0 if "_0.txt" in fname_in else 1.5708
  shift = (0, 0)

  # Move parameters matching the raster spacings in the file:
  wd_cont = 0.40; mp_cont = move_parms.make(wd_cont, acc, csp_trace, 0.0)
  wd_fill = 0.40; mp_fill = move_parms.make(wd_fill, acc, csp_trace, 0.0)
  wd_link = 0.40; mp_link = move_parms.make(wd_link, acc, csp_trace, 0.0)

  rd = open(fname_in, 'r')
  smptol = 0.0
  fixdir = True
  OCRS, OPHS, OLKS, CTS, Z = txt_read.read(rd, mp_cont, mp_fill, mp_link, angle, shift, fixdir, smptol)
  rd.close()

  nsc = solve(OPHS)
  
  return nsc

# --------------------------------------------------------------- #

def found_keys(fname_in):
  split_list = fname_in.split('_')
  key_part = split_list[0].split('/')[-1]
  key_angle = split_list[2].replace('.txt', '')
    
  return key_part, key_angle

# --------------------------------------------------------------- #

def main(infolder):
  infolder  = "tests/in/"
  dict_scanline = dict()

  lst = os.listdir(infolder)
  lst.sort()

  for file in lst:
    if '.txt' in file:
      key_part, key_angle = found_keys(file)
      
      fname_in = infolder + file
      nsc, x, y = read_txt(fname_in)

      if key_part not in dict_scanline.keys():
        dict_scanline[key_part] = dict()

      dict_scanline[key_part][key_angle] = dict()
      dict_scanline[key_part][key_angle]['nsc'] = nsc
      dict_scanline[key_part][key_angle]['x'] = round(x, 1)
      dict_scanline[key_part][key_angle]['y'] = round(y, 1)
  
  return dict_scanline