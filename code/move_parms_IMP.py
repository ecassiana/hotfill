# Implementation of module {move_parms}.

import move_parms
import sys
from math import sqrt, sin, cos, floor, ceil, inf, nan, pi

class Move_Parms_IMP:
  
  def __init__(self, wd, ac, sp, ud):
    self.wd = wd
    self.ac = ac
    self.sp = sp
    self.ud = ud
  # ----------------------------------------------------------------------
    
def make(wd, ac, sp, ud):
  assert wd >= 0 and wd < +inf
  assert ud >= 0 and ud < +inf
  assert (wd == 0) or (ud == 0)
  assert ac > 0 
  assert sp > 0 and sp < +inf
  return move_parms.Move_Parms(wd, ac, sp, ud)
  # ----------------------------------------------------------------------
  
def make_for_jumps(parms):
  wd = 0
  ac = parms['acceleration']
  sp = parms['job_jump_speed']
  ud = parms['extrusion_on_off_time']
  return make(wd, ac, sp, ud)
   # ----------------------------------------------------------------------
 
def make_for_contours(parms):
  wd = parms['contour_trace_width']
  ac = parms['acceleration']
  sp = parms['job_contour_speed']
  ud = 0
  return make(wd, ac, sp, ud)
  # ----------------------------------------------------------------------
  
def make_for_fillings(parms):
  wd = parms['solid_raster_width']
  ac = parms['acceleration']
  sp = parms['job_filling_speed']
  ud = 0
  return make(wd, ac, sp, ud)
  # ----------------------------------------------------------------------

def width(mp): 
  return mp.wd
  # ----------------------------------------------------------------------

def dynamics(mp): 
  return mp.ac, mp.sp, mp.ud
  # ----------------------------------------------------------------------
  
def is_jump(mp):
  return mp.wd == 0
  # ----------------------------------------------------------------------
  
def ud_penalty(mp):
  return mp.ud
  # ----------------------------------------------------------------------

def transition_penalty(mp0, mp1):
  if is_jump(mp0) and not is_jump(mp1):
    ac, sp, ud = dynamics(mp0)
    tud = ud
  elif not is_jump(mp0) and is_jump(mp1):
    ac, sp, ud = dynamics(mp1)
    tud = ud
  else:
    tud = 0
  return tud
  # ----------------------------------------------------------------------

def nozzle_travel_time(dpq, dpm, mp):
  ac, sp, ud = dynamics(mp)

  # Figure out max speed {sp_max} during the move, and the
  # acceleration/deceleration and cruise, times and distances
  # {ac_time,ac_dist,cr_time,cr_dist}:

  ac_time = 0 if ac == inf else sp/ac  # Time accelerating or decelerating to{vel}.
  if ac != +inf and sp*ac_time >= dpq:
    # Not enough space to reach cruise speed:
    ac_time = sqrt(dpq/ac)
    ac_dist = dpq/2
    sp_max = ac*ac_time
    cr_dist = 0
    cr_time = 0
  else:
    # Enough distance to reach cruise speed:
    ac_dist = 0 if ac == +inf else sp*ac_time/2 # Distance while accelerating or decelerating.
    sp_max = sp
    cr_dist = dpq - 2*ac_dist # Cruise distance.
    cr_time = cr_dist/sp

  # if dpm != None and dpm <= 0:
  #   sys.stderr.write("ac_dist = %8.2f cr_dist = %8.2f\n" % (ac_dist, cr_dist))
  #   sys.stderr.write("ac_time = %8.2f cr_time = %8.2f\n" % (ac_time, cr_time))
  #   sys.stderr.write("sp_max = %8.2f\n" % sp_max)

  # Compute the passage time {ttot}:
  if dpm == None or dpm >= dpq:
    # Passage is at end of move:
    ttot = ac_time + cr_time + ac_time
  elif dpm > ac_dist + cr_dist:
    # Passage is during deceleration:
    ttot = ac_time + cr_time
    dcm = dpm - (ac_dist + cr_dist)  # Distance after start of decel.
    delta = sp_max*sp_max - 2*dcm*ac
    # sys.stderr.write("%8.2f %8.2f %8.2f %8.2f\n" % (dcm, sp_max, ac, delta))
    tcm = (sp_max - sqrt(delta))/ac # Extra time to passage.
    ttot += tcm
  elif dpm >= ac_dist:
    # Passage is during cruising phase:
    ttot = ac_time
    dam = dpm - ac_dist # Distance after acceleration.
    tam = dam/sp # Extra time to passage.
    ttot += tam
  elif dpm > 0:
    # Passage is during acceleration:
    ttot = 0
    tpm = sqrt(2*dpm/ac) # Extra time to passage.
    ttot += tpm
  else:
    # Passage is at very start:
    ttot = 0
  return ttot
  # ----------------------------------------------------------------------

def show(wr, pref, mp, suff):
  if pref != None: wr.write(pref)
  wr.write("wd: %5.3f mm  ac: %8.3f mm/s^2  sp: %8.3f mm/s  ud: %5.3f s" % ((width(mp),)+dynamics(mp)))
  if suff != None: wr.write(suff)
  return
  # ----------------------------------------------------------------------

