#! /usr/bin/python3
# Test program for the module {move_parms}.

import move_parms
import hacks
import job_parms
import rn
import pyx
import sys
from math import sqrt, sin, cos, floor, ceil, inf, nan, pi

parms = job_parms.slow()

mp_jump = move_parms.make_for_jumps(parms)
mp_cont = move_parms.make_for_contours(parms)
mp_fill = move_parms.make_for_fillings(parms)

def test_nozzle_travel_time():
  sys.stderr.write("--- testing {nozzle_travel_time} ---\n")

  pfile = open("tests/out/move_parms_TST_times.dat", "w")  # Timing plot file.

  def write_time_dist_plot(jmp,cnt,dpq):
    # Tests nozzle_travel_time for jumps if {jmp} is true,
    # for contours if {jmp} is false and {cnt} is true, or for
    # filling if {jmp} and {cnt} are both false.
    #
    # In case of jump, includes the nozzle up/down time at each end, from {parms}.
    ac = parms['acceleration']
    if jmp:
      wd = 0
      sp = parms['job_jump_speed']
      ud = parms['extrusion_on_off_time']
      mp1 = mp_jump
      assert move_parms.is_jump(mp1)
    elif cnt:
      wd = parms['contour_trace_width']
      sp = parms['job_contour_speed']
      ud = 0
      mp1 = mp_cont
      assert not move_parms.is_jump(mp1)
    else:
      wd = parms['solid_raster_width']
      sp = parms['job_filling_speed']
      ud = 0
      mp1 = mp_fill
      assert not move_parms.is_jump(mp1)
    mp = move_parms.make(wd,ac,sp,ud)
    move_parms.show(sys.stderr, "mp = ", mp, "\n")
    
    wd1 = move_parms.width(mp1)
    ac1,sp1,ud1 = move_parms.dynamics(mp1)
    assert wd1 == wd and ac1 == ac and sp1 == sp and ud1 == ud

    # Plot passage time as a function of position along axis:
    nd = 200
    for id in range(nd+1):
      dpm = (dpq*id)/nd
      tpm = ud + move_parms.nozzle_travel_time(dpq, dpm, mp) 
      if id == nd: tpm += ud
      pfile.write("%8.4f %8.4f\n" % (dpm, tpm))
    pfile.write("\n")

  for jmp, cnt in (False, False), (False, True), (True, False):
    write_time_dist_plot(jmp, cnt, 5.0)
    write_time_dist_plot(jmp, cnt, 0.3)

  pfile.close()
  return
  # ----------------------------------------------------------------------

def test_transition_penalty():
  sys.stderr.write("--- testing {transition_penalty,ud_penalty} ---\n")

  tud = parms['extrusion_on_off_time']
  assert move_parms.ud_penalty(mp_jump) == tud
  assert move_parms.ud_penalty(mp_cont) == 0
  assert move_parms.transition_penalty(mp_cont, mp_jump) == tud
  assert move_parms.transition_penalty(mp_cont, mp_jump) == tud
  assert move_parms.transition_penalty(mp_cont, mp_fill) == 0
  assert move_parms.transition_penalty(mp_jump, mp_jump) == 0
  return
  # ----------------------------------------------------------------------

test_transition_penalty()
test_nozzle_travel_time()

