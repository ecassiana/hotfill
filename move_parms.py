# Module to represent width and dynamic parameters of moves.
# Last edited on 2021-10-17 11:40:38 by stolfi

import move_parms_IMP; from move_parms_IMP import Move_Parms_IMP

class  Move_Parms(Move_Parms_IMP):
  # An object {mp} of this class stores some parameters for
  # a move or set of moves. They include the
  # {width(mp)}, the /cruise speed/ {speed(mp)}, the 
  # /acceleration/ at the beginning of the move (and deceleration
  # at the end), and the /nozzle up-down time/ at each end.
  #
  #   * The /nominal width/ {wd}, in millimeters. It is such that
  #   rasters of that width spaced {wd} apart would provide a solid fill
  #   of uniform thickness. This parameter is zero if and only if the
  #   move is a jump.
  #
  #   * The /acceleration/ {ac}, in mm/s^2, of the nozzle at the
  #   beginning of the move, and also the deceleration at the end.
  #
  #   * The /cruise speed/ {sp}, in mm/s, the maximum linear speed of
  #   the nozzle during the execution of the trace.
  #
  #   * The /up-down time/ {ud}, in seconds, being the time required for
  #   nozzle to be raised and/or for the filament to be retracted at the 
  #   start of the jump if the previous move was a trace, and also for 
  #   the inverse operation at end of the jump if the next move is a trace.
  #   This parameter should be nonzero only for jumps.
  #
  # If move (trace or jump) is long enough, the nozzle is assumed to
  # speed up from rest initially with acceleration {+ac} until reaching
  # the cruise speed {sp}, then continue at that speed until near the
  # end, then slow down with acceleration {-ac} until coming to a rest
  # at the endpoint of the move.
  # 
  # The parameters {wd} and {ud} must be non-negative, and the values
  # for {ac} and {sp} must be positive. The acceleration {ac} may be
  # infinite.
  # 
  # Other parameters may be added in the future, such as the desired
  # temperature of the nozzle during extrusion, the filament material or
  # color, or thickness of the trace -- if those parameters cannot be
  # assumed uniform over the whole slice.
  #
  # Putting these parameters in a separate object can save a lot of
  # space, since the same {Move_Parms} object can be shared by 
  # thousands of moves.
  pass
 
def make(wd, ac, sp, ud):
  # Creates a {Move_Parms}  object with the given attributes. 
  return move_parms_IMP.make(wd, ac, sp, ud)

def make_for_jumps(parms):
  # Creates a {Move_Parms}  object with suitable attributes for jumps,
  # as specified in the job parameter dictionary {parms}. See {job_parms.py}.
  return move_parms_IMP.make_for_jumps(parms)

def make_for_contours(parms):
  # Creates a {Move_Parms}  object with suitable attributes for traces in contours,
  # as specified in the job parameter dictionary {parms}. See {job_parms.py}.
  return move_parms_IMP.make_for_contours(parms)

def make_for_fillings(parms):
  # Creates a {Move_Parms}  object with suitable attributes for traces in in fillings,
  # as specified in the job parameter dictionary {parms}. See {job_parms.py}.
  return move_parms_IMP.make_for_fillings(parms)

def width(mp):
  # Returns the nominal width parameter {wd} of the {Move_Parms}  object {mp}.
  return move_parms_IMP.width(mp)

def dynamics(mp):
  # Returns three results: the acceleration {ac}, the cruise speed {sp},
  # and the trace/jump transition penalty {ud}, as stored in the
  # {Move_Parms} object {mp}.
  return move_parms_IMP.dynamics(mp)

def ud_penalty(mp):
  # Returns trace/jump transition penalty {ud}, as stored in the
  # {Move_Parms} object {mp}. 
  return move_parms_IMP.ud_penalty(mp)

def is_jump(mp):
  # Returns {True} if and only if {mp} describes a jump -- that is,
  # if {width(mp) == 0}.
  return move_parms_IMP.is_jump(mp)

def transition_penalty(mp0, mp1):
  # If one of the {Move_Parms} objects {mp0,mp1} describes a 
  # trace and the other describes a jump, returns the
  # trace/jump transition time penalty that would apply
  # if moves with those parameters were consecutive elements of a tool-path. Otherwise 
  # returns 0.
  return move_parms_IMP.transition_penalty(mp0, mp1)

def nozzle_travel_time(dpq, dpm, mp):
  # If {dpm} is {None}, computes the time that the nozzle takes to
  # travel distance {dpq}, including acceleration and deceleration, when
  # executing a move with parameters {mp}. 
  # 
  # Does NOT include the extra time for raising and lowering of the
  # nozzle and/or retracting and re-feeding the filament when
  # transitioning between tracing and jumping.
  #
  # If {dpm} is not {None}, returns instead the time for the nozzle to
  # travel the initial distance {dpm} (which should be between 0 and
  # {dpq} inclusive) along that move.
  return move_parms_IMP.nozzle_travel_time(dpq, dpm, mp)

def show(wr, pref, mp, suff):
  # Writes to file {wr} the parameters stored in {mp}, in a
  # human-readable format, prefixed by the string {pref} (if not {None})
  # and followed by {suff} (ditto). Does not print a newline (unless
  # included in {pref} or {suff}).
  return move_parms_IMP.show(wr, pref, mp, suff)
  
