# Implementation of module {move_example}
# Last edited on 2021-03-21 09:54:08 by jstolfi

import move_example
import move
import move_parms
import hacks
import rn
import pyx 
from math import nan, inf, sqrt
import sys

def rectangle_rasters(plo, axis, n, sz, step, mp_trace):
  TRS = []
  for itr in range(n):
    if axis == 0:
      p = (plo[0], plo[1] + itr*step)
      q = (p[0] + sz, p[1])
    elif axis == 1:
      p = (plo[0] + itr*step, plo[1])
      q = (p[0], p[1] + sz)
    else:
      assert False, "invalid {axis}"
    mv = move.make(p, q, mp_trace)
    move.set_name(mv, "T%d" % itr)
    TRS.append(mv)
  return TRS
  # ----------------------------------------------------------------------

def rectangle_links(plo, axis, n, sz, step, mp_trace):
  assert axis == 0 or axis == 1, "invalid {axis}"
  LKS0 = []; LKS1 = []
  for itr in range(n-1):
    for side in range(2):
      a = ( side*sz, itr*step )
      b = ( side*sz, a[1] + step)
      p0 = (plo[0] + a[axis], plo[1] + a[1-axis])
      p1 = (plo[0] + b[axis], plo[1] + b[1-axis])
      mv = move.make(p0, p1, mp_trace)
      move.set_name(mv, "L" + "ab"[side] + str(itr))
      (LKS0, LKS1)[side].append(mv)
  return LKS0, LKS1
  # ----------------------------------------------------------------------

def rectangle_jumps(plo, axis, n, sz, step, mp_jump):
  assert axis == 0 or axis == 1, "invalid {axis}"
  JMS0 = []; JMS1 = []
  for itr in range(n-1):
    for side in range(2):
      a = ( (1-side)*sz, itr*step )
      b = ( side*sz, a[1] + step)
      q0 = (plo[0] + a[axis], plo[1] + a[1-axis])
      p1 = (plo[0] + b[axis], plo[1] + b[1-axis])
      mv = move.make(q0, p1, mp_jump)
      move.set_name(mv, "J" + "ba"[side] + str(itr))
      (JMS0, JMS1)[side].append(mv)
  return JMS0, JMS1
  # ----------------------------------------------------------------------

def misc_A(mp_trace1, mp_trace2, mp_jump):

  p1 = (1,1)
  p2 = (3,1)
  p3 = (2,3)
  p4 = (4,4)
  p5 = (1,4)

  tr12 = move.make(p1,p2, mp_trace1); move.set_name(tr12, "T0")
  tr34 = move.make(p3,p4, mp_trace2); move.set_name(tr34, "T1")
  tr55 = move.make(p5,p5, mp_trace1); move.set_name(tr55, "T2")
  
  jm23 = move.make(p2,p3, mp_jump); move.set_name(jm23, "J0")

  MVS = [ tr12, jm23, tr34, tr55, ]

  return MVS
  # ----------------------------------------------------------------------
  
def misc_B(mp_trace, mp_jump):

  pa0 = (1,1)
  pa1 = (4,1)
  pa2 = (5,2)
  pa3 = (2,2)
    
  tra01 = move.make(pa0, pa1, mp_trace); move.set_name(tra01, "Ta0")
  tra12 = move.make(pa1, pa2, mp_trace); move.set_name(tra12, "Ta1")
  tra23 = move.make(pa2, pa3, mp_trace); move.set_name(tra23, "Ta2")
  tra30 = move.make(pa3, pa0, mp_trace); move.set_name(tra30, "Ta3")

  pb0 = (6,6)
  pb1 = (5,3)
  trb01 = move.make(pb0, pb1, mp_trace); move.set_name(trb01, "Tb0")
  
  pc0 = (3,6)
  pc1 = (4,3)
  trc01 = move.make(pc0, pc1, mp_trace); move.set_name(trc01, "Tc0")
  
  jmb1c0 = move.make(pb1, pc0, mp_jump); move.set_name(jmb1c0, "Jc0")
  jmc1b0 = move.make(pc1, pb0, mp_jump); move.set_name(jmc1b0, "Jc1")

  pd0 = (0,3)
  pd1 = (1,6)
  pd2 = (2,3)
  trd01 = move.make(pd0, pd1, mp_trace); move.set_name(trd01, "Td0")
  trd12 = move.make(pd1, pd2, mp_trace); move.set_name(trd12, "Td1")
  
  TRS = [ tra01, tra12, tra23, tra30, trb01, trc01, trd01, trd12, ]
  JMS = [ jmb1c0, jmc1b0, ]
  
  return TRS, JMS
  # ----------------------------------------------------------------------

def misc_C(mp_trace, mp_jump):
  
  p00 = (1,3)
  p01 = (3,3)

  p10 = (2,4)
  p11 = (5,4)

  p20 = (4,3)
  p21 = (5,1)

  dp3 = (sqrt(5)/2, 0)
  p30 = tuple(rn.add(p20, dp3))
  p31 = tuple(rn.add(p21, dp3))

  p40 = (5,5)
  p41 = (6,5)

  p50 = (3,6)
  p51 = (4,6)

  p60 = (3,2)
  p61 = p60
  
  p70 = (2,7)
  p71 = (4,7)
 
  p80 = (5,7)
  p81 = (7,7)
  
  p90 = (3,8)
  p91 = (6,8)

  if mp_trace == None:
    TRS = None
  else:
    mv0 = move.make(p00, p01, mp_trace); move.set_name(mv0, "Ta0")
    mv1 = move.make(p10, p11, mp_trace); move.set_name(mv1, "Tc0")
    mv2 = move.make(p20, p21, mp_trace); move.set_name(mv2, "Tb0")
    mv3 = move.make(p30, p31, mp_trace); move.set_name(mv3, "Tc1")
    mv4 = move.make(p40, p41, mp_trace); move.set_name(mv4, "Tb2")
    mv5 = move.make(p50, p51, mp_trace); move.set_name(mv5, "Tb1")
    mv6 = move.make(p60, p61, mp_trace); move.set_name(mv6, "Ta1") # Zero-length move.
    mv7 = move.make(p70, p71, mp_trace); move.set_name(mv7, "Td0")  
    mv8 = move.make(p80, p81, mp_trace); move.set_name(mv8, "Td1")  
    mv9 = move.make(p90, p91, mp_trace); move.set_name(mv9, "Te0")  
    TRS = [ mv0, mv1, mv2, mv3, mv4, mv5, mv6, mv7, mv8, mv9, ]

  if mp_jump == None:
    JMS = None
  else:
    # Make jumps between those moves:
    jm06 = move.make(p01, p60, mp_jump); move.set_name(jm06, "Ja0")
    jm25 = move.make(p20, p50, mp_jump); move.set_name(jm25, "Jb0")
    jm54 = move.make(p51, p40, mp_jump); move.set_name(jm54, "Jd0")
    jm13 = move.make(p11, p30, mp_jump); move.set_name(jm13, "Jc0")
    jm78 = move.make(p71, p80, mp_jump); move.set_name(jm78, "Jb1")
    JMS = [ jm06, jm25, jm54, jm13, jm78, ]

  return TRS, JMS
  # ----------------------------------------------------------------------
