# Procedure to create examples of moves.

import move_example_IMP
import move_parms
import pyx

def rectangle_rasters(plo, axis, n, sz, step, mp_trace):
  # Returns a list {TRS} with {n} raster line traces, all with the
  # parameters record {mp_trace}. If {axis} is 0 the raster traces will
  # be horizontal, oriented from left to right and listed from bottom to
  # top. If {axis} is 1 the raster traces will be vertical, oriented
  # from bottom to top and listed from left to right.
  #
  # In any case, the raster traces will be {sz} mm long and spaced
  # {step} mm apart, and the first raster will start at {plo}. Their
  # endpoints will fit a rectangle with corners {plo} and {plo +
  # (sz,(n-1)*step)} or {plo + ((n-1)*step, sz)}, depending on {axis}.
  return move_example_IMP.rectangle_rasters(plo, axis, n, sz, step, mp_trace)
  
def rectangle_links(plo, axis, n, sz, step, mp_trace):
  # Returns two lists {LKS0,LKS1} of traces, which are the links needed
  # to connect the traces of {rectangle_rasters} or their reversals into
  # a single continuous path.
  # 
  # Specifically, The list {LKS0} contains a link trace connecting the
  # initial point of each raster (except the last one) to the initial
  # point of the next raster. The list {LKS1} is similar but connects
  # the final points instead. All links will have the same parameter
  # record {mp_trace}.
  return move_example_IMP.rectangle_links(plo, axis, n, sz, step, mp_trace)
  
def rectangle_jumps(plo, axis, n, sz, step, mp_jump):
  # Returns two lists {JMS0,JMS1} of jumps, which are the links needed
  # to connect the traces of {rectangle_rasters} or their reversals into
  # a single path.
  # 
  # Specifically, the list {JMS0} contains a jump that connects the final point of
  # each raster (except the last one) to the intial point of the next
  # one. The list {JMS1} is similar but connects the initial point of
  # each raster to the final point of the next one.
  return move_example_IMP.rectangle_jumps(plo, axis, n, sz, step, mp_jump)

def misc_A(mp_trace1, mp_trace2, mp_jump):
  # Returns a list {MVS} of four {Move} objects.
  # Includes two traces with parameters {mp_trace1} and {mp_trace2},
  # one jump with parameters {mp_jump} connecting them,
  # and one isolated zero-length trace with parameters {mp_trace1}.
  #
  #  Move   Type    pini  pfin  obs
  #  ------ -----   ----- ----- ------------------
  #  MVS[0] trace1  (1,1) (3,1)
  #  MVS[1] jump    (3,1) (2,3) continues {MVS[0]}
  #  MVS[2] trace2  (2,3) (4,4) continues {MVS[1]}
  #  MVS[3] trace1  (1,4) (1,4) isolated dot
  #
  return move_example_IMP.misc_A(mp_trace1, mp_trace2, mp_jump)
  
def misc_B(mp_trace, mp_jump):
  # Returns a list {TRS} of eight {Move} objects with parameters
  # {mp_trace}, and a list {JMS} of two {Move} objects with parameters
  # {mp_jump}. If either parameter is {None}, the corresponding result
  # will be {None}.
  #
  # If the nominal width of mp_trace is 1.0 mm, several pairs of
  # traces will be adjacent.
  #
  #  Move   Type   pini  pfin  obs
  #  ------ -----  ----- ----- ------------------
  #  TRS[0] trace  (1,1) (4,1) horizontal
  #  TRS[1] trace  (4,1) (5,2) diagonal, links {TRS[0]} to {TRS[2]}
  #  TRS[2] trace  (5,2) (2,2) horizontal adj to {TRS[0]}
  #  TRS[3] trace  (2,2) (1,1) diagonal, links {TRS[2]} to {TRS[0]}
  #  TRS[4] trace  (6,6) (5,3) diagonal
  #  TRS[5] trace  (3,6) (4,3) diagonal
  #  TRS[6] trace  (0,3) (1,6) diagonal
  #  TRS[7] trace  (1,6) (2,3) diagonal, extends {TRS[6]}
  #
  #  JMS[0] jump   (5,3) (3,6) diagonal, joins {TRS[4]} to {TRS[5]}
  #  JMS[1] jump   (4,3) (6,6) diagonal, joins {TRS[5]} to {TRS[4]} 
  return move_example_IMP.misc_B(mp_trace, mp_jump)

def misc_C(mp_trace, mp_jump):
  # Returns a list {TRS} of ten {Move} objects with parameters
  # {mp_trace}, and a list {JMS} of five {Move} objects with parameters
  # {mp_jump}. If either parameter is {None}, the corresponding result
  # will be {None}.
  #
  # If the nominal width of mp_trace is 1.0 mm, several pairs of
  # traces will be adjacent.
  #
  # Traces {TRS[k]}:
  #
  #   k name  pini  pfin obs
  #   - ---- ----- ----- ---------------
  #   0  Ta0 (1,3) (3,3) horizontal                             
  #   1  Tc0 (2,4) (5,4) horizontal, adj to {TRS[0]}            
  #   2  Tb0 (4,3) (5,1) diagonal                               
  #   3  Tc1 (*,3) (@,1) diagonal, adj to {TRS[2]}              
  #   4  Tb2 (5,5) (6,5) horizontal                             
  #   5  Tb1 (3,6) (4,6) horizontal                             
  #   6  Ta1 (3,2) (3,2) isolated dot                           
  #   7  Td0 (2,7) (4,7) horizontal                             
  #   8  Td1 (5,7) (7,7) horizontal, collinear with {TRS[7]}    
  #   9  Te0 (3,8) (6,8) horizontal adj to {TRS[7]} and {TRS[8]}
  #
  # Jumps {JMS[k]}:
  #
  #   k name  pini  pfin obs
  #   - ---- ----- ----- ---------------
  #   0  Ja0 (3,3) (3,2) joins {TRS[0]} to {TRS[6]}, vertical  
  #   1  Jb0 (4,3) (3,6) joins {rev(TRS[2])} to {TRS[5]}       
  #   2  Jd0 (4,6) (5,5) joins {TRS[5]} to {TRS[4]}            
  #   3  Jc0 (5,4) (*,3) joins {TRS[1]} to {TRS[3]}            
  #   4  Jb1 (4,7) (5,7) joins {TRS[7]} to {TRS[8]}, horizontal
  # 
  # The coordinates "*" and "@" are {4+a} and {5+a} where 
  # {a} is such that the axes of the traces {TRS[2]} and {TRS[3]}
  # are 1 mm apart.
  return move_example_IMP.misc_C(mp_trace, mp_jump)
