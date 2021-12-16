#! /usr/bin/python3
# Test program for the module {rootray}
# Last edited on 2021-10-03 18:35:45 by stolfi

import rootray
from math import inf
import sys

fa = (+1, [ 0.234, 1.000, 1.0000001, 5,  10,   100000000, 100000002 ])
sys.stderr.write("fa = %s\n" % str(fa))

fb = (-1, [ 0.3, 0.4, 0.5, 5, 6, 7, 8, 9, 10,    100000002, 100000004 ])
sys.stderr.write("fb = %s\n" % str(fb))

def test_minstep():
  sys.stderr.write("--- testing {rootray.minsep} ---------------\n")

  for x1,x2 in (
      ( 1.000, 1.000 + 1.0e-7 ), 
      ( 1.0e8, 1.0e8 + 2 ),
    ):
      sep = rootray.min_sep(x1,x2)
      sys.stderr.write("minsep(%21.15e, %21.15e) = %21.15e\n" % (x1,x2,sep))
  return
  # ----------------------------------------------------------------------

def test_validate():
  sys.stderr.write("--- testing {rootray.validate} ---------------\n")

  rootray.validate(fa)
  rootray.validate(fb)
  # ----------------------------------------------------------------------

def test_complement():
  sys.stderr.write("--- testing {rootray.complement} ---------------\n")

  fc1 = rootray.complement(fa)
  sys.stderr.write("fa  = %s\n" % str(fa))
  sys.stderr.write("fc1 = %s\n" % str(fc1))
  rootray.validate(fc1)
  assert fc1[0] == -1
  assert fc1[1] == fa[1]
  return
  # ----------------------------------------------------------------------

def test_union():
  sys.stderr.write("--- testing {rootray.union} ---------------\n")

  fua =  (+1, [0.2564, 0.6666])
  fub =  (+1, [0.6923, 0.7435])

  fu1 = rootray.union(fua, fub)
  sys.stderr.write("fu1 = %s\n" % str(fu1))
  rootray.validate(fu1)
  assert fu1 == (+1, [0.2564, 0.6666, 0.6923, 0.7435])

  fu2 = rootray.union(fua, fua)
  sys.stderr.write("fu2 = %s\n" % str(fu2))
  rootray.validate(fu2)
  assert fu2 == (+1, [0.2564, 0.6666])

  fu3 = rootray.union(fua, rootray.complement(fua))
  rootray.validate(fu3)
  assert fu3 == (-1, [])

  fu4 = rootray.union(fua, rootray.vacuum())
  rootray.validate(fu4)
  assert fu4 == fua

  fu5 = rootray.union(fub, rootray.plenum())
  rootray.validate(fu5)
  assert fu5 == (-1, [])

  fuc = (+1, [-2, +2])
  fud = (+1, [-1, +3])

  fu8 = rootray.union(fuc, fud)
  sys.stderr.write("fu8 = %s\n" % str(fu8))
  assert fu8 == (+1, [-2, +3 ])
  return
  # ----------------------------------------------------------------------

def test_intersection():
  sys.stderr.write("--- testing {rootray.intersection} ---------------\n")

  fia =  (-1, [0.2564, 0.6666])
  fib =  (-1, [0.6923, 0.7435])

  fi1 = rootray.intersection(fia, fib)
  sys.stderr.write("fi1 = %s\n" % str(fi1))
  rootray.validate(fi1)
  assert fi1 == (-1, [0.2564, 0.6666, 0.6923, 0.7435])

  fi2 = rootray.intersection(fib, fib)
  rootray.validate(fi2)
  assert fi2 == fib

  fi3 = rootray.intersection(fia, rootray.complement(fia))
  rootray.validate(fi3)
  assert fi3 == (+1, [])

  fi4 = rootray.intersection(fia, rootray.plenum())
  rootray.validate(fi4)
  assert fi4 == fia

  fi5 = rootray.intersection(fib, rootray.vacuum())
  rootray.validate(fi5)
  assert fi5 == (+1, [])

  fic = (+1, [-2, +2])
  fid = (+1, [-1, +3])
  
  fi8 = rootray.intersection(fic, fid)
  sys.stderr.write("fi8 = %s\n" % str(fi8))
  assert fi8 == (+1, [-1, +2 ])
  return
  # ----------------------------------------------------------------------

def test_difference():

  sys.stderr.write("--- testing {rootray.difference} ---------------\n")

  fda =  (-1, [0.2564, 0.6666])
  fdb =  (+1, [0.6923, 0.7435])

  fd1 = rootray.difference(fda, fdb)
  sys.stderr.write("fd1 = %s\n" % str(fd1))
  rootray.validate(fd1)
  assert fd1 == (-1, [0.2564, 0.6666, 0.6923, 0.7435])

  fd2 = rootray.difference(fdb, fdb)
  rootray.validate(fd2)
  assert fd2 == (+1, [])

  fd3 = rootray.difference(fda, rootray.complement(fda))
  rootray.validate(fd3)
  assert fd3 == fda 

  fd4 = rootray.difference(fda, rootray.vacuum())
  rootray.validate(fd4)
  assert fd4 == fda

  fd5 = rootray.difference(fdb, rootray.plenum())
  rootray.validate(fd5)
  assert fd5 == (+1, [])

  fdc = (+1, [-2, +2])
  fdd = (+1, [-1, +3])
  
  fd8 = rootray.difference(fdc, fdd)
  sys.stderr.write("fd8 = %s\n" % str(fd8))
  assert fd8 == (+1, [-2, -1 ])
  return
  # ----------------------------------------------------------------------

def test_clip_times():
  sys.stderr.write("--- testing {rootray.clip_times} ---------------\n")

  fxa =  (-1, [0.2564, 0.6666])
  fxb =  (+1, [0.6923, 0.7435])

  fx1 = rootray.clip_time(fxa, 00.00, +0.68)
  rootray.validate(fx1)
  assert fx1 == fxa

  fx2 = rootray.clip_time(fxa, -1.00, +0.50)
  rootray.validate(fx2)
  assert fx2 == (-1, [0.2564,] )

  fx3 = rootray.clip_time(fxa, +0.50, +1.00)
  rootray.validate(fx3)
  assert fx3 == (+1, [0.6666,])

  fx4 = rootray.clip_time(fxa, 1.0e10, 2.0e10)
  rootray.validate(fx4)
  assert fx4 == (-1, [])

  fx6 = rootray.clip_time(fxb, 00.00, +0.68)
  rootray.validate(fx6)
  assert fx6 == (+1, [])
  
  return
  # ----------------------------------------------------------------------
  
test_minstep()
test_validate()
test_complement()
test_union()
test_intersection()
test_difference()
test_clip_times()
