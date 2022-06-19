# Implementation of module {rootray}.

import rn
from math import hypot, inf
import sys

def min_sep(t1,t2):
  e1 = 1.0e-8*t1
  e2 = 1.0e-8*t2
  eps = 1.0e-100 + hypot(e1,e2)
  assert eps != +inf, "not applicable to infinities"
  return eps

def vacuum():
  return (+1,[])
  
def plenum():
  return (-1,[])
  
def complement(f):
  assert (type(f) is tuple or type(f) is list) and len(f) == 2
  s, T = f
  return (-s, T)
  
def union(f1, f2):
  # sys.stderr.write("\n")
  # sys.stderr.write("--- enter union ---\n")
  # sys.stderr.write("f1 = %s f2 = %s\n" % (str(f1),str(f2)))
  assert (type(f1) is tuple or type(f1) is list) and len(f1) == 2
  s1, T1 = f1
  assert (type(f2) is tuple or type(f2) is list) and len(f2) == 2
  s2, T2 = f2
  
  # Merge root times removing buried ones:
  n1 = len(T1); k1 = 0; t1 = T1[0] if n1 >0  else +inf
  n2 = len(T2); k2 = 0; t2 = T2[0] if n2 >0  else +inf
  
  s = -1 if s1 < 0 or s2 < 0 else +1
  T = []
  tant = -inf
  scur = s
  while t1 < +inf or t2 < +inf:
    # At this point we have created in {T} the root-ray up to time {tcur=tant+eps} 
    # and the status of {r(tcur)} relative to the union is given by {scur}.
    # Also the status of {r1(tcur)} and {r2(tcur)} is {s1} and {s2}, respectively.
    # The next root times of {f1} and {f2} after {tant} are {t1} and {t2}
    # wich are {T1[k1]} and {T2[k2]} if those elements exist.
    
    # sys.stderr.write("\n")
    # sys.stderr.write("tant = %12.8f s1 = %+2d s2 = %+2d scur = %+2d\n" %(tant,s1,s2,scur))
    assert scur < 0 or s1 > 0 or s2 > 0

    # Identify the next event:
    tnext = min(t1, t2)
    assert tnext > tant
    
    # Identify the status of the two root-rays and of the union after the next event:
    if tnext == t1: s1 = -s1;
    if tnext == t2: s2 = -s2;
    snext = -1 if s1 < 0 or s2 < 0 else +1

    # sys.stderr.write("tnext = %12.8f s1 = %+2d s2 = %+2d snext = %12.8f \n" % (tnext,s1,s2,snext))
    
    # See if union changed:
    if snext != scur:
      # Append {tnext} and flip the status:
      scur = snext
      T.append(tnext)

    # In any case, advance to {tnext+eps}
    tant = tnext
    
    # Advance the two operands:
    if tnext < +inf:
      if t1 == tnext: 
        assert k1 < n1; k1 = k1 + 1; t1 = T1[k1] if k1 < n1 else +inf
      if t2 == tnext:
        assert k2 < n2; k2 = k2 + 1; t2 = T2[k2] if k2 < n2 else +inf
  # sys.stderr.write("union: s = %d T = %s\n" % (s, T))
  return (s, remove_near_dups(T))
  
def remove_near_dups(T):
  # Returns a copy of {T} removing all pairs of 
  # consecutive elements that are not well-separated.
  # sys.stderr.write("remove_near_dups: T = %s\n" % str(T)) 
  n = len(T)
  Tclean = []
  k = 0 # Next element of {T} to consider
  tant = -inf # Last element added to {Tclean}
  while k < n-1:
    # At this point {T[k]} and {T[k+1]} exist, and all elements of {T}
    # before {T[k]} have been added to {Tclean}, and are well separated
    # from {T[k]}.
    ta = T[k]; tb = T[k+1]
    # sys.stderr.write("remove_near_dups: k= %d ta = %21.15e tb = %21.15e\n" % (k, ta, tb)) 
    assert -inf < ta and ta < tb and tb < +inf
    if tb - ta < min_sep(ta, tb):
      # Eliminate {T[k]} and {T[k+1]}:
      k = k + 2
    else:
      # Safe to add {T[k]}:
      if tant != -inf: assert ta - tant > min_sep(tant, ta)
      Tclean.append(ta)
      tant = ta
      k = k + 1
  
  # Add the last unpaired element if any:
  if k < n: 
    ta = T[k]
    if tant != -inf: assert ta - tant > min_sep(tant, ta)
    Tclean.append(ta); k = k + 1
  assert k == n
  return Tclean
  
def intersection(f1, f2):
  # sys.stderr.write("\n")
  # sys.stderr.write("--- enter intersection ---\n")
  # sys.stderr.write("f1 = %s f2 = %s\n" % (str(f1),str(f2)))
  return complement(union(complement(f1), complement(f2)))
 
def difference(f1, f2):
  return complement(union(complement(f1), f2))
  
def clip_time(f, tmin, tmax):
  assert (type(f) is tuple or type(f) is list) and len(f) == 2
  assert tmin < tmax
  s, T = f
  Tnew = []
  tant = -inf
  for t in T:
    assert tant < t and t < +inf
    if t >= tmax:
      break
    elif t > tmin:
      Tnew.append(t)
    else:
      s = -s
    tant = t
  return (s, Tnew)
  
def validate(f):
  assert (type(f) is tuple or type(f) is list) and len(f) == 2
  s,T = f
  assert type(s) is int and (s == -1  or s == +1)
  assert type(T) is list or type(T) is tuple
  tant = -inf
  for t in T:
    assert type(t) is int or type(t) is float
    assert t > tant and t < +inf
    if tant != -inf: assert t - tant > min_sep(tant, t)
