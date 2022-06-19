# Tools for ray-tracing csg models in arbitrary spaces with arbitrary probing trajectories.

import rootray_IMP

# A /root-ray/ is a pair {(s,T)} where {s} is a sign {+1} or {-1}, and
# {T} is a list of zero or more /root times/ (floats), that describes
# the state of a point {r(t)} that travels along some continuous trajectory
# (usually a straight line) with respect to some figure {F} in some space.
#
# The sign {s} of the root-ray tells the position of {r(t)} relative to
# {F} as {t} tends to {-oo}: {+1} means in {F}'s exterior and {-1}
# means in {F}'s interior. Each root time in {T} is a value of {t} where
# {r(t)} crosses from the interior of {F} to the exterior, or
# vice-versa.  The root times in {T} are always strictly increasing.
#
# This module assumes that {r(t)} cannot be on the boundary of {F} for
# more than a discrete set of values of {t}. It also assumes that 
# {r(t)} cannot touch the boundary without crossing it.  This will be 
# referred to as the /transversality/ condition. 
#
# These procedures in this module compute the root-ray {f} of a
# trajectory {r} relative to a figure {F}, which is defined by a Boolean
# operation on zero or more figures {F1,F2,...}, from the root-rays
# {f1,f2,...} of {r} relative to those figures.
#
# Whenever a root-ray procedure, as defined in the documentation string,
# would violate this condition, it will implicitly perturb the
# trajactory by a tiny amount so that the contact with the boundary is
# eliminated or reduced to a finite number of discrete isolated
# crossings.
#
# In fact, since the crossing times are computed with floating-point
# arithmetic, two crossing times that are too close together are assumed
# to be due to roundoff errors in the data or in the computation, and
# eliminated from the root-ray. As a result, extremely thin parts or
# gaps in the figure may disappear in the the root-ray, in an inconsistent
# and unpredictable way.
  
def vacuum():
  # Returns a root-ray {(+1,[])} that describes a trajectory always in the exterior
  # of the figure.  It is the neutral element for {union} and the absorbing element for 
  # {intersection}.
  return rootray_IMP.vacuum()

def plenum():
  # Returns a root-ray {(-1,[])} that describes a trajectory always in the interior
  # of the figure.  It is the absorbing element for {union} and the neutral element for 
  # {intersection}.
  return rootray_IMP.plenum()
  
def complement(f):
  # Returns a root-ray {g} for the complement {G = \RR^3\setminus F} of a figure {F},
  # given the root-ray {f} of the latter.  Namely, if {f=(s,T)}, returns {(-s,T)}.
  return rootray_IMP.complement(f)
  
def union(f1, f2):
  # Given the root-rays {f1,f2} of a trajectory {r} relative to figures {F1,F2},
  # returns the root-ray of {r} relative to the union {F1 \cup F2}.
  # Note that two parts of {F1} and {F2} which are separated by a very thin gap
  # may or may not become welded together.
  return rootray_IMP.union(f1, f2)
  
def intersection(f1, f2):
  # Given the root-rays {f1,f2} of a trajectory {r} relative to figures {F1,F2},
  # returns the root-ray of {r} relative to the intersection of {F1 \cap F2}.
  # Note that a very thin layer that is created by barely overlapping
  # parts of {F1} and {F2} may or may not disapear.
  return rootray_IMP.intersection(f1, f2)

def difference(f1, f2):
  # Given the root-rays {f1,f2} of a trajectory {r} relative to figures {F1,F2},
  # returns the root-ray of {r} relative to the set difference of {F1 \setminus F2}.
  # It is equivalent to {intersection(f1,complement(f2))}. Note that 
  # a part of {F2} that is inside {F1} but very close to its boundary
  # may or may not punch a hole there.
  return rootray_IMP.difference(f1, f2)

def min_sep(t1,t2):
  # Returns the minimum separation between finite root times {t1} and {t2}.
  # It is monotonic in {abs(t1)} and {abs(t2)}.
  return rootray_IMP.min_sep(t1,t2)

def clip_time(f, tmin, tmax):
  # Given a root-ray {f} for some trajectory {r(t)}, returns only the section where
  # the parameter {t} is limited to the range {(tmin _ tmax)}.   Requires {tmin < tmax}. 
  #
  # Namely, if {f=(s,T)}, returns {(s',T')} where {T'} contains only the
  # root times in that range, and {s'} is {s} or {-s} depending on whether
  # an even or odd number of times was removed at the front of {T}.
  #
  # This operation is equivalent to computing {op1(op2(f,gmin),gmax)} where 
  # {op1} and {op2} are {intersection} or {union}, and the root-rays
  # {gmin} and {gmax} have a single root time each, respectively 
  # at {tmin} and {tmax}, with suitable signs.
  return rootray_IMP.clip_time(f, tmin, tmax)

def validate(f):
  # Runs consistency tests on {f}, including minimum root separation.
  return rootray_IMP.validate(f)
