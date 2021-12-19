# Main functions of the HotFill algorithm.
# Last edited on 2021-11-25 01:16:15 by stolfi

import hotfill_IMP
  
def solve(OPHS, mp_jump, maxband, quick):
  # This procedure tries to build a tool-path for a solid fill from a
  # list {OPHS} of raster elements, respecting the cooling time limits
  # of the contacts attached to them, while trying to minimize the total
  # fabrication time.
  #
  # The results will be:
  #
  #  * {fph}, a single path that includes every trace of {OPH}
  #    or its reverse, exactly once; plus all applicable links;
  #    plus zero or more jumps, with {Move_Parms} parameters {mp_jump}.
  #
  #  * The type {z_best} and the cutline indices {i_best,j_best}
  #    of the last bandpath in {fph}.
  #
  #  * The list {BPHS} of the bandpaths that were concatenated
  #    to make the resulting fullpath, intercalated with the connectors
  #    that were inserted between them. Each connector will be either
  #    an oriented path or a {Move} object which is a jump.
  #
  #  * The tableaus {BTCV} and {TUK} of the dynamic programming.
  #
  # The parameter {OPHS} should be a list of oriented paths, each
  # consisting of a single horizontal trace and oriented left to right.
  # These raster elements should lie on a number of uniformly spaced
  # horizontal scan-lines. All raster traces are assumed to have the
  # same {MoveParms} record, and the nominal width is assumedto be equal
  # to the scanline spacing.
  # 
  # Each raster element {oph} in {OPHS} should have a set
  # {path.get_links(oph)} with zero or more link paths that end at
  # {pini(oph)}; and a set {path.get_links(rev(oph))} with zero or more
  # link paths that end at {pfin(oph)}. The reverse of each of those
  # link paths should be a link path of some other raster trace, usually
  # on a scan-line adjacent to that o {oph}.
  #
  # The procedure {path.get_contacts(oph,isd)} should return the set of
  # all relevant contacts whose side {isd} is the trace of the raster
  # element {oph}. Side 0 of each contact must be the lowest of the two
  # traces. The cooling time limits of all those contacts should have
  # been set to a value greater than 0 (possibly {+oo}).
  #
  # The resulting path will be a concatenation of /bandpaths/, 
  # where each bandpath uses all the rasters in a set of contuguous
  # scan-lines.  The procedure uses dynamic programming to find the
  # partition of the scan-lines into bands that provides a 
  # tool-path that is valid (respectd all contact cooling constraints)
  # with hopefully small fabtime.
  #
  # The bandpaths will be limited to {maxband} scanlines.
  # 
  # If {quick} is false, uses greedy bandpath heuristic, whose running
  # time is quadratic on the number of rasters in the band. If {quick}
  # is true, uses another heuristic with has linear running time, but
  # may (or may not) have bigger fabtime.
  return hotfill_IMP.solve(OPHS, mp_jump, maxband, quick)

def describe_and_check_input(wr, OPHS, mp_jump, maxband, quick):
  # Prints main attributes of {solve}'s input data to 
  # file {wr} in human-readable form.
  return hotfill_IMP.describe_and_check_input(wr, OPHS, mp_jump, maxband, quick)

# THE DYNAMIC PROGRAMING TABLEAUX

# The procedure uses two tables {BTCV} and {TUK} with indices {[z][i][j]}
# where {z} is 0 or 1, and {i,j} are in {0..nsc}.
# 
# The entry {BTCV[z][i][j]} will always be {None} if {i >= j}.
# Otherwise, once it has been computed it will be a quintuple
# {(bph,CTS_lo,TCVS_lo,CTS_hi,TCVS_hi)} where:
# 
#   * {bph} is the bandpath {bph[z][i][j]}.
#
#   * {CTS_lo} is a list of the contacts on cut-line {i}.
#
#   * {TCVS_lo} is a list of {tcov(bph,ct)} for those contacts.
#
#   * {CTS_hi} is a list of the contacts on cut-line {j}.
#
#   * {TCVS_hi} is a list of {tcov(rev(bph),ct)} for those contacts.
# 
# The lists {CTS_lo,CTS_hi} are sorted left-to-right, and the lists {TCVS_lo,TCVS_hi}
# are in the same order.
# 
# The entry {BTCV[z][i][j]} may be {(None,None,None,None,None)} to mean that the
# proceduree could not or would not build a valid bandpath of type {z}
# for the band {(i,j)}.
#
# The entry {TUK[z][i][j]} will always be {None} if {i >= j}.  Othwerwise,
# once it has been computed it will be a triple {(T,u,k)} where:
#
#   * {T} is the fabtime of the fullpath {fph[z][i][j]}.
# 
#   * {u} and {k} are such that {fph[z][i][j]} is {fph[u][k][i]} concatenated with {bph[z][i][j]}.
#
# The entry {TUK[z][i][j]} may be {(None,None,None)} to mean that the
# procedue could not build a valid fullpath ending with a {z}-type
# bandpath for the band {(i,j)}. In particular, if {BTCV[z][i][j]} is
# {(None,None,None,None,None)}, then {TUK[z][i][j]} is {(None,None,None)} .

def recover_fullpath(BTCV, TUK, z, i, j, mp_jump):
  # Obtains the best {z,i,j}-fullpath from the tableaus {BTCV} and {TUK}.
  # Also returns a list {BPHS} of the bandpaths,
  # links, and jumps that compose that fullpath.
  return hotfill_IMP.recover_fullpath(BTCV, TUK, z, i, j, mp_jump)
