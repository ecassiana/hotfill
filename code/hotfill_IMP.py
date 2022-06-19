# Implementation of the module {hotfill}.

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

def solve(OPHS, mp_jump, maxband, quick):
 
  write_time_plot = False # True to write a file with times of {bandpath.build}.

  nph = len(OPHS) # Number of raster elements.
  assert nph >= 1
  
  # Separate the rasters by scan-line:
  xdir = (1,0)
  ydir = (0,1)
  ystep,yphase = raster.get_spacing_and_phase(OPHS, xdir, ydir)
  OPHS = raster.sort_by_scanline(OPHS, xdir, ydir, ystep, yphase)
  SCS = raster.separate_by_scanline(OPHS, xdir, ydir, ystep, yphase)
  nsc = len(SCS) # Number of scan-lines.
  
  sys.stderr.write("{hotfill.solve}: found %d rasters in %d scanlines\n" % (nph,nsc))
  
  # Allocate the tableaus:
  BTCV = [ None, None]
  TUK = [ None, None ] 
  for z in 0,1:
    BTCV[z] = [ None ]*(nsc + 1)
    TUK[z] = [ None ]*(nsc + 1)
    for i in range(nsc+1):
      BTCV[z][i] = [ None ]*(nsc + 1)
      TUK[z][i] = [ None ]*(nsc + 1)
    
  if write_time_plot: wr_times = open("tests/out/bandpath_times.txt", "w")  # File for bandpath times.

  # Fill the tableaus. The order is not critical as long as 
  # every entry {[u][k][i]} is computed before any entry {[z][i][j]}.
  sys.stderr.write("max band width = %d scanlines\n" % maxband)
  for j in range(1,nsc+1):
    imin = 0 if j <= maxband else j - maxband
    for i in range(imin, j):
      for z in 0,1:
        if write_time_plot: tstart = time.clock_gettime(time.CLOCK_THREAD_CPUTIME_ID)
        
        BTCV[z][i][j] = bandpath.build(OPHS, SCS[i:j], z, mp_jump, quick)

        if write_time_plot:
          tcpu = time.clock_gettime(time.CLOCK_THREAD_CPUTIME_ID) - tstart
          nph_band = 0
          for SS in SCS[i:j]: nph_band += len(SS)
          nsc_band = j - i
          wr_times.write("%5d %5d %10.6f\n" % (nph_band, nsc_band, tcpu))

        TUK[z][i][j] = min_fullpath(TUK, BTCV, z, i, j, mp_jump, maxband)
          
  if write_time_plot: wr_times.close()
  
  # Find the _best complete fullpath:
  z_best = None
  i_best = None
  T_best = +inf
  for i in range(nsc):
    for z in 0,1:
      TUK_zim = TUK[z][i][nsc]
      if TUK_zim != None:
        T_zim, u_zim, k_zim = TUK_zim
        assert T_zim != None
        if T_zim < T_best:
          z_best = z; i_best = i; T_best = T_zim

  # Extract the _best complete fullpath {fph_best}:
  j_best = nsc
  if i_best == None:
    sys.stderr.write("{hotfill.solve} failed\n")
    fph_best = None
    BPHS_best = None
  else:
    sys.stderr.write("{hotfill.solve} succeeded\n");
    sys.stderr.write("  fabtime = %.3f s\n" % T_best)
    fph_best, BPHS_best = recover_fullpath(BTCV, TUK, z_best, i_best, j_best, mp_jump)
    assert abs(path.fabtime(fph_best) - T_best) < 1.0e-8
    # sys.stderr.write("  Rcoolmax = %.3\n", contact.rcoolmax(fph_best,)

  return fph_best, z_best, i_best, j_best, BPHS_best, BTCV, TUK
  # ----------------------------------------------------------------------

def recover_fullpath(BTCV, TUK, z, i, j, mp_jump):
  
  bph_zij, CTS_lo_zij, TCVS_lo_zij, CTS_hi_zij, TCVS_hi_zij = BTCV[z][i][j]
  
  if i == 0:
    fph_zij = bph_zij
    BPHS_zij = [bph_zij,]
  else:
    T_zij, u_zij, k_zij = TUK[z][i][j]
    fph_uki, BPHS_uki = recover_fullpath(BTCV, TUK, u_zij, k_zij, i, mp_jump)
    use_links = True
    fph_zij = path.concat([fph_uki,bph_zij,], use_links, mp_jump)
    cph_i = path.get_connector(fph_uki, bph_zij, use_links, mp_jump)
    BPHS_zij = BPHS_uki + [cph_i, bph_zij,]
  sys.stderr.write("  band [%d][%03d][%03d] (%3d scanlines)\n" % (z,i,j,j-i))
  return fph_zij, BPHS_zij
  # ----------------------------------------------------------------------

def min_fullpath(TUK, BTCV, z, i, j, mp_jump, maxband):
  # On input, {TUK} and {BTCV} must be the {hotfill} tableaus. The
  # entries {BTCV[u][k][i]} and {TUK[u][k][i]} must have been computed
  # for all {u} in {0,1} and all {k} in {0..i-1}.
  #
  # In addition, {BTCV[z][i][j]}, if not {(None,None,None)}, must be {(bph,TCVS_lo,TCVS_hi)}
  # where {bph} is a {z}-type bandpath for the band {(i,j)}, and {TCVS_lo,TCVS_hi}
  # are its cover times for the contacts on cut-lines {i} and {j}, as
  # described in {hotfill.solve}.
  #
  # The procedure returns a triple {(T,u,k)}, where {T} is the fabtime
  # of the best valid {(i,j)}-fullpath {fph[z][i][j]} that ends with the
  # bandpath {BTCV[z][i][j]}, and {u,k} are such that {fph[z][i][j]} is
  # the concatenation of {fph[u][k][i]} and {bph[z][i][j]}.
  #
  # If {i == 0}, the path {fph[z][i][j]} is just the bandpath {bph[z][i][j]}.
  # Then {T} is {path.fabtime(bph[z][i][j])}, and {u,k} are {None}.
  #
  # If there is no such path (in particular, if {bph[z][i][j]} is {None})
  # returns {(None,None,None)}.
  #
  
  bph_zij, CTS_lo_zij, TCVS_lo_zij, CTS_hi_zij, TCVS_hi_zij = BTCV[z][i][j]
      
  if bph_zij == None:
    TUK_zij = (+inf, None, None)
  elif i == 0:
    TUK_zij = (path.fabtime(bph_zij), None, None)
  else:
    # Try all fullpaths that end with the bandpath {bph_zij}:
    TB_zij = path.fabtime(bph_zij)
    u_best = None
    k_best = None
    TF_best = +inf
    for u in 0,1:
      kmin = max(0,i-maxband)  
      for k in range(kmin,i):
        TUK_uki = TUK[u][k][i]
        if TUK_uki != None:
          TF_uki, uu_uki, kk_uki = TUK_uki
          if TF_uki != +inf:
            bph_uki, CTS_lo_uki, TCVS_lo_uki, CTS_hi_uki, TCVS_hi_uki = BTCV[u][k][i]
            use_links = True
            TC = path.connection_time(bph_uki, bph_zij, use_links, mp_jump)
            if check_valid(CTS_hi_uki, TCVS_hi_uki, TC, CTS_lo_zij, TCVS_lo_zij):
              TF_cand = TF_uki + TC + TB_zij
              if TF_cand < TF_best:
                u_best = u; k_best = k; TF_best = TF_cand
    TUK_zij = (TF_best, u_best, k_best)
  return TUK_zij
  # ----------------------------------------------------------------------
  
def check_valid(CTS0, TCVS0, TC, CTS1, TCVS1):
  # Returns {True} iff the concatenation of a {u,k,i}-fullpath {fph} and
  # a {z,i,j}-bandpath {bph} would be valid.
  #
  # Assumes that internal contacts in the two paths are valid, so it
  # needs to check only the validity of the contacts between them;
  # namely, the contacts on the cut-line {i}.
  #
  # The lists {CTS0} and {CTS1} must have equal, and hold the contacts
  # on cut-line {i}, in left-to-right order. 
  #
  # The list {TCVS0} must have the cover times of those contacts by the
  # path {rev(fph)}. The list {TCVS1} have the cover times of those
  # contacts by the path {bph}.
  #
  # The parameter {TC} must be the fabtime of the connector (link path
  # or jump) between {fph} and {bph}.
  
  nct = len(CTS0)
  assert len(CTS1)  == nct
  assert len(TCVS0) == nct
  assert len(TCVS0) == nct
  
  for ict in range(nct):
    ct = CTS0[ict]
    assert CTS1[ict] == ct
    tcv0 = TCVS0[ict]
    tcv1 = TCVS1[ict]
    tcool = tcv0 + TC + tcv1
    if tcool > contact.tcool_limit(ct): return False
  return True
  # ----------------------------------------------------------------------
 
def describe_and_check_input(wr, OPHS, mp_jump, maxband, quick):
  wr.write("### INPUT DATA ######################################################################\n")
  wr.write("maxband = %d\n" % maxband)
  wr.write("quick = %s\n" % quick)
  wr.write("!! Not implemented !!\n")
  wr.write("#####################################################################################\n")
  return
  # ----------------------------------------------------------------------

