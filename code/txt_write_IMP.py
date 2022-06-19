# Implementation of the module {txt_write}.

import txt_write

import block
import move
import move_parms
import path
import raster
import contact
import hacks
import rn

import sys
from math import sqrt, sin, cos, atan2, log, exp, floor, ceil, inf, nan, pi

def write(wr, OCRS, OPHS, Z, angle, shift):

  ncr = len(OCRS)
  nph = len(OPHS)

  # Simple parameters:
  wr.write("Z%.6f\n" % Z)
  wr.write("N%d\n" % nph)
  wr.write("K%d\n" % ncr)
  
  # Make copy {OPHS_srt} of {OPHS} sorted by scan-line order:
  assert nph >= 1
  ystep, yphase = raster.get_spacing_and_phase(OPHS, (1,0), (0,1))
  OPHS_srt = raster.sort_by_scanline(OPHS, (1,0), (0,1), ystep, yphase)
  
  # Count distinct contacts and link paths:
  CTS = set()
  LKS = set()
  for oph in OPHS:
    for idr in range(2):
      for olki in path.get_links(path.spin(oph, idr)):
        lki, dri = path.unpack(olki)
        LKS.add(lki)
    for isd in range(2):
      for ct in path.get_contacts(oph,isd) :
        CTS.add(ct)
  nct = len(CTS)
  nlk = len(LKS)

  sys.stderr.write("writing %d contours, %d rasters, %d links, %d contacts\n" % (ncr,nph,nlk,nct)) 

  xdir = ( cos(angle), sin(angle) )
  ydir = ( -xdir[1], +xdir[0] )
  write_rasters(wr, OPHS_srt, nph, xdir, ydir, shift)
  if ncr > 0: write_contours(wr, OCRS, ncr, xdir, ydir, shift)
  if nct > 0 or nlk > 0: write_contacts_and_links(wr, OPHS_srt, nct, nlk, xdir, ydir, shift, ystep, yphase)
  wr.flush()
  return
  # ----------------------------------------------------------------------

def encode_point(p, xdir, ydir, shift, sep):
  # Returns a string with the coordinates of point {p}, as two floats
  # separated by the string {sep}, after subtracting the vector {shift} and
  # mapping the result by the linear map that takes {(1,0),(0,1)} to
  # {xdir,ydir}.
  ps = rn.sub(p, shift)
  pr = rn.mix(ps[0], xdir, ps[1], ydir)
  txt = ("%.6f%s%.6f" % (pr[0], sep, pr[1]))
  return txt
  # ----------------------------------------------------------------------

def write_rasters(wr, OPHS, nph, xdir, ydir, shift):
  # Write the 'R' lines for the rasters in {OPHS}, assumed to be in scanline order..
  # As side effect, will reorient the rasters in {OPHS} to point left to right.
  # There should be {nph} such lines.

  assert nph == len(OPHS)
  wr.write("### LINHAS DE RASTER \n")
  sys.stderr.write("  ### RASTERS\n")
  
  nwrite_R = 0
  for iph in range(nph):
    oph = OPHS[iph]
    assert path.nelems(oph) == 1
    p = path.pini(oph); 
    q = path.pfin(oph); 
    # Write the raster oriented left to right:
    xp = rn.dot(p, xdir)
    xq = rn.dot(q, xdir)
    if xp > xq:
      p,q = q,p; xp,xq = xq,xp
      OPHS[iph] = path.rev(OPHS[iph])
      rbit = 1      
    else:
      rbit = 0
    igr = path.get_group(oph)
    assert igr != None and igr >= 0
    wr.write("R%d," % iph); nwrite_R += 1
    wr.write(encode_point(p, xdir, ydir, shift, ","))
    wr.write(",")
    wr.write(encode_point(q, xdir, ydir, shift, ","))
    wr.write(",%d,%d\n" % (rbit, igr))

  assert nwrite_R == nph
  return
  # ----------------------------------------------------------------------

def write_contacts_and_links(wr, OPHS, nct, nlk, xdir, ydir, shift, ystep, yphase):
  # Writes an 'L' line to {wr} for every pair of 
  # rasters in {OPHS} that are in consecutive scanlines
  # and make a contact, as defined by {path.get_contacts}.
  # 
  # Also writes in that line the geometry of the link paths
  # between those two rasters, as defined by {path.get_links}.
  # Assumes that {OPHS} is in scanline order and the
  # scanline positions along {ydir} have separation {ystep}
  # and phase {yphase}.
  #
  # There should be {nct} "L" lines written, with {nlk} links (excluding "None"s) on them.

  debug = False
  
  nwrite_L = 0  # Number of "L" lines written.
  nwrite_ct = 0 # Number of contacts on those lines.
  nwrite_lk = 0 # Number of non-"None" link paths on those lines.
  
  # Write links and contacts between adjacent rasters:
  wr.write("### LINHAS ADJACENTES \n")
  sys.stderr.write("  ### CONTACTS AND LINKS\n")
  SCS = raster.separate_by_scanline(OPHS, (1,0), (0,1), ystep, yphase)
  CTS_found = []
  nsc = len(SCS) # Number of scanlines (including empty ones)
  for isc in range(nsc-1):
    SC0 = SCS[isc]   # List of indices of rasters in scanline {isc}.
    SC1 = SCS[isc+1] # List of indices of rasters in scanline {isc+1}.
    if debug:
      sys.stderr.write("  ====================================\n")
      sys.stderr.write("  scanlines %d %s and %d %s\n" % (isc,str(SC0),isc+1,str(SC1)))
    for iph0 in SC0:
      for iph1 in SC1:
        nwr, ct, nlk01 = write_L_line(wr, OPHS, iph0, iph1, xdir, ydir, shift);
        nwrite_L += nwr
        if ct != None: 
          nwrite_ct += 1
          CTS_found.append(ct)
        nwrite_lk += nlk01
    if debug:
      sys.stderr.write("  ====================================\n")
        
  sys.stderr.write("  expected to write %d L-lines with %d contacts and %d links\n" % (nct, nct, nlk))

  sys.stderr.write("  wrote %d L-lines with %d contacts and %d links\n" % (nwrite_L,nwrite_ct,nwrite_lk))
  
  assert nwrite_ct == len(CTS_found)
  assert nwrite_ct == nct
  assert nwrite_lk == nlk
  assert nwrite_L == nct

  return CTS_found

  # ----------------------------------------------------------------------

def encode_link(olk, xdir, ydir, shift):
  # Encodes the geometry of the link path {olk} as a 
  # list of points spearated by ';', each point 
  # being two coordinates separated by '&'.
  # If {olk} is {None}, returns the string 'None' instead.
  if olk == None:
    xlk = "None"
  else:
    nmv = path.nelems(olk)
    xlk = ""
    for ipt in range(nmv+1):
      mvi = path.elem(olk, min(ipt,nmv-1))
      pti = move.pini(mvi) if ipt < nmv else move.pfin(mvi)
      xlk += (";" if ipt > 0 else "") + encode_point(pti, xdir, ydir, shift, "&")
  return xlk
  # ----------------------------------------------------------------------
    
def write_L_line(wr, OPHS, iph0, iph1, xdir, ydir, shift):
  # Checks whether there is a contact between raster elements
  # {OPHS[iph0]} and {OPHS[iph1]}. The elements are assumed to be
  # single-trace oriented paths on consecutive scanlines, from lower
  # to upper. They may have opposite orientations.
  # 
  # If there is a contact (found through {path.get_contacts}), writes the
  # corresponding 'L' line. If there are link paths between the two
  # elements, writes their geomery on that line, too.
  #
  # Returns the number of "L"-lines (0 or 1), the contact {ct} (or
  # {None}), and the number {nlk01} of links written to the file.
  # 
  # If there is no contact, does not write anything, and returns {None}.
  # There should be no link paths between the two lines in that case.

  debug = False
  verbose = False

  if verbose:
    sys.stderr.write("  ---- %d %d ---------------------------------\n" % (iph0,iph1))

  # Get the upper and lower rasters as path and moves:
  oph0 = OPHS[iph0];  assert path.nelems(oph0) == 1; mv0, dr0 = move.unpack(path.elem(oph0,0)) # Lower.
  oph1 = OPHS[iph1];  assert path.nelems(oph1) == 1; mv1, dr1 = move.unpack(path.elem(oph1,0)) # Upper.

  if debug:
    for isd in range(2):
      ophi = (oph0,oph1)[isd]
      path.show(sys.stderr, ("    oph%d:" % isd), ophi, "\n", True, 0,0,0)

  # Get the link paths {OLKS01} between the two rasters, if any:
  OLKS01 = path.get_all_connecting_links(oph0, oph1)

  # Find the contacts bwtween the two path elements:
  CTS0u = path.get_contacts(oph0,0) # Contacts on upper side of {oph0}.
  CTS1d = path.get_contacts(oph1,1) # Contacts on lower side of {oph1}.
  CTS01 = CTS0u & CTS1d  # Contacts between the two rasters.

  if debug and (len(CTS0u) > 0 or len(CTS1d) > 0):
    sys.stderr.write("~"*120 + "\n")
    CTS0d = path.get_contacts(oph0,1) # Contacts on lower side of {oph0}.
    CTS1u = path.get_contacts(oph1,0) # Contacts on upper side of {oph1}.
    sys.stderr.write("  R%d (%d lower, %d upper contacts)\n" % (iph0,len(CTS0d),len(CTS0u)))
    contact.show_list(sys.stderr, "", tuple(CTS0d), None)
    contact.show_list(sys.stderr, "", tuple(CTS0u), None)
    sys.stderr.write("  R%d (%d lower, %d upper contacts)\n" % (iph1,len(CTS1d),len(CTS0u)))
    contact.show_list(sys.stderr, "", tuple(CTS1d), None)
    contact.show_list(sys.stderr, "", tuple(CTS1u), None)
    sys.stderr.write("  %d common contacts\n" % len(CTS01))
    contact.show_list(sys.stderr, "", tuple(CTS01), None)
    sys.stderr.write("~"*120 + "\n")

  nwr = 0 # Number of "L"-lines written (0 or 1).
  nlk = 0 # Number of links written (0, 1, or 2).
  nct = 0 # Number of contacts written (0 or 1)

  if len(CTS01) == 0:
    # No contact between the rasters.
    ct = None
    assert len(OLKS01) == 0, "links between rasters with no contacts."
  elif len(CTS01) == 1:
    # There is a contact between {oph0} and {oph1}:
    [ct] = CTS01 # Get that lone contact in {ct}.
    nct += 1
    if debug: contact.show(sys.stderr, "ct = ", ct, "\n", 0)

    # Paranoia:
    mv_ct = [None,None]
    for isd in range(2):
      mv_ct[isd] = contact.side_move(ct, isd)
      SDPi = list(contact.get_side_paths(ct, isd))
      if debug:
        sys.stderr.write("paths on side %s:\n" % isd)
        for phj, drj, imvj in SDPi:
          ophj = path.spin(phj, drj)
          path.show(sys.stderr, "  ", ophj, (" imv = %d\n" % imvj), True, 0,0,0) 
      assert len(SDPi) == 1
      phi, dri, imvi = SDPi[0]
      omvi_ph = path.elem(path.spin(phi,dri), imvi)
      mvi_ph, dri_ph = move.unpack(omvi_ph)
      assert mvi_ph == mv_ct[isd]
    assert mv_ct[0] != mv_ct[1]
    assert (mv_ct[0] == mv0 and mv_ct[1] == mv1) or (mv_ct[0] == mv1 and mv_ct[1] == mv0)

    assert len(OLKS01) <= 2, "too many links between two rasters"
    nlk += len(OLKS01)
    XLKS01 = ["None", "None"]
    for ilk in range(len(OLKS01)):
      olki = OLKS01[ilk]
      XLKS01[ilk] = encode_link(olki, xdir, ydir, shift)
    if debug:
      for ilk in range(len(XLKS01)):
        sys.stderr.write("    xlk%d = %s\n" % (ilk,XLKS01[ilk]))

    # Write the contact and links:
    wr.write("L%d,L%d,%s,%s\n" % (iph0,iph1, XLKS01[0],XLKS01[1]))
    nwr += 1
  else:
    assert False, "more than one contact between two rasters"

  if verbose:
    sys.stderr.write("  --- nwr = %d nlk = %d nct = %d ---\n" % (nwr,nlk,nct))
    assert nwr <= 1
    assert nct == nwr
    assert nlk <= 2*nwr

  return nwr, ct, len(OLKS01)
  # ----------------------------------------------------------------------

def write_contours(wr, OCRS, ncr, xdir, ydir, shift):
  # Write to {wr} the contours in {OCRS}, as a bunch of 'C' lines.
  # There should be {ncr} such lines.

  assert ncr == len(OCRS) # Number of contour components.

  wr.write("### CONTORNOS \n")
  sys.stderr.write("  ### CONTOURS\n")

  for icr in range(ncr):
    ocr = OCRS[icr]
    for ipt in range(path.nelems(ocr)):
      p = move.pini(path.elem(ocr, ipt))
      wr.write("C%d,%s\n" % (icr, encode_point(p, xdir, ydir, shift, ",")))
  return
  # ----------------------------------------------------------------------
