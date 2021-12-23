#! /usr/bin/python3
# Last edited on 2021-11-25 18:17:49 by stolfi

import txt_write
import txt_read
import job_parms
import move_parms
import move
import path
import block
import raster_example
import hacks; from hacks import fbomb
import contact
import input_data

import rn

import sys
from math import sqrt, sin, cos, log, exp, floor, ceil, inf, nan, pi

parms = job_parms.typical_js()
mp_jump = move_parms.make_for_jumps(parms)

# Some arbitrary dynamics parameters:
ac = parms['acceleration'] 
sp = parms['max_extrusion_speed'] 

# Move parameters matching the raster spacings in the file:
wd_cont = 0.75; mp_cont = move_parms.make(wd_cont, ac, sp, 0.0)
wd_fill = 1.00; mp_fill = move_parms.make(wd_fill, ac, sp, 0.0)
wd_link = 0.50; mp_link = move_parms.make(wd_link, ac, sp, 0.0)

parms['solid_raster_width'] = wd_fill
parms['contour_trace_width'] = wd_cont

zstep = parms['slice_thickness']

tol = 1.5e-5 # Because most numbers are written with 6 decimal fraction digits.

def compare_contour_lists(OCRS_org, OCRS_rdb):
  
  sys.stderr.write("comparing contour lists...\n")

  debug = False

  def show_nesting(tag, ocr):
    path.show(sys.stderr, ("  children of %s contour " % tag), ocr, "\n", False, 0,0,0)
    CHILDREN = path.inner_contours(ocr)
    if CHILDREN == None:
      sys.stderr.write("    None\n")
    else:
      path.show_list(sys.stderr, "    ", CHILDREN, None, False, False)
    sys.stderr.write("\n")
    return
    # ....................................................................

  ncr = len(OCRS_org)
  # Compare the original and readback contours:
  assert len(OCRS_rdb) == ncr
  # sys.stderr.write("OCRS_rdb = %s\n" % str(OCRS_rdb))
  for icr in range(ncr):
    if debug: sys.stderr.write("~"*100 + "\n")
    ocri_org = OCRS_org[icr]
    ocri_rdb = OCRS_rdb[icr]
    if debug:
      show_nesting("org", ocri_org)
      show_nesting("rdb", ocri_rdb)
    path.compare(ocri_org, ocri_rdb, tol, die=True)
    # Check if nesting information is the same:
    for jcr in range(ncr):
      ocrj_org = OCRS_org[jcr]
      ocrj_rdb = OCRS_rdb[jcr]
      ne_org = path.contour_nesting(ocri_org, ocrj_org)
      ne_rdb = path.contour_nesting(ocri_rdb, ocrj_rdb)
      assert ne_org == ne_rdb
    if debug: sys.stderr.write("~"*100 + "\n")
  return
  # ----------------------------------------------------------------------

def path_key(oph):
  # A "fingerprint" of a path that is hopefully independent
  # of orientation and roundoff errors:

  p = path.pini(oph); pp = (floor(0.01*p[0]/tol + 0.5), floor(0.01*p[1]/tol + 0.5)) 
  q = path.pfin(oph); qq = (floor(0.01*q[0]/tol + 0.5), floor(0.01*q[1]/tol + 0.5))
  if pp < qq:
    return (pp, qq)
  else:
    return (qq, pp)
  # ....................................................................

def contact_key(ct):
  # A "fingerprint" of a contact that is hopefully independent
  # of orientation and roundoff errors:
  p, q = contact.endpoints(ct)
  m = rn.mix(0.5,p, 0.5,q)
  mm = (floor(0.01*m[0]/tol + 0.5), floor(0.01*m[1]/tol + 0.5)) 
  return mm
  # ....................................................................
      
def compare_path_lists(OPHS_org, OPHS_rdb, rasters):
  # Compares the original and readback list of paths:

  debug = False
  
  title = "rasters" if rasters else "link"
  sys.stderr.write("comparing %s lists...\n" % title)

  nph = len(OPHS_org)

  if debug:
    sys.stderr.write("OPHS_org: (%d)\n" % len(OPHS_org))
    path.show_list(sys.stderr, "  ", OPHS_org, None, True, False)
    sys.stderr.write("OPHS_rdb: (%d)\n" % len(OPHS_rdb))
    path.show_list(sys.stderr, "  ", OPHS_rdb, None, True, False)
  
  OPHS_org = sorted(OPHS_org, key = path_key)
  OPHS_rdb = sorted(OPHS_rdb, key = path_key)
  
  if debug:
    sys.stderr.write("OPHS_org sorted: (%d)\n" % len(OPHS_org))
    path.show_list(sys.stderr, "  ", OPHS_org, None, True, False)
    sys.stderr.write("OPHS_rdb sorted: (%d)\n" % len(OPHS_rdb))
    path.show_list(sys.stderr, "  ", OPHS_rdb, None, True, False)
  
  assert len(OPHS_rdb) == nph

  for iph in range(nph):
    ophi_org = OPHS_org[iph]
    ophi_rdb = OPHS_rdb[iph]
    if debug:
      sys.stderr.write("comparing paths {OPHS_org[%d]}, {OPHS_rdb[%d]}\n" % (iph,iph))
      path.show(sys.stderr, "  ophi_org = ", ophi_org, "\n", False, 0,0,0)
      path.show(sys.stderr, "  ophi_rdb = ", ophi_rdb, "\n", False, 0,0,0)
    # The native orientations of the paths may still be opposite.
    # And writing and reading back may have removed points that are too close.
    key_org = path_key(ophi_org)
    key_rdb = path_key(ophi_rdb)
    if debug:
      sys.stderr.write("ophi_org = %s\n" % str(key_org))
      sys.stderr.write("ophi_rdb = %s\n" % str(key_rdb))
    assert key_org == key_rdb, "links do not match"

    if rasters:
      # Check if group information is the same:
      assert path.get_group(ophi_org) == path.get_group(ophi_rdb) 

      # Check if link lists are the same:
      for krv in range(2):
        ophik_org = path.spin(ophi_org, krv)
        ophik_rdb = path.spin(ophi_rdb, krv)
        LKS_org = path.get_links(ophik_org)
        LKS_rdb = path.get_links(ophik_rdb)
        if debug:
          sys.stderr.write("  krv = %d LKS_org (%d) = %s\n" % (krv,len(LKS_org),str(LKS_org)))
          sys.stderr.write("  krv = %d LKS_rdb (%d) = %s\n" % (krv,len(LKS_rdb),str(LKS_rdb)))
        assert len(LKS_org) == len(LKS_rdb)
        # ??? Should check if link paths are the same. ???
  return
  # ----------------------------------------------------------------------
    
def compare_contact_lists(CTS_org, CTS_rdb):
  # Compares the original and readback contacts:
  
  debug = False
  
  sys.stderr.write("comparing contact lists...\n")
  
  if debug:
    sys.stderr.write("CTS_org: (%d)\n" % len(CTS_org))
    contact.show_list(sys.stderr, "  ", CTS_org, None)
    sys.stderr.write("CTS_rdb: (%d)\n" % len(CTS_rdb))
    contact.show_list(sys.stderr, "  ", CTS_rdb, None)
  
  CTS_org = sorted(CTS_org, key = contact_key)
  CTS_rdb = sorted(CTS_rdb, key = contact_key)
  
  if debug:
    sys.stderr.write("CTS_org sorted: (%d)\n" % len(CTS_org))
    contact.show_list(sys.stderr, "  ", CTS_org, None)
    sys.stderr.write("CTS_rdb sorted: (%d)\n" % len(CTS_rdb))
    contact.show_list(sys.stderr, "  ", CTS_rdb, None)

  nct = len(CTS_org)
  assert len(CTS_rdb) == nct
  for ict in range(nct):
    # Check if contact lists are the same:
    cti_org = CTS_org[ict]
    cti_rdb = CTS_rdb[ict]
    key_org = contact_key(cti_org)
    key_rdb = contact_key(cti_rdb)
    if debug:
      sys.stderr.write("cti_org = %s\n" % str(key_org))
      sys.stderr.write("cti_rdb = %s\n" % str(key_rdb))
    assert key_org == key_rdb, "contacts do not match"
  return
  # ----------------------------------------------------------------------
   
def do_show_data(title, OCRS, OPHS, OLKS, CTS, Z):
  # Prints the given contour paths {OCRS}, fill paths {OPHS}, link paths {OLKS}, and contacts {CTS},
  # as well as the nominal Z-coordinate {Z}.
  
  sys.stderr.write("### %s #########################################################################\n" % title)
  sys.stderr.write("Z = %.3f:\n" % Z)
  sys.stderr.write("CONTOURS (%d):\n" % len(OCRS))
  path.show_list(sys.stderr, "  ", OCRS, None, False, True)
  sys.stderr.write("RASTERS (%d):\n" % len(OPHS))
  path.show_list(sys.stderr, "  ", OPHS, None, True, False)
  sys.stderr.write("LINKS (%d):\n" % len(OLKS))
  path.show_list(sys.stderr, "  ", OLKS, None, True, True)
  sys.stderr.write("CONTACTS (%d):\n" % len(CTS))
  contact.show_list(sys.stderr, "  ", CTS, None)
  return 
  # ----------------------------------------------------------------------

def do_test_write(outname, OCRS_org, BCS_org, Z_org, mp_cont, mp_fill, mp_link, angle, shift, debug):

  outfolder = "tests/out/"
  
  OPHS_org = []
  OLKS_org = set()
  CTS_org = set()
  for ibc in range(len(BCS_org)):
    bc = BCS_org[ibc]
    # Take only the first choice of each block:
    oph = block.choice(bc, 0)
    # Make sure it is a single trace:
    assert path.nelems(oph) == 1
    assert not move.is_jump(path.elem(oph,0))
    path.set_group(oph, ibc)
    OPHS_org.append(oph)
    for olk in path.get_links(oph) + path.get_links(path.rev(oph)):
      lk, drlk = path.unpack(olk)
      OLKS_org.add(lk)
    for ct in set.union(path.get_contacts(oph, 0), path.get_contacts(oph, 1)):
      CTS_org.add(ct)
  OLKS_org = list(OLKS_org)
  CTS_org = list(CTS_org)
  
  ncr_org = len(OCRS_org)
  nph_org = len(OPHS_org)
  nlk_org = len(OLKS_org)
  nct_org = len(CTS_org)
  sys.stderr.write("original data: %d contours, %d rasters, %d links, %d contacts\n" % (ncr_org,nph_org,nlk_org,nct_org))

  if debug: do_show_data("original data", OCRS_org, OPHS_org, OLKS_org, CTS_org, Z_org)

  fname_out = outfolder + outname + ".txt"

  sys.stderr.write("writing %s ...\n" % fname_out)
  wr = open(fname_out, "w")
  txt_write.write(wr, OCRS_org, OPHS_org, Z_org, angle, shift)
  wr.close()
  
  sys.stderr.write("reading back %s ...\n" % fname_out)
  rd = open(fname_out, "r")
  smptol = 0
  fixdir = False
  OCRS_rdb, OPHS_rdb, OLKS_rdb, CTS_rdb, Z_rdb = txt_read.read(rd, mp_cont, mp_fill, mp_link, angle, shift, fixdir, smptol)
  rd.close()
  
  if debug: do_show_data("read-back data", OCRS_rdb, OPHS_rdb, OLKS_rdb, CTS_rdb, Z_rdb)

  assert len(OLKS_rdb) == len(OLKS_org)
      
  # Check the links (paranoia):
  OLKS_chk = set()  # Set of all {Path} objects of links:
  for ophi in OPHS_rdb:
    for ophir in ophi, path.rev(ophi):
      OLKSir = path.get_links(ophir)
      for olkirj in OLKSir:
        lkirj, drirj = path.unpack(olkirj)
        OLKS_chk.add(lkirj)
  OLKS_chk = list(OLKS_chk) # Convert from set to list.
  assert len(OLKS_rdb) == len(OLKS_chk)
  
  compare_contour_lists(OCRS_org, OCRS_rdb)
  compare_path_lists(OPHS_org, OPHS_rdb, rasters=True)
  compare_path_lists(OLKS_org, OLKS_rdb, rasters=False)
  compare_contact_lists(CTS_org, CTS_rdb)
  
  # Compare the original and readback {Z_org} coordinate:
  assert abs(Z_org - Z_rdb) <= 1.5e-4 # Written with 5 decimals by {RP3}.

  # ??? Should validate output -- contours, links, contacts, etc. ???
  
  BCS = [ block.from_paths((oph, path.rev(oph))) for oph in OPHS_org ]

  CLRS = hacks.trace_colors(len(BCS), None)
  wd_axes = 0.15*wd_fill
  deco = (len(BCS) < 100)
  links = True
  contacts = True
  input_data.plot_input(outfolder + outname, OCRS_org, BCS, None, CLRS, wd_axes,deco, links, contacts)  
  return
  # ----------------------------------------------------------------------

def test_write_from_read(partname, islice, angle, shift):

  sys.stderr.write("--- testing {write} from {txt_read.read} ---\n")

  infolder = "tests/in/2021-05-15-elis/" + partname + "/"
  sys.stderr.write("infolder = %s\n" % infolder)
  sys.stderr.write("partname =  %s\n" % partname)
  sys.stderr.write("islice =  %s\n" % islice)
  sys.stderr.write("angle =  %s\n" % angle)
  
  tag = partname + "_" + ("%03d" % islice)
  
  fname_in = infolder + tag + ".txt"

  # Some arbitrary dynamics parameters:
  ac = parms['acceleration'] 
  sp = parms['max_extrusion_speed'] 

  # Move parameters matching the raster spacings in the file:
  wd_cont = 2.00; mp_cont = move_parms.make(wd_cont, ac, sp, 0.0)
  wd_fill = 3.00; mp_fill = move_parms.make(wd_fill, ac, sp, 0.0)
  wd_link = 2.00; mp_link = move_parms.make(wd_link, ac, sp, 0.0)

  sys.stderr.write("reading %s ...\n" % fname_in)
  rd = open(fname_in, "r")
  smptol = 0
  fixdir = False # Orient the raster elements as in the file.
  OCRS_org, OPHS_org, OLKS_org, CTS_org, Z_org = txt_read.read(rd, mp_cont, mp_fill, mp_link, angle, shift, fixdir, smptol)
  rd.close()
  
  outname = "txt_write_TST_read_" + tag
  BCS_org = [ block.from_paths([oph,]) for oph in OPHS_org ]
  do_test_write(outname, OCRS_org, BCS_org, Z_org, mp_cont, mp_fill, mp_link, angle, shift, True)
  return
  # ----------------------------------------------------------------------
  
def test_write_synthetic(dsname, variant):

  sys.stderr.write("--- testing {write} of synthetic data ---\n")

  sys.stderr.write("dsname =  %s\n" % dsname)
  sys.stderr.write("variant =  %s\n" % variant)

  tag = dsname + "_" + variant

  # Some arbitrary dynamics parameters:
  ac = parms['acceleration'] 
  sp = parms['max_extrusion_speed'] 

  # Move parameters matching the raster spacings in the file:
  wd_cont = 0.75; mp_cont = move_parms.make(wd_cont, ac, sp, 0.0)
  wd_fill = 1.00; mp_fill = move_parms.make(wd_fill, ac, sp, 0.0)
  wd_link = 0.50; mp_link = move_parms.make(wd_link, ac, sp, 0.0)
  
  OCRS_org, BCS_org, o_org, Z_org = input_data.make_synthetic(dsname, variant, mp_cont, mp_fill, mp_link, mp_jump)
  
  outname = "txt_write_TST_syn_" + tag
  angle = pi/6
  shift = (-60,30)
  do_test_write(outname, OCRS_org, BCS_org, Z_org, mp_cont, mp_fill, mp_link, angle, shift, True)
  return
  # ----------------------------------------------------------------------

test_write_synthetic(dsname = "paper_fig_B", variant = "rasters")

test_write_synthetic(dsname = "patch_array", variant = "1x1x7no")

test_write_synthetic(dsname = "patch_array", variant = "1x1x27no")
test_write_synthetic(dsname = "patch_array", variant = "2x1x27no")
test_write_synthetic(dsname = "patch_array", variant = "3x1x27no")

test_write_synthetic(dsname = "patch_array", variant = "2x1x27is")

test_write_from_read(partname = "chain_link_2", islice = 2, angle = pi/2, shift = (0,60))

