#! /usr/bin/python3
# Test program for module {raster_regroup}.

import raster_regroup
import block
import path
import move
import move_parms
import contact
import contact_example
import raster
import raster_example
import hacks
import job_parms

import color
import rn
import pyx
import sys
from math import sqrt, sin, cos, floor, ceil, inf, nan, pi

parms = job_parms.typical()
parms['solid_raster_width'] = 1.0 
job_parms.write(sys.stderr, None, parms, None)

mp_jump = move_parms.make_for_jumps(parms)
mp_cont = move_parms.make_for_contours(parms)
mp_fill = move_parms.make_for_fillings(parms)

wd_fill = move_parms.width(mp_fill)
wd_cont = move_parms.width(mp_cont)

def plot_by_group(OPHS, islands, iplot, tag):
  # Plots the paths colored according to their groups. There must be no {None} groups.
  # Plot files are "tests/out/raster_regroup_TST_is{islands}_{iplot}_{tag}.{ext}"
  # where {ext} is "png", "jpg", "eps".
  
  # Find number of goups {ngr}:
  maxgr = -1
  for oph in OPHS: 
    igr = path.get_group(oph)
    assert igr != None, "group index is {None}" 
    maxgr = max(maxgr, igr)
  ngr = maxgr + 1
  sys.stderr.write("found %d non-empty groups\n" % ngr)
  
  # Generate colors for groups:
  CLRS_gr = hacks.trace_colors(ngr, None)
  
  # Assign path colors according to group:
  CLRS_ph = [ CLRS_gr[path.get_group(oph)] for oph in OPHS ]
  
  # Plot them:
  rwd = 0.80
  wd_axes = 0.05*wd_fill
  grid = False
  deco = False
  pretag = ("is%s" % "FT"[islands]) + "_" + ("%02d" % iplot)
  fname = "tests/out/raster_regroup_TST_" + pretag + "_" + tag
  path.plot_to_files(fname, OPHS, CLRS_ph, rwd, wd_axes, grid, deco)
  return
  # ----------------------------------------------------------------------

def test_regroups(islands):
  sys.stderr.write("--- testing {merge_all,by_contours,split_at_forks,split_by_size} ---\n")
  
  xdir = (1, 0); ydir = (0, 1) # Defined by {patch_array}
  
  cols = 13 # Good for islands and without islands.
  rows = 2
  nx = 3
  ny = 2 + 7 + 2  
  mp_link = mp_fill
  OCRS, OPHS, OLKS, CTS = raster_example.patch_array(cols, rows, nx,ny, islands, mp_cont, mp_fill, mp_link)
  ncr = len(OCRS)  # Number of contour paths.
  nph = len(OPHS)  # Number of raster fill elements.
  nlk = len(OLKS)  # Number of link paths between those rasters.
  nct = len(CTS)   # Number of contacts.

  # Make sure that the rasters are sorted by scanline:
  ystep, yphase = raster.get_spacing_and_phase(OPHS, xdir, ydir)
  OPHS = raster.sort_by_scanline(OPHS, xdir, ydir, ystep, yphase)
  
  # Make sure all group indices are defined:
  for oph in OPHS: assert path.get_group(oph) != None
  
  sys.stderr.write("  ... raster fill elements ...\n")
  path.show_list(sys.stderr, "  ", OPHS, None, True, True)
  
  # Save original groups:
  GRS = [ path.get_group(oph) for oph in OPHS ]
  
  npl =  0 # Number of plots created.
  
  sys.stderr.write("  ... original groups ...\n")
  plot_by_group(OPHS, islands, npl, "orginal")
  npl += 1
  
  sys.stderr.write("  ... {merge_all} ...\n")
  for k in range(nph): path.set_group(OPHS[k], GRS[k])
  raster_regroup.merge_all(OPHS)
  plot_by_group(OPHS, islands, npl, "merge_all")
  npl += 1
  
  sys.stderr.write("  ... {by_contours} ...\n")
  for k in range(nph): path.set_group(OPHS[k], GRS[k])
  raster_regroup.by_contours(OPHS, OCRS)
  plot_by_group(OPHS, islands, npl, "by_contour")
  npl += 1
  
  sys.stderr.write("  ... {split_at_forks} from original ...\n")
  for k in range(nph): path.set_group(OPHS[k], GRS[k])
  raster_regroup.split_at_forks(OPHS)
  plot_by_group(OPHS, islands, npl, "ORG_split_at_forks")
  npl += 1
  
  sys.stderr.write("  ... {split_at_forks} after {merge_all} ...\n")
  raster_regroup.merge_all(OPHS)
  raster_regroup.split_at_forks(OPHS)
  plot_by_group(OPHS, islands, npl, "MRG_split_at_forks")
  npl += 1
  
  # Calling {split_at_forks} should give the same result after {by_contours} or {merge_all}. 
  
  max_lines = 3
  
  sys.stderr.write("  ... {split_by_size} from original ...\n")
  for k in range(nph): path.set_group(OPHS[k], GRS[k])
  raster_regroup.split_by_size(OPHS, ydir, max_lines)
  plot_by_group(OPHS, islands, npl, "ORG_split_by_size")
  npl += 1
  
  sys.stderr.write("  ... {split_by_size} after {merge_all} ...\n")
  raster_regroup.merge_all(OPHS)
  raster_regroup.split_by_size(OPHS, ydir, max_lines)
  plot_by_group(OPHS, islands, npl, "MRG_split_by_size")
  npl += 1
  
  sys.stderr.write("  ... {split_by_size} after {by_contours} ...\n")
  raster_regroup.by_contours(OPHS, OCRS)
  raster_regroup.split_by_size(OPHS, ydir, max_lines)
  plot_by_group(OPHS, islands, npl, "OCRS_split_by_size")
  npl += 1
  
  sys.stderr.write("  ... {split_by_size} after {split_at_forks} ...\n")
  for k in range(nph): path.set_group(OPHS[k], GRS[k])
  raster_regroup.split_at_forks(OPHS)
  raster_regroup.split_by_size(OPHS, ydir, max_lines)
  plot_by_group(OPHS, islands, npl, "SPL_split_by_size")
  npl += 1
  
  sys.stderr.write("  ... {split_by_size} after {merge_all,split_at_forks} ...\n")
  for k in range(nph): path.set_group(OPHS[k], GRS[k])
  raster_regroup.merge_all(OPHS)
  raster_regroup.split_at_forks(OPHS)
  raster_regroup.split_by_size(OPHS, ydir, max_lines)
  plot_by_group(OPHS, islands, npl, "MRG_SPL_split_by_size")
  npl += 1

  sys.stderr.write("  ... {split_by_size} after {by-contour,split_at_forks} ...\n")
  for k in range(nph): path.set_group(OPHS[k], GRS[k])
  raster_regroup.by_contours(OPHS, OCRS)
  raster_regroup.split_at_forks(OPHS)
  raster_regroup.split_by_size(OPHS, ydir, max_lines)
  plot_by_group(OPHS, islands, npl, "CRS_SPL_split_by_size")
  npl += 1

  # ??? More cases? ???

  return
  # ----------------------------------------------------------------------

test_regroups(False)
test_regroups(True)
 

