# /usr/bin/python3
# Test program for module {paper_figures_B}

import paper_figures_B
import paper_example_B
import move
import move_parms
import path
import contact
import txt_write
import job_parms
import hacks
import pyx
import rn
import sys
from math import sqrt, sin, cos, floor, ceil, inf, nan, pi

# Figures and number of sub-figures:

def test_plot_figure():
  sys.stderr.write("--- testing {make}---\n")

  fignums = (
    
    ( "cover", "basic" ),
    ( "paths", "basic"),
    ( "moves", "dir"),
    ( "moves", "rev"),
    ( "input", "rasct"),
    ( "input", "links"),
    ( "input", "graph"),
    ( "zigzag", "nonalt"),
    ( "zigzag", "alt"),
    ### ( "blocks", "genblocks"),
    ### ( "blocks", "genchoices"),
    ### ( "blocks", "gengraph"),
    ### ( "blocks", "ctrblocks"),
    ### ( "blocks", "ctrchoices"),
    ( "hotfill", "sol50000"),
    ( "hotfill", "sta03500"),
    ( "hotfill", "sol03500"),
    ( "hotfill", "sol02300"),
    ( "hotfill", "sol01700"),
    ( "hotfill", "tab03500"),
    ( "cold", "path"),
    ### ( "rivers", "whole"),
    ### ( "rivers", "subs"),
  )
  
  for fig, subfig in fignums:
    sys.stderr.write("--- testing {plot_figure} fig = '%s' = subfig %s---\n" % (fig,subfig))
    name = "tests/out/paper_figures_B_TST_p_" + fig + "_" + subfig
    paper_figures_B.plot_figure(name, fig, subfig)
  return
  # ----------------------------------------------------------------------

test_plot_figure()
