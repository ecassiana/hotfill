# Some input data sets for illustrations of the paper - variant "B".
# Last edited on 2021-11-09 16:50:35 by stolfi

import paper_figures_B_IMP

def plot_figure(fname, fig, subfig):
  # Plots the figure {fig}, sub-figure {subfig} for the paper, for the
  # sample slice generated by {paper_example_B.make}. Creates
  # three files called "{fname}.{ext}" where {ext} is "png", "jpg",
  # "eps".
  #
  # The plot uses the paths of {OPHS} connected (or not) in different
  # ways, using some or all of the links, depending on the parametr
  # {fig} and {subfig}, as follows:
  #
  #   "input", 0: Shows the contour paths {OCRS}, the raster paths of
  #   {OPHS}, oriented left to right, and the contacts of {CTS}.
  #
  #   "input", 1: Shows the links of {OLKS}.
  #
  #   "input", 2: Shows the contact graph as described by {VGS} and {EGS}.
  #
  #   "zigzag", 0: Shows the rasters {OPHS}, oriented left to right and
  #   sorted in scanline order, connected by jumps into a single path --
  #   the non-alternating scan-line-order path. Normally it should not
  #   use any links. Also shows a single contact from {CTS}, the one
  #   with possibly max cooling time.
  #
  #   "zigzag", 1: Shows the rasters {OPHS} connected by jumps or links
  #   into a single path -- the alternating scan-line-order path. Also
  #   shows a single contact from {CTS}, the one with possibly max
  #   cooling time.
  #
  #   "cold", 0: Shows the rasters {OPHS} connected by jumps or links
  #   into a single path which should use the minimum number of jumps,
  #   as could be produced by RP3 or slic3r. Also shows one or more
  #   contacts from {CTS} with with possibly excessive cooling time.
  #
  #   "rivers", 0: Shows the rasters {OPHS} as single traces, all oriented left to right and
  #   sorted in scanline order, color-coded according to the partition into rivers.
  #   Also shows the essential cut-lines.
  #
  #   "rivers", 1: Shows the rasters {OPHS} as single traces, all oriented left to right and
  #   sorted in scanline order, color-coded according to the partition into SUB-rivers.
  #   defined by the essential cut-lines. Also shows the essential cut-lines.
  #
  #   "canon", 0: Shows the rasters of {OPHS} combined into two paths,
  #   one a canonical {icu}-path and another one a canonical
  #   {(icu,jcu)}-bandpath. Also shows the essential and non-essential
  #   cut-lines relevant for those paths. Also shows the contacts of
  #   {CTS} that lie betwen the two, on cut-line {icu}.
  #
  #   "blocks", 0: Shows the paths of {OPHS} joined into choices of
  #   blocks that could be produced by {HotFill} for use by
  #   {HotPath}, namely the canonical snakes and their alternatives
  #   with reversed orientation and/or reverse order. The same raster
  #   traces are shared by all choices of each blocks, sometimes
  #   reversed. Shows the cut-lines used to define those canonical
  #   bands. Does not show contacts or links.
  #
  #   "blocks", 1: Shows a set of arbitrary blocks as could be used
  #   as input to {HotPath} directly (not generated by {HotPath).
  #   Also shows the links of {OLKS} that can be used to connect
  #   those blocks. Also shows the subset of the input contacts that lie between 
  #   between different blocks.   One of the blocks is highlighted.
  #
  #   "blocks", 2: Shows only the four choices of the highlighted block in
  #   the plot ("blocks",1).
  #
  # Where needed, the plot will assume that the contour is drawn with
  # traces of nominal width {wdc} and the filling rasters with nominal
  # width {wdf}.
  paper_figures_B_IMP.plot_figure(fname, fig,subfig)
