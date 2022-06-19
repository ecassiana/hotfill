# Tools for plotting input datasets of the {HotPath} algorithm.

import input_data_IMP
import path

def make_synthetic(dsname, variant, mp_cont, mp_fill, mp_link, mp_jump):
  # Returns a synthetic input data set for {HotPath} algorithms,
  # consisting of a lists of contour paths {OCRS}, a list of blocks {BCS},
  # a starting point {o}, and the nominal {Z} coordinate of the slice.
  # 
  # The paths in {OCRS} will be contour paths that are not included in
  # the blocks {BCS}. Their traces will have parameters {mp_cont}. The
  # traces in the blocks of {BCS} will have parameters {mp_fill},
  # {mp_cont}, and {mp_link}, depending on their nature.
  #
  # The procedure may also generate link paths that can be used to
  # connect endpoints of choices of the blocks. The traces in these
  # paths will have parameters {mp_link}. They can be obtained with
  # {path.get_links} applied to the block choices.
  #
  # The procedure may also generate contacts between the traces of
  # choices of distinct blocks. They can be obtained with
  # {path.get_contacts} applied to the block choices.
  return input_data_IMP.make_synthetic(dsname, variant, mp_cont, mp_fill, mp_link, mp_jump)

def plot_input(fname, OCRS, BCS, o, CLRS, wd_axes,deco, links, contacts):
  # Plots the input data for the HTPP problem.
  # Assumes that {OCRS} is a list of paths that forms the slice's contour,
  # and {BCS} is a list of blocks.
  #
  # If {links} is true, will plot the links between blocks, obtained by
  # {path.get_links} from the choices of {BCS}.
  #
  # If {contacts} is true, will plot any contacts that are attached to
  # the choices of {BCS}, obtained by {path.get_contacts}.
  #
  # The paraneter {CLRS} must be a list of colors to use for the traces of the 
  # choices of {BCS}. If it has only one elelemnt, all traces are painted 
  # with that color.  Otherwise it must have the same length as {BCS}, and the traces in each
  # block will be painted with the corresponding color. 
  #
  # If {deco} is true also plots the matter footprint and draws
  # axes, dots, and arrowheads on the traces.
  input_data_IMP.plot_input(fname, OCRS, BCS, o, CLRS, wd_axes,deco, links, contacts)

def show_check_and_plot(fname, OCRS, OPHS, OLKS, CTS, o, ydir, wd_axes, deco):
  # Prints out the oriented paths in the lists {OCRS} and {OPHS} and
  # plots them to files "{fname}.{ext}" where {ext} is "jpg", "png",
  # etc. Assumes that the paths in {OCRS} are contours, and those in
  # {OPHS} comprise the filling
  #
  # Also checks whether the link paths attached to the paths in
  # {OPHS}with {path.get_links} and the contacts attached with
  # {path.get_contacts,contact.get_side_paths} are consistent with the
  # lists {OLKS} of oriented link paths and {CTS} of {Contat} objects,
  # respectively.
  #
  # If {ydir} is not {None}, also checks whether the sides of each contact
  # are properly ordered by theri projections on {ydir}.
  #
  # If {deco} is true also plots the matter footprint and draws
  # axes, dots, and arrowheads on the traces.
  input_data_IMP.show_check_and_plot(fname, OCRS, OPHS, OLKS, CTS, o, ydir, wd_axes, deco)
