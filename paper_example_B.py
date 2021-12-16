#! /usr/bin/python3
# Implementation of module {paper_example_B}
# Last edited on 2021-10-21 02:07:49 by stolfi

import paper_example_B_IMP

def make_turkey(mp_cont, mp_fill, mp_link, mp_jump):
  # Returns results {OCRS,OPHS,OLKS,CTS,VGS,EGS} where
  # 
  #   {OCRS} is a list of paths that comprise the contour of a simple
  #   slice.
  #
  #   {OPHS} is a list of paths that comprise a solid raster fill for the slice.
  #
  #   {OLKS} is a list of paths that are links tof the paths in {OPHS}.
  #   
  #   {CTS} is a list of contacts betwen the filling paths.
  #
  #   {VGS} is a list of points that are vertices of the contact graph.
  #
  #   {EGS} is a list of segments that a are the edges of the contact graph.
  #
  # The list {OCRS} always has the outline paths. Each componet of the outline
  # is a separate closed path, oriented so that the interior of the slice is
  # on the left margin of the outline. 
  #
  # The list {OPHS} has the raster traces {TRS}, as one-move paths,
  # oriented left to right.
  #
  # The list {OLKS} has all the potential link paths for the rasters in
  # {OPHS}, oriented bottom to top. The list {CTS} has all the contacts
  # between those rasters. {VGS} and {EGS} are the vertices and edhes of
  # those rasters and contacts. All other lists other than {OCRS} are
  # {None}.
  # 
  # The parameters {mp_cont,mp_fill,mp_jump} are used for the
  # outer contour traces, the filling traces (rasters and links), and the 
  # filling jumps, respectively.
  return paper_example_B_IMP.make_turkey(mp_cont, mp_fill, mp_jump)

def make_simple_moves(mp_fill, mp_jump):
  # Returns a list with a single trace and a single jump, with the specified
  # parameter records.  They will be approximatey horizontal and with the same
  # endpoinst.
  #
  # For compatibility with other procedures in this module, they are
  # packaged as oriented paths with a single move each.
  return paper_example_B_IMP.make_simple_moves(mp_fill, mp_jump)

def make_simple_path(mp_fill, mp_jump):
  # Returns a list with a single oriented path with seven moves, including
  # two jumps. 
  return paper_example_B_IMP.make_simple_path(mp_fill, mp_jump)

def make_simple_contacts(mp_fill):
  # Returns a list {OPHS} of a few paths, and a list {CTS} with a few contacts
  # between them
  return paper_example_B_IMP.make_simple_contacts(mp_fill)

def make_simple_cover(mp_fill, mp_jump):
  # Returns a list {OPHS} of two paths, and a list {CTS} of two contacts,
  # such that {OPHS[0]} covers ide 0 of {CTS[0]} and both sides of {CTS[1]}.
  return paper_example_B_IMP.make_simple_cover(mp_fill, mp_jump)
