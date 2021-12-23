# Some example sets of raster filling elements for tests and illustrations.
# Last edited on 2021-10-12 19:03:36 by stolfi

import raster_example_IMP

# MISCELLANEOUS RASTER FILLS FOR TESTING

def patch_array(cols, rows, nx,ny, islands, mp_cont,mp_fill,mp_link): 
  # Returns
  #
  #   {OCRS}, a list of contours;
  #
  #   {OPHS}, a list of filling rasters;
  #
  #   {OLKS}, a list of link paths between those rasters;
  #
  #   {CTS}, a list of contacts between those rasters.
  #
  # Each element of {OPHS} is a {Path} object that consists of a single
  # /raster/ (horizontal trace) with parameters {mp_fill}. Both trace
  # and path are oriented left to right. The rasters are sorted
  # in /scanline order: by increasing {Y} coordinate, and 
  # increasing {X} coordinate in each scanline.
  #
  # The rasters are arranged into a two-dimensional array of rectangular
  # raster stacks (/patches/), with long single rasters (/beams/)
  # before, between, and after the rows of that array. The array has
  # {rows} rows and {cols} columns.
  #
  # Every patch has {ny}, {ny-2}, or {ny-4} rasters. The trace enpoints
  # of each patch lie on the boundary of a rectangle whose width is
  # {nx*wdf} and whose heights is atmost {ny*wdf}, where {wdf} is
  # {move_parms.width(mp_fill)}. The number {nx} must be positive,
  # and {ny} must be at least 5.
  #
  # If {islands} is true, {cols} must be congruent to 1 modulo 4. Inn
  # column 0 of the array, every patch has {ny} rasters and makes
  # contach with the two adjacent beams. In column 1, every patch makes
  # contact with the lower beam, but has height {ny-2} so it does not
  # touch the upper beam. In column 2, the patches have {ny-4} rasters
  # and don't touch either beam. In column 3, the patches have {ny-2}
  # rasters too, but touch the upper beam only. The pattern then
  # repeats, every 4 array columns.
  #
  # If {islands} is false, cols must be congruent to 1 modulo 3. The
  # arrangement is the same, except that the island patch (column 2) is
  # omitted, and the pattern repeats every 3 columns.
  #
  # The procedure also creates a set of /link paths/, whose traces have
  # parameters {mp_link}, that can be used to connect the filler
  # elements. Each link connects two endpoints of two rasters in
  # adjacent scan-lines that are vertically aligned. The links that
  # connect to a filler path {oph} can be accessed through
  # {path.get_links(oph)} and {path.get_links(rev(oph))}.
  # 
  # Each raster will be assigned {path.set_group} to a separate
  # group. The group indices will be consecutive, starting at zero.
  #
  # The contacts in the list {CTS} will also be attached to the 
  # raster paths of {OPHS} and can be obtained with {path.get_contacts}.
  # For each contact {ct}, the move {contact.side_move(ct,0)} will be
  # below the contact, and {contact.side_move(ct,1)} will be above it.
  #
  # The contours in the list {OCRS} will use parameters {mp_cont}.
  # They will be properly oriented (ccw for islands, cw for holes) and 
  # their nesting will be available through {path.inner_contours} and
  # {path.outer_contours}
  return raster_example_IMP.patch_array(cols, rows, nx,ny, islands, mp_cont,mp_fill,mp_link)

def rasters_A(mp_fill, xdir, ydir, yphase, eps):
  # Returns a list of nine raster filling elements parallel to direction
  # {xdir}, in scanline order. The traces have parameters {mp_fill}
  # and lie on four scanlines whose mean spacing is close to the common
  # width {wd} of those traces, with scanline phase {yphase} (mm).
  # The rasters are slightly perturbed from the ideal scanline positions
  # by up to {eps} millimeters.
  return raster_example_IMP.rasters_A(mp_fill, xdir, ydir, yphase, eps)

def rasters_B(nph, xdir, mp_trace, mp_link):
  # Returns, a list {TRS} of {nph+2} traces, a list {PHS} of {nph}
  # raster-like {Path} objects, a list {LKS} of {2*(nph+1)} link paths
  # connecting the endpoints of those traces, and a list {CTS} of {np+1}
  # contacts between them.
  # 
  # All traces will have parameters {mp_trace}, will be horizontal and
  # with the same X span, will be oriented left to right, and will be
  # separated in Y by the width of {mp_trace}. The paths {PHS[0.nph-1]}
  # use the traces {TRS[1..nph]}, respectively; one trace each, in its
  # native orientation.
  # 
  # Each link path will consist of two traces with parameters {mp_link},
  # bent in the middle. These link paths (and their reverses) are
  # attacked to the paths of {PHS} so that they can be obtained with
  # {path.get_links}.
  # 
  # Each contact {CTS[ict]} has traces {TRS[ict]} and {TRS[ict+1]} as
  # sides 0 and 1, respectively.
  return raster_example_IMP.rasters_B(nph, xdir, mp_trace, mp_link)
