# Reading slice elements from an RP3-generated text file.
# Last edited on 2021-11-25 18:13:36 by stolfi

import txt_read_IMP

def read(rd, mp_cont, mp_fill, mp_link, angle, shift, fixdir, smptol):
  # Reads from the file handle {rd} a description of the contours, 
  # raster filling elements, contacts, and link paths of a single object slice
  # with a simple raster-type solid fill.
  #
  # Returns 
  #
  #    a list {OCRS} of contour paths
  #
  #    a list {OPHS} of filling paths
  #
  #    a list {OLKS} of link paths between those filling paths
  #
  #    a list {CTS} of contacts between the filling elements.
  # 
  #    the nominal {Z} coordinate of the slice.  
  # 
  # ??? Describe the file format ???
  #
  # The contour elements (defined by the lines with 'C' code) are closed
  # paths, properly oriented, with trace parameters {m_cont}. The
  # procedure {path.compute_contour_nesting} is called to compute and
  # save their nesting relations.
  # 
  # The filling elements must be parallel raster traces whose axis makes
  # the specified {angle} (in radians, counterclockwise) with the
  # horizontal.  The words "left", "right", "above" and "below" in this 
  # spec are relative to the coordinate system rotaded by that {angle}.
  #
  # Each filling element ('R' code) is converted to a single {Move}
  # object with parameters {mp_fill}, oriented left-to-right, and then
  # to a single-move {Path} object. The {OPHS} result is the list of
  # those {Path} objects, in the order they occur in the file.
  #
  # The {Move} is always oriented from left to right. If {fixdir} is
  # false, the {Path} will use the {Move} in the orientation specified
  # in the file by the original order of the endpoints and by the 'rbit'
  # field of its 'R' line. If {fixdir} is true, the {Path} will use the
  # move in its native (left-to-right) orientation, ignoring the order
  # of the endpoints in the file and the 'rbit' flag.
  # 
  # The procedure saves in each {Path} object of {OPHS} the 'group'
  # index from the respective 'R' line, with {path.set_group}. This
  # index is supposed to indicate a user-given grouping of rasters into
  # blocks such as the "continuous raster sequences" created by RP3.
  # 
  # The lines with 'L' code are converted to contacts and link paths
  # between the filling elements. Each link path will consist of one or
  # more traces with the same parameters {mp_link}. These links are
  # attached to the raster paths with {path.add_link} so that they can
  # be obtained with {path.get_links} and {path.connecting_link}.
  #
  # If {smptol} is positive, the link paths are simplified with
  # {path.simplify} and that tolerance parameter. The siplification will
  # preserve the endpoits of the link. If {smptol} is zero, the link
  # paths will be as given in the file, including any zero-length
  # traces.
  #
  # The {Contact} object {ct}, if any, and its two adjacent paths
  # {oph[0]} (below) and {oph[1]} (above) are also associated through
  # {contact.add_side_path(ct,isd,oph[isd])},
  # {contact.add_side_path(ct,isd,rev(oph[isd]))}, and
  # {path.add_contact(oph[isd],isd,ct)}, for {isd} in {0..1}.
  #
  # In the returned data sets, all points read from the file will be
  # rotated by {-angle} radians about the origin, and then translated by
  # the vector {shift} (a pair of {float}s). Note that rotation may
  # change the sign of some coordinates. A suitable {shift} may be used
  # to ensure that all coordinates in the resulting paths and contacts
  # are non-negative, if so desired. Thus, in the returned data set, the
  # rasters will be horizontal and sorted by increasing Y then
  # increasing X.
  # 
  return txt_read_IMP.read(rd, mp_cont, mp_fill, mp_link, angle, shift, fixdir, smptol)

