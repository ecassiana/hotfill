# Writing slice elements to an RP3-style text file.
# Last edited on 2021-10-10 12:56:30 by stolfi

import txt_write_IMP

def write(wr, OCRS, OPHS, Z, angle, shift):
  # Writes to file handle {wr} a description of a single object slice at
  # nominal height {Z}, consisting of the contours {OCRS}, the filling
  # elements {OPHS}. Also writes the link paths and contacts associated to
  # the paths in {OPHS}.
  # 
  # ??? Describe the file format ???
  #
  # Each element {ocr} in the list {OCRS} (written as a line with the 'C'
  # code) must be a /contour/, a closed path without jumps that follows
  # the boundary of the slice. The contours must be properly oriented
  # and must have the same trace parameters.
  # 
  # Each element {oph} in the filling list {OPHS} (written as a line
  # with 'R' code) must be a /raster/, an oriented path with a single
  # trace move. The axis of each raster must be horizontal and oriented
  # from left to right, and the rasters must be sorted in normal
  # scan-line order (increasing Y then increasing X).
  #
  # Each raster {oph} of {OPHS} is written out with the first endpoint
  # being the one that would be the leftmost one in the urotated view,
  # independently of its orientation. However, its orientation with
  # respect to the left-to-right one is recorded in the 'rbit' field of
  # the line. The procedure also writes in the 'R' line the 'group'
  # index assocated to {oph} by {path.get_group}.
  # 
  # The procedure also writes a line with 'L' code for each pair of
  # rasters that share either a contact or a connecting link path. The
  # links that end at each endpoint are obtaned with
  # {path.get_links(oph)} and {path.get_links(path.rev(oph))}. All the
  # links are supposed to have the same move parameters as the rasters.
  # The contact information is obtained from {path.get_contacts} and
  # {contact.get_side_paths(ct)}; the results of these procedures two
  # must be consistent with each other and with {OPHS}.
  #
  # All the traces in {OPHS} must have the same parameters, including
  # the same nominal width that should be equal to the spacing of the
  # scanlines. Otherwise, the {MoveParms} records of the moves in
  # {OCRS}, {OPHS}, and {OLKS} are ignored. So are al the move, path,
  # and contact names, and the cooling time limit {tcool_limit(ct)} of
  # contacts in {CTS}.
  #
  # Before being written to the file, every point will be geometrically
  # transformed by subtracting the 2-vector {shift} from it and then
  # rotating it by {angle} about the origin.
  return txt_write_IMP.write(wr, OCRS, OPHS, Z, angle, shift)
