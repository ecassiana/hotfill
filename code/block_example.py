# Tools to create some example blocks.

import block_example_IMP

# SIMPLE BLOCKS:

def single_raster(xlo, xhi, y, mp_trace):
  # Returns a block that is a single raster line with endpoints
  # {(xlo,y)} and {(xhy, y)} and parameters {mp_trace}, with two choices
  # -- the two orientations.
  # 
  #   k name n paths 
  #   - ---- - --------------
  #   0 Bsr  2 Ph,~Ph
  return block_example_IMP.single_raster(xlo, xhi, y, mp_trace)

def raster_rectangle(plo, nx, ny, hor, ver, alt, mp_trace, mp_jump):
  # Returns a block whose choices (2, 4, or 8) are paths of
  # raster traces filling a rectangle, connected by link traces or jumps.
  # 
  # All moves in the block will share the {Move_Parms} records
  # {mp_trace} or {mp_jump}, as appropriate. The spacing of the rasters
  # will be the trace width {wd} implied by {mp_trace}. The ENDPOINTS of
  # the rasters will span a rectangle of width {(nx-1)*wd} and height
  # {(ny-1)*wd}, with lower left corner {plo}.
  #
  # At least one of the bools {hor} and {ver} must be true, and the
  # integers {nx,xy} must be positive. If the bool {hor} is true, {nx}
  # must be at least 2, and the choices of the block will include two
  # paths that share the same {ny} horizontal rasters of length
  # {(nx-1)*wd}. If {ver} is true, {ny} must be at least 2, and the
  # choices will include two paths that share the same {nx} vertical
  # rasters of length {(ny-1)*wd}. In either casse, the paths differ on
  # the orientation of the first raster, and therefore on the set of
  # connectors used.
  #
  # If the bool {alt} is false, the rasters will have the same
  # orientation and will be connected by jumps. If {alt} is true, they
  # will have alternating orientations and will be connected by link
  # traces.
  #
  # In any case, if there are at least two rasters, the reversal of each
  # path {ph} above will also be included in the block, right after
  # {ph}. Thus, depending on the parameters, the block may have 2, 4, or
  # 8 choices.
  #
  # Note that the nominal sausages and decorations (arrows, dots, etc.)
  # will overflow that rectangle by {wd/2} or more in every direction.
  # The counts {nx} and {ny} must be at least 2.
  return block_example_IMP.raster_rectangle(plo, nx, ny, hor, ver, alt, mp_trace, mp_jump)

def spiral_rectangle(plo, szx, szy, mp_trace):
  # Returns a block with 4 choices, which are spral paths filling
  # the rectangle with corners {plo} and {plo + (szx,szy)},
  # as in {path_example.spiral_rectangle}. All choices will turn
  # counterclockwise, starting at the four corners of the rectangle
  # and ending somewhere in the middle.
  return block_example_IMP.spiral_rectangle(plo, szx, szy, mp_trace)

def onion(nch, ctr, Rc, mp_cont, Rf, mp_fill, mp_jump):
  # Returns a block that has {nch} choices, each being a path created by
  # {path_example.onion} with those parameters, but with {nch} values of {phase}
  # equally spaced on the circle.
  return block_example_IMP.onion(nch, ctr, Rc, mp_cont, Rf, mp_fill, mp_jump)
  
# SAMPLE BLOCKS FOR THE UNIT TEST PROGRAM

def misc_A(mp_trace, mp_jump):
  # Returns a list {BCS} of three blocks whose choices are created by {path_example.misc_E}.
  # Also returns a list {OPHS} of the paths and lists {TRS} and {JMS}
  # of the traces and moves.
  # 
  # Blocks ({n} = number of choices):
  #
  #   k name n paths 
  #   - ---- - --------------
  #   0 BA0  4 Pa,~Pa,Pb,~Pb
  #   1 BA1  2 Pd,Pe
  #   2 BA2  1 ~Pc
  # 
  # Those paths have a total of ten traces and five internal jumps. 
  # The coordinates are independent of the width of the moves.
  # Note that choices 0 and 1 of block {BCS[0]} share the same moves, and
  # ditto for choices 2 and 3. 
  return block_example_IMP.misc_A(mp_trace, mp_jump)

def misc_B(mp_trace, mp_jump):
  # Returns a list {BCS} of three blocks whose choices are created by {path_example.misc_E}.
  # Also returns the list {OPHS} of the paths and the lists {TRS} and {JMS}
  # of their traces and moves.
  #
  # This example is similar to {misc_A} but the paths are joined into blocks in different ways.
  #
  # The blocks and their choices are:
  #
  #   block  n choices
  #   ------ - ----------
  #   BCS[0] 5 OPHS[0], ~OPHS[0], OPHS[4], ~OPHS[4], OPHS[2]
  #   BCS[1] 2 OPHS[1], ~OPHS[1]
  #   BCS[2] 1 OPHS[3]
  #
  # where {n} is the number of choices, and {~ph} means {path.rev(ph)}.
  return block_example_IMP.misc_B(mp_trace, mp_jump)
 
def misc_C(mp_trace):
  # Returns a list {BCS} of seven blocks, each a serpentine raster fill of a rectangle.
  # Five are arranged in a horizontal row, and the other two are
  # above and below the seecond of these five blocks.
  return block_example_IMP.misc_C(mp_trace)
 
def misc_D(mp_trace, contacts):
  # Returns a list {BCS} of blocks and list {CTS} of contacts.
  # 
  # The list {BCS} will have four blocks, each a serpentine raster fill
  # of a rectangle. There are two wide blocks at top and bottom, and two
  # narrower "roads", side by side, between the two. The bottom block
  # and the roads have four choices, all with horizontal traces, while
  # the top block has eight choices, half with horizontal traces and
  # half with vertical traces.
  #
  # If {contacts} is true, creates and returns in the list {CTS} some
  # contacts between the blocks. They are attached to the paths with
  # {path.add_contacts}. Otherwise the list {CTS} is empty.
  #
  return block_example_IMP.misc_D(mp_trace, contacts)
 
def misc_E(mp_trace, mp_jump):
  # Returns a list {BCS} of two blocks, each a raster fill of a
  # rectangle, side by side, each with four choices. The left one has
  # horizontal rasters, in alternating directions, connected by link
  # traces, while the right one has vertical rasters, in the same
  # direction, connected by jumps.
  return block_example_IMP.misc_E(mp_trace, mp_jump)

def misc_G(mp_cont,mp_fill,mp_jump):
  # Returns a list {BCS} of two blocks, a list of the three paths {PHS}
  # used in the blocks (not including their reversals), and two lists {TRS0,TRS1}
  # of the trace {Move} objects used in those paths (without repetitions
  # or reversals). 
  #
  # The paths comprise a few raster moves on three scan lines, with
  # various contacts between them. Block {BCS[0]}, whose traces are in
  # {TRS[0]}, uses paths {PHS[0]} and {PHS[2]} in both orientations.
  # Block {BCS[1]}, whose traces are in {TRS1}, uses path {PHS[1]} in
  # both orientations.
  # 
  # The blocks:
  # 
  #   k name n paths 
  #   - ---- - --------------
  #   0 BG0  4 PG0,~PG0,PG2,~PG2
  #   1 BG1  2 PG1,~PG1
  #
  return block_example_IMP.misc_G(mp_cont,mp_fill,mp_jump)

# TOOLS

def raster_raster_block_contact(bc0, bc1):
  # ???
  return block_example_IMP.raster_raster_block_contact(bc0, bc1)
