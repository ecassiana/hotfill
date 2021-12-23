# Tools and types to handle blocks -- precomputed subpaths of a tool-path.
# Last edited on 2021-10-29 14:20:22 by stolfi

import block_IMP; from block_IMP import Block_IMP

class Block(Block_IMP):
  # An object {bc} of class {Block} represents a /block/, a collection
  # of /alternative paths/ or /choices/ that can be used to create some
  # specific of part the slice. If the block {bc} has {n} alternatives,
  # the function {choice(bc,ich)} below returns the alternative path with
  # index {ich} in {0..n-1}.
  #
  # In a simple example, let {P1,...,Pm} be a set {m >= 2} raster lines
  # lying on consecutive scan lines, each making significant contact only
  # with the previous and following ones; except that {P1} may have
  # contacts with zero or more than one rasters below it, and {Pm} may
  # have contacts with zero or more than one above it. We may decide that
  # those rasters are to be fabricated consecutively, in alternating
  # directions. There are four ways to do that: one can start from {P1} or
  # from {Pm}, and the starting raster line may be fabricated from left to
  # right or from right to left.  
  #
  # These four alternative paths are actually two distinct paths {Q0}
  # and {Q1} and their reversals. While all these four paths share the
  # same {Move} objects that represent the raster traces, {Q0} and {Q1}
  # also have different /connectors/ -- extra jumps, traces, or paths
  # that connect one end of one raster line to the startof the next one.
  # Therefore {Q0} and {Q1} will usually have different execution times
  # and different relative covering times for external contacts.
  #
  # By joining those rasters into a block, the user is specifying that
  # exactly one of these four choices -- {Q0}, {Q1}, {rev(Q0)}, or
  # {rev(Q1)} -- should be inserted into the final tool-path, as a a
  # single stretch of consecutive moves.
  #
  # As another example, a component of the slice's contour -- such as a
  # hole or an ``island'' -- can be represented as a block of {m} paths.
  # Each path covers the same trajectory in the same direction, but starts
  # and ends at a different point.
  #
  # Joining several elementary paths of the input data set into blocks
  # greatly simplifies the path optimization problem, by reducing the
  # number of choices from {2^N*N!} for {N} paths to {K^M} for {M} blocks,
  # each with a geometric average of {K} choices (including reversals when
  # applicable). 
  #
  # Joining paths into blocks saves computing time also by making it
  # unnecessary to recompute the cover and cooling times of internal
  # contacts (between traces of the same alternative of a block) since
  # their cooling times are fixed, and thus cooling time limit
  # constraints can be checked as the block is created.
  #
  # On the other hand, by limiting the choices of the heuristic, this
  # optimization may prevent it from finding the optimal path. It may even
  # make the problem unsolvable, if a block is so large that doing it all
  # in one go would inevitably violate some cooling constraints elsewhere.
  # Still, if the blocks are chosen judiciously, this negative aspect may
  # be minimized. In fact, this optimization may result in a better final
  # path, because it would prevent the heuristic from finding some expensive
  # paths and stopping there.
  #
  # ASSEMBLING THE BLOCKS
  #
  # In general, one should consider joining several elements of the
  # original problem into a block only if they are close to each other
  # and mostly far from the rest of the slice; so that a good path is
  # very likely to fabricate them consecutively, once it starts with any
  # of them.
  #
  # For simple blocks -- like sequences of adjacent raster lines, or the
  # traces that make up a closed contour line -- the choices that are
  # worth considering can be identified and constructed by ad-hoc
  # procedures. In more complex situations, one may have to use the
  # heuristic itself to choose the alternative paths.
  #
  # For instance, if the block is to consist of the contours of six
  # small islands arranged in a circle, somewhat distant from the rest
  # of the slice, it may be desirable to precompute some or all of the
  # {6*5/2 = 15} optimal tool-paths that trace the six contours in one
  # go, begining and ending at the outermost point of each pair of
  # contours. That pre-processing may still be better than letting each
  # island be a separate input block, and letting the heuristic
  # discover by trial and error that those contours are indeed better
  # fabricated consecutively and *repeadly* try to find the optimal order
  # among all the {6! = 120} possibilities.
  #
  # A {Block} object also has a mutable /name/ attribute, that is used
  # for debugging and documentation. It must be either a string or {None}.
  #
  pass

def from_paths(OPHS):
  # Creates and returns a {Block} object {bc} whose alternatives are the 
  # oriented paths in the list or tuple {OPHS}.  Specifically, {choice(bc,ich)}
  # will be {OPHS[ich]} for {ich=0,1,...,n-1} where {n=len(OPHS)}.
  #
  # The list must have at least one path. It does not make sense for a
  # choice to be a trivial path (with zero moves), or to begin or end
  # with a jump.
  #
  # Note that if a path {oph} its to be considered in its two possible
  # directions, both {oph} and {path.rev(oph)} must be included in
  # the list {OPHS}.  This only makes sense if the path's initial
  # and final point are well-spaced, or there will be cooling constraints
  # between some of its traces and other blocks.
  return block_IMP.from_paths(OPHS)
  
def nchoices(bc):
  # Returns the number of choices in the {Block} object {bc}.
  return block_IMP.nchoices(bc)
  
def choice(bc, ich):
  # Returns the oriented path which is the alternative of index {ich} 
  # of the block {bc}. The index {ich} must be in
  # {0..nch-1} where {nch=nchoices(bc)}.
  return block_IMP.choice(bc, ich) 
  
def avg_choices(BCS):
  # Returns the geometric average of the number of choices 
  # of the blocks in the list {BCS}.
  return block_IMP.avg_choices(BCS)
  
def bbox(BCS, links, contacts):
  # The parameter {BCS} should be a list or tuple of {Block} objects.
  # Returns the bounding box of all choices of all those blocks,
  # as a pair of points {(plo, phi)} that its lower left and upper right
  # corners. If the list {BCS} is empty, return {None}.
  #
  # The bounding box is the smallest axis-aligned rectangle that
  # contains all the endpoints of all moves of all choices of those
  # blocks. 
  #
  # If {links} is true, the box will include also all endpoints of all
  # moves of all link paths associated to those choices and their
  # reversals by {path.get_links}. If {contacts} is true, also includes
  # the endpoints of all contacts associated to those choices by
  # {path.get_contacts}.
  #
  # Note that the box may not contain decoration items such as move
  # axes, dots, and arrowheads, as well as the sausages of traces and
  # tics and arrows of contacts.
  #
  return block_IMP.bbox(BCS, links, contacts)

def barycenter(bc):
  # Returns the barycenter of the block, defined as the barycenter of
  # all its traces (excluding jumps).
  #
  # Each trace is assumed to be a rod of uniform linear density proportional 
  # its width.
  return block_IMP.barycenter(bc)
  
def moves(bc):
  # Returns a list of all the {Move} objects that occur in all the choices
  # of the {Block} {bc}, without orientation bits and without repetitions.
  return block_IMP.moves(bc)

def find_choice_with_move(bc, omv):
  # Given a {Block} object {bc} and an oriented move {omv}, returns the
  # index of the first choice of {bc} contains that move or its reverse.
  # If {bc} has no such choice, returns {None}.
  return block_IMP.find_choice_with_move(bc, omv)

def find_block_with_move(BCS, omv):
  # Given a list {BCS} of {Block} objects and an oriented move {omv},
  # returns the index in {BCS} of the first block that has some choice
  # that contains that move or its reverse. If there is no such block in
  # {BCS}, returns {None}.
  return block_IMP.find_block_with_move(BCS, omv)

def min_fabtime(bc):
  # Returns the minimum of {path.fabtime(oph)} for all choices {oph} of 
  # block {bc}.
  return block_IMP.min_fabtime(bc)
 
def min_tot_fabtime(BCS):
  # Given the list{BCS} of blocks, returns an estimate of the minimum time needed to
  # fabricate them, even if the best alternative could be chosen in each
  # block.
  return block_IMP.min_tot_fabtime(BCS)
  
# PLOTTING

def plot_to_files(fname, BCS, CLRS, rwd, wd_axes, matter, links, contacts):
  # Plots the blocks of {BCS} into a file called {fname}.
  #
  # The argument {BCS} must be a list or tuple of {Block} objects. The
  # {CLRS} argument must be a list or tuple of {pyx.color.rgb} objects.
  #
  # The plot will be a stack of {nch} sub-plots where {nch} is the
  # maximum of {block.nchoices(bc)} for all {bc} in {BCS}.
  # The sub-plot of index {ich} will show choice {ich} of every block, 
  # where that choice exists.  The bottom sub-plot (index 0) will be displaced by 
  # the 2-vector {dp}; if {dp} is {None}, the proceedure assumes {(0,0)}
  # (no displacement). Each subsequent sub-plot will be
  # displaced vertically by an additional {ystep} millimeters.
  #
  # Each sub-plot is plotted with with {plot_standard} below, with the same
  # {CLRS}, {wd_axes}, and {matter} parameters and some default choices for
  # the style parameters ({axes}, {dots}, {arrows}, {links}, and {contacts}).  A millimeter
  # grid will be painted in the background.
  #
  # Then the procedure writes the plot (with all sub-plots) to files
  # "{fname}.{ext}" where {ext} is "eps", "png", and "jpg".
  block_IMP.plot_to_files(fname, BCS, CLRS, rwd, wd_axes, matter, links, contacts)
 
def plot_standard(c, BCS, ich, dp, CLRS, rwd, wd_axes, axes, dots, arrows, matter, links, contacts):
  # The arguments {BCS} and {CLRS} must be lists or tuples of {Block} objexts
  # and {pyx.color.rgb} objects, respectively.
  # 
  # Plots onto the {pyx} canvas {c} the choice {ich} of all blocks of
  # {BCS}. The plot will be displaced by the 2-vector {dp}; if {dp} is
  # {None}, the proceedure assumes {(0,0)} (no displacement).
  #
  # If {links} is true, also plots the link paths of the choices, 
  # obtained {path.get_links}, with a fixed color.
  #
  # If {contacts} is true, also plots the contacts whose sides are on each choice, 
  # obtained {path.get_contacts}, with a fixed color.
  #
  # If a block does not have a choice with index {ich}, it will be omitted
  # from the plot. However, if {matter} is true, the block's /matter
  # shadow/ -- the combined extent of material of all choices -- 
  # will still be painted.
  #
  # If {CLRS} has a single element, uses that color for all trace
  # sausages. If {CLRS} is {None}, uses some default color Otherwise
  # {CLRS} must have the same length as {BCS}, and all traces of block
  # {BCS[k]} will be painted wit color {CLRS[k]}. Each trace sausage is
  # plotted with width {rwd*wd} where {wd} is the trace's nominal width.
  #
  block_IMP.plot_standard(c, BCS, ich, dp, CLRS, rwd, wd_axes, axes, dots, arrows, matter, links, contacts)
  
def plot_matter_shadow(c, BCS, dp, links):
  # Plots to the {pyx} canvas {c} the matter shadow of every block {bc}
  # in the list {BCS}, namely the union of the estimated matter extent
  # of all traces in all choices.
  #
  # If {links} is true, also plots the matter shadow of the links of all choices.
  block_IMP.plot_matter_shadow(c, BCS, dp, links)

# DEBUGGING

def validate(bc):
  # Runs some consistency checks on the {Block} object {bc}.  Aborts with assert failure
  # if it detects any errors.
  block_IMP.validate(bc)

def has_name(bc):
  # Returns {True} if and only {bc}'s name atribute is not {None}.
  return block_IMP.has_name(bc)

def get_name(bc):
  # Given a {Block} object {bc}, returns its name attribute, if not {None}. 
  # If that attribute is {None}, returns "C?" instead.
  return block_IMP.get_name(bc)

def set_name(bc, name):
  # Saves the string {name} as the name attrbute of the {Block} object {bc}.
  block_IMP.set_name(bc, name)

def tag_names(BCS, tag):
  # Prepends the string {tag} (if not {None}) to the names of all
  # {Block} objects in {BCS}. The objects had better be all distinct.
  block_IMP.tag_names(BCS, tag)

def show(wr, pref, bc, suff, paths, wna, wnc):
  # Writes to file {wr} a readable description of the block {bc}.
  #
  # The description includes:
  #
  #   * the name of the block, as returned by {get_name}
  #   * the number of choices {nch}
  #
  # If {paths} is true, the description also includes the list of the
  # names of the paths that comprise the block {bc}; namely, the
  # sequence {path.get_name(block.choice(bc,ich))} for {ich} in
  # {range(nch)}.
  # 
  # The name is padded to {wna} columns and left-aligned. The number of
  # choices is padded to {wnc} columns. The whole description uses one
  # line, is preceded by {pref} (if not {None}, and followed by {suff}
  # (if not {None}). It is NOT terminated by a newline, unless provided
  # by {suff}.
  block_IMP.show(wr, pref, bc, suff, paths, wna, wnc)

def show_list(wr, pref, BCS, suff, paths):
  # Writes to file {wr} a readable description of the {Block} objects in
  # the list {BCS}.
  #
  # The output consists of a header and then the description of each
  # block, one per line, as produced with {show} with the given
  # parameters {paths}. Each line is prefixed by {pref} and the block
  # index in {BCS}, and followed by {suff} and a newline.
  block_IMP.show_list(wr, pref, BCS, suff, paths)
