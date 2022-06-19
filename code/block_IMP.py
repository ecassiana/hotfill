# Implementation of the module {block}.

import block
import path
import move
import move_parms
import contact
import job_parms
import hacks
import rn 
import pyx
import sys
from math import sqrt, floor, ceil, exp, log, sin, cos, acos, pi, nan,  inf

class Block_IMP:
  # A {Block_IMP} object {bc} has a list {bc.OPHS} of {nch} one or more oriented
  # paths that are the choices.
  
  def __init__(self, OPHS ):
    self.OPHS = OPHS
    self.name = None
    # Fields for the heuristic:
    self.used_choice = None  # Which choice was used in {Q}, or {None}
    self.seams = None        # Seams to which this block is side of.
    self.contacts = None     # List of contacts whose sides are in this block.

def from_paths(OPHS):
  assert type(OPHS) == list or type(OPHS) == tuple
  assert len(OPHS) >  0
  bc = block.Block(tuple(OPHS))
  return bc

def nchoices(bc):
  assert isinstance(bc, block.Block)
  n = len(bc.OPHS)
  return n
  
def choice(bc, ich):
  assert isinstance(bc, block.Block)
  n = len(bc.OPHS)
  assert ich >= 0 and ich < n
  och = bc.OPHS[ich] # Choice of the elements before reversal.
  return och

def avg_choices(BCS):
  slog = 0
  for bc in BCS:
    nch = block.nchoices(bc)
    assert nch >=1
    slog += log(nch)
  avg = exp(slog/len(BCS))
  return avg

def bbox(BCS, links, contacts):
  B = None
  for bc in BCS:
    B = rn.box_join(B, path.bbox(bc.OPHS))
    if links:
      for oph in bc.OPHS:
        OLKS = list(path.get_links(oph)) + list(path.get_links(path.rev(oph)))
        B = rn.box_join(B, path.bbox(OLKS))
    if contacts:
      for oph in bc.OPHS:
        CTS = list(path.get_contacts(oph, 0)) + list(path.get_contacts(oph, 1))
        B = rn.box_join(B, contact.bbox(CTS))
  return B
  # ----------------------------------------------------------------------

def barycenter(bc):
  wpsum = (0,0) # Sum of all points, weighted by move length.
  wsum = 0 # Sum of weights (move length x width).
  nch = block.nchoices(bc)
  for ich in range(nch):
    och = block.choice(bc, ich)
    for imv in range(path.nelems(och)):
      omv = path.elem(och, imv)
      if not move.is_jump(omv):
        p = move.pini(omv)
        q = move.pfin(omv)
        w = rn.dist(p, q)*move.width(omv)
        wpsum = rn.mix3(1,wpsum, w/2,p, w/2,q)
        wsum += w
  assert wsum > 0, "no traces in block"
  pbar = rn.scale(1/wsum, wpsum)
  return pbar
  # ----------------------------------------------------------------------

def moves(bc):
  MVS = []
  for ich in range(nchoices(bc)):
    oph = choice(bc,ich)
    for ich in range(path.nelems(oph)):
      omv = path.elem(oph,ich)
      mv, dr = move.unpack(omv)
      MVS.append(mv)
  MVS = list(set(MVS))
  return MVS
  # ----------------------------------------------------------------------
    
def find_choice_with_move(bc, omv):
  for ich in range(nchoices(bc)):
    oph = choice(bc, ich)
    if path.find_move(oph, omv) != None:
      return ich
  return None
  #----------------------------------------------------------------------

def find_block_with_move(BCS, omv):
  nbc = len(BCS)
  for ich in range(nbc):
    bck = BCS[ich]
    if find_choice_with_move(bck, omv) != None:
      return ich
  return None
  # ----------------------------------------------------------------------

def min_fabtime(bc):
  assert isinstance(bc, block.Block)
  n = len(bc.OPHS)
  mex = +inf
  for kbc in range(n):
    ophk = bc.OPHS[kbc]
    texk = path.fabtime(ophk)
    if texk < mex: mex = texk
  return mex
  # ----------------------------------------------------------------------

def min_tot_fabtime(BCS):
  mex = 0
  for bc in BCS:
    mex += block.min_fabtime(bc)
  return mex
  # ----------------------------------------------------------------------

def plot_to_files(fname, BCS, CLRS, rwd, wd_axes, matter, links, contacts):

  # Find the max number of choices {nch}:
  nch = 0
  for bc in BCS: nch = max(nch, nchoices(bc))
  
  # Get a plotting bounding box of all blocks:
  B = bbox(BCS, True, True)

  c, szx,szy = hacks.make_canvas(hacks.round_box(B,0.5), None, None, True, True, 1, nch)
  ystep = szy
   
  axes = True
  dots = True
  arrows = True
  for ich in range(nch):
    dp = (0, ich*ystep)
    plot_standard \
      ( c, BCS, ich, dp, CLRS=CLRS, rwd=rwd, wd_axes=wd_axes,
        axes=axes, dots=dots, arrows=arrows, matter=matter, links=links, contacts=contacts
      )
  hacks.write_plot(c, fname)
  return
  # ----------------------------------------------------------------------

def plot_standard(c, BCS, ich, dp, CLRS, rwd, wd_axes, axes, dots, arrows, matter, links, contacts):

  assert type(BCS) is list or type(BCS) is tuple
  nbc = len(BCS)
  if CLRS == None: 
    CLRS = [ pyx.color.rgb(0.050, 0.800, 0.000), ] # Default trace color.
  else:
    assert type(CLRS) is list or type(CLRS) is tuple
  nclr = len(CLRS)
  assert nclr == 1 or nclr == nbc
  if dp == None: dp = (0,0)

  if matter:
    # Plot the matter shadows of all block choices and links:
    plot_matter_shadow(c, BCS, dp, links)

  if links:
    # Plot the links of choice {ich} of every block, underneath the choices:
    rwd_lk = 0.75 * rwd
    plot_links_of_choices(c, BCS, ich, dp, rwd_lk, wd_axes)

  # Plot the choice {ich} of every blcos:
  rwd_tr = rwd
  plot_choices(c, BCS, ich, dp, CLRS, rwd_tr, wd_axes, axes, dots, arrows)

  if contacts:
    # Plots the contacts with choice {ich} of evey block:
    plot_contacts_of_choices(c, BCS, ich, dp, wd_axes)
  return
  # ----------------------------------------------------------------------

def plot_choices(c, BCS, ich, dp, CLRS, rwd, wd_axes, axes, dots, arrows):
  # Plots the choice with index {ich} of each block in {BCS},
  # with its respective color from the list {CLRS}. Does not plot the 
  # matter shadows, links, or contacts.

  matter = False

  for ibc in range(len(BCS)):
    bc = BCS[ibc]
    if ich < block.nchoices(bc):
      # sys.stderr.write("  plotting it\n")
      clr_ph = CLRS[0] if len(CLRS) == 1 else CLRS[ibc]
      ophi = block.choice(bc, ich)
      path.plot_standard \
        ( c, [ophi,], dp, None, [clr_ph,], rwd=rwd, wd_axes=wd_axes,
          axes=axes, dots=dots, arrows=arrows, matter=matter
        )
  return
  # ----------------------------------------------------------------------

def plot_links_of_choices(c, BCS, ich, dp, rwd, wd_axes):
  # Plots the links that connect to the endpoints of choice {ich} of
  # each block. Will not plot any decoration (axes, dots, arrows,
  # matter).
  #
  # Some links will be plotted twice, but it seems hard to avoid that
  # and it is invisible anyway.

  clr_lk = pyx.color.rgb( 1.000, 0.050, 0.400 )

  axes = False
  dots = False
  arrows = False
  matter = False
 
  for bc in BCS:
    if ich < block.nchoices(bc):
      ophi = block.choice(bc, ich)
      for oph in ophi, path.rev(ophi):
        OLKS = list(path.get_links(oph))
        path.plot_standard( 
          c, OLKS, dp, None, [clr_lk,], rwd=rwd, wd_axes=wd_axes,
          axes=axes, dots=dots, arrows=arrows, matter=matter )
  return
  # ----------------------------------------------------------------------

def plot_contacts_of_choices(c, BCS, ich, dp, wd_axes):
  # Plots the contacts that have sides on of choice {ich}
  # of each block. Some contacts will be plotted twice, but it seems
  # hard to avoid that and it is invisible anyway.
  
  clr_ct = pyx.color.rgb( 1.000, 0.300, 0.000 )
  
  wd_ct = wd_axes
  sz_tic = 0.0
  arrow = False

  for bc in BCS:
    if ich < block.nchoices(bc):
      oph = block.choice(bc, ich)
      for isd in range(2):
        CTS = path.get_contacts(oph, isd)
        for ct in CTS:
         contact.plot_single(c, ct, dp, clr_ct, None, 0, wd=wd_ct, sz_tic=sz_tic, arrow=arrow)

  return
  # ----------------------------------------------------------------------

def plot_matter_shadow(c, BCS, dp, links):
  rwd_matter = 1.13
  cmatter = hacks.matter_color(None)  # Material footprint color.
  for bc in BCS:
    for ich in range(block.nchoices(bc)):
      ophi = block.choice(bc, ich)
      path.plot_layer \
        ( c, ophi, dp, jmp=False, clr=cmatter, 
          rwd=rwd_matter, wd=0, dashed=False, wd_dots=0, sz_arrows=0 
        )
      if links:
        for olkij in path.get_links(ophi) + path.get_links(path.rev(ophi)):
          path.plot_layer \
            ( c, olkij, dp, jmp=False, clr=cmatter, 
              rwd=rwd_matter, wd = 0, dashed=False, wd_dots=0, sz_arrows=0
            )  
  return
  # ----------------------------------------------------------------------

def validate(bc):
  n = len(bc.OPHS)
  assert n >= 1
  # sys.stderr.write("--- validating bc = %s n = %d ----------\n" % (str(bc), n))
  # Validate the list of choices:
  for och in bc.OPHS:
    path.validate(och)
    ph, dr = path.unpack(och) # Typechecking.
    nmv = path.nelems(ph)
    assert nmv > 0
    assert not move.is_jump(path.elem(ph,0))
    assert not move.is_jump(path.elem(ph,nmv-1))
  return
  # ----------------------------------------------------------------------

def has_name(bc):
  return bc.name != None
  # ----------------------------------------------------------------------

def get_name(bc):
  assert isinstance(bc, block.Block)
  name = bc.name
  if name == None: name = "B?"
  return name
  # ----------------------------------------------------------------------

def set_name(bc, name):
  assert type(name) is str
  bc.name = name
  return
  # ----------------------------------------------------------------------

def tag_names(BCS, tag):
  if tag != None and tag != "":
    assert type(tag) is str
    for bc in BCS:
      bc.name = tag + get_name(bc)
  return
  # ----------------------------------------------------------------------

def show(wr, pref, bc, suff, paths, wna, wnc):
  if pref != None: wr.write(pref)
  wr.write("%-*s" % (wna,get_name(bc)))
  nch = nchoices(bc)
  wr.write(" %*d" % (wnc, nch))
  if paths:
    wr.write(" ")
    for ich in range(nch):
      oph = choice(bc, ich)
      phname = path.get_name(oph)
      if ich > 0: wr.write(",")
      wr.write(phname)
  if suff != None: wr.write(suff)
  return
  # ----------------------------------------------------------------------
  
def show_list(wr, pref, BCS, suff, paths):
  assert type(BCS) is list or type (BCS) is tuple
  nbc = len(BCS)
  if nbc == 0: return
  wna = 4 # Width of "name" column; min 4 because of the header.
  wnc = 1 # Width of "number of paths" column.
  for bc in BCS: 
    wna = max(wna, len(get_name(bc)))
    wnc = max(wnc, len(str(nchoices(bc))))
  wix = len(str(nbc-1)) # Num digits in index.
  
  wr.write("\n")

  # Write header:
  wr.write("%*s%*s %-*s %*s paths \n" % (len(pref),'',wix,"k",wna,"name",wnc,"n"))
  wr.write("%*s%s %s %s --------------\n" % (len(pref),'',"-"*wix,"-"*wna,"-"*wnc))
  
  # Write blocks:
  for kbc in range(len(BCS)):
    bc = BCS[kbc]
    if pref != None: wr.write(pref)
    wr.write("%*d " % (wix,kbc))
    show(wr, None, bc, suff, paths, wna, wnc)
    wr.write("\n")

  wr.write("\n")
  return 
  # ----------------------------------------------------------------------
  
