#! /usr/bin/python3
# Implementation of module {paper_figures_B}
# Last edited on 2021-11-25 22:02:53 by stolfi

import paper_figures_B
import paper_example_B
import move
import hotfill
import move_example
import move_parms
import contact
import path
import path_example
import rootray_shape
import txt_write
import rootray
import raster
import block
import hacks

import rn

import pyx
import random
import sys
from math import sqrt, hypot, sin, cos, atan2, floor, ceil, inf, nan, pi

def plot_figure(fname, fig, subfig):
  
  title = "### figure %s sub-figure %s ###" % (fig,subfig)
  sys.stderr.write("#"*80 + "\n")
  sys.stderr.write(title + "#"*(80-len(title)) + "\n")
  sys.stderr.write("#"*80 + "\n")
  sys.stderr.write("\n")
  
  style = make_style_dict()
  color = make_color_dict()
  
  parms_fname = "tests/out/paper_figures_B_TST_p_fig_printer_parms.tex"
  write_text_parms(parms_fname, style)

  txt_fname = "tests/out/paper_figures_B_TST_p_turkey.txt"
  write_txt_file(txt_fname, style)
  
  if fig == "paths":
    c = plot_figure_paths(subfig, style,color)
  elif fig == "moves":
    c = plot_figure_moves(subfig, style,color)
  elif fig == "cover":
    c = plot_figure_cover(subfig, style,color)
  elif fig == "input":
    c = plot_figure_input(subfig, style,color)
  elif fig == "zigzag":
     c = plot_figure_zigzag(subfig, style,color)
  elif fig == "cold":
     c = plot_figure_cold(subfig, style,color)
  elif fig == "rivers":
     c = plot_figure_rivers(subfig, style,color)
  elif fig == "blocks":
     c = plot_figure_blocks(subfig, style,color)
  elif fig == "hotfill":
     c = plot_figure_hotfill(subfig, style,color)
  else:
     assert False, ("figure '%s ' not implemented.\n" % fig)

  if c != None:
    hacks.write_plot(c, fname)
  else:
    sys.stderr.write("!! figure '%s subfigure %s' not plotted.\n" % (fig,subfig))

  sys.stderr.write("#"*80 + "\n")
  sys.stderr.write("\n")

  return
  # ----------------------------------------------------------------------

def write_txt_file(fname, style):

  # Get the contours, fillers, links, contacts, and graph of the "turkey" part:
  mp_cont, mp_fill, mp_link, mp_jump = make_move_parms(style);
  OCRS,OPHS,OLKS,CTS,VGS,EGS = paper_example_B.make_turkey(mp_cont, mp_fill,mp_link,mp_jump) 

  wr = open(fname, "w")
  Z = 10.0
  angle = 0
  shift = (0,0)
  for oph in OPHS: path.set_group(oph,0)
  txt_write.write(wr, OCRS, OPHS, Z, angle, shift)
  wr.close()

  # ----------------------------------------------------------------------

def plot_figure_paths(subfig, style,color):
  # Plots the figure that shows a smple path {OPHS[0]} with labels "pini", etc.
  # and its reverse.
  
  assert subfig == "basic", ("invalid subfig %s" % subfig)
  
  # Get the paths:
  mp_cont, mp_fill, mp_link, mp_jump = make_move_parms(style);
  OPHS = paper_example_B.make_simple_path(mp_fill, mp_jump)
  
  assert len(OPHS) == 1
  oph = OPHS[0]
  nmv = path.nelems(oph) 

  B = path.bbox([oph,]) # Bounding box of all move endpoints.

  Xstep = B[1][0] - B[0][0] + 8.0  # Displacement betwwen the two versions of the path.
  
  # Compute the figure's bounding box {Bfig}:
  Bfig = (B[0], rn.add(B[1], (Xstep,0))) # Bounding box for both versions of the path.

  # Add space at bottom and top for labels:
  mrg0 = (3.5, 2.0) # Extra margin for labels at left and bottom
  mrg1 = (3.5, 2.4) # Extra margin for labels at right and top.
  Bfig = rn.box_expand(Bfig, mrg0, mrg1)

  # # Widen {Bfig} symmetrically to standard "math figure" widtdh:
  # Bfig = widen_box_for_math_figure(Bfig, style)
 
  autoscale = False
  c = make_figure_canvas(Bfig, autoscale, style,color)

  for rev in False, True:
    dp = (int(rev)*Xstep, 0)
    
    ophr = path.rev(oph) if rev else oph

    # Plot the matter shadow:
    plot_trace_matter(c, [ophr,], dp, style,color)

    # Plot the path:
    rwdf = style['rwd_fill']
    clr = color['fill'] 
    deco = True
    plot_paths(c, [ophr,], dp, [clr,], rwdf,deco, style,color)
    
    # Plot the endpoints of the paths in black and larger:
    plot_path_endpoints(c, [ophr,], dp, style,color)
    
    # Labeled points ({xxxP}) and label displacements ({xxxD}) on the unreversed path:
    piniP = path.pini(oph); piniD = (-1.2, +1.2)
    pfinP = path.pfin(oph); pfinD = (-1.0, +1.2)
    pmovP = [ rn.mix(0.5,move.pini(mvk),0.5,move.pfin(mvk)) for k in range(nmv) for mvk in (path.elem(oph,k),) ]
    pmovD = [ 
      (-3.0, -0.9),                # P[0]
      (-1.0, -1.8 -0.3*int(rev)),  # P[1]
      (-1.6, -1.8 -0.3*int(rev)),  # P[2]
      (+0.8 + 0.2*int(rev), -0.9), # P[3]
    ]
    pmidP = rn.mix(0.5,piniP, 0.5,pfinP); pmidD = (-0.3, +1.0)

    txphr = r"P" if not rev else r"{\mkern-0.7\thinmuskip\overleftarrow{P}\mkern-0.7\thinmuskip}"
    
    txpini = r"$\mathop{p_{\su{INI}}}(%s)$" % txphr
    txpfin = r"$\mathop{p_{\su{FIN}}}(%s)$" % txphr
    txpmid = r"$%s$" % txphr
    txmovS = [ (r"$%s[%d]$" % (txphr,k)) for k in range(nmv) ]

    mfsize = style['mathfsize']
    plot_math_label(c, txpmid, pmidP, pmidD, dp, mfsize, style)
    if not rev:
      plot_math_label(c, txpini, piniP, piniD, dp, mfsize, style)
      plot_math_label(c, txpfin, pfinP, pfinD, dp, mfsize, style)
      for k in range(nmv):
        plot_math_label(c, txmovS[k], pmovP[k], pmovD[k], dp, mfsize, style)
    else:
      plot_math_label(c, txpfin, piniP, piniD, dp, mfsize, style)
      plot_math_label(c, txpini, pfinP, pfinD, dp, mfsize, style)
      for k in range(nmv):
        plot_math_label(c, txmovS[nmv-1-k], pmovP[k], pmovD[k], dp, mfsize, style)
  
  return c
  # ----------------------------------------------------------------------

def plot_figure_moves(subfig, style,color):
  # If {subfig} == "dir" shows a simple move {OPHS[0]} and a simple jump {OPHS[1]}, with labels "pini", etc.
  # If {subfig} == "rev" plots the reversed moves instead.
  
  assert subfig == "dir" or subfig == "rev", ("invalid subfig %s" % subfig)
  
  # Get the two single-move paths:
  mp_cont, mp_fill, mp_link, mp_jump = make_move_parms(style);
  OPHS = paper_example_B.make_simple_moves(mp_fill, mp_jump)

  assert len(OPHS) == 2

  B = path.bbox(OPHS) # Bounding box of endpoints of both moves (undisplaced).
  Xstep = B[1][0] - B[0][0] + 11.0  # Displacement betwwen the two versions of the path.
  
  # Compute the figure's bounding box {Bfig}:
  Bfig = (B[0], rn.add(B[1], (Xstep,0))) # Bounding box for both versions of the path.
  
  # Add space for labels at top and bottom:
  mrg0 = (5.0,1.5) # Extra margin for labels at left and bottom
  mrg1 = (5.0,1.5) # Extra margin for labels at right and top.
  Bfig = rn.box_expand(Bfig, mrg0, mrg1)
  
  # # Widen to standars "math figure" width:
  # Bfig = widen_box_for_math_figure(Bfig, style)
  
  autoscale = False
  c = make_figure_canvas(Bfig, autoscale, style,color)
    
  rev = (subfig == "rev")

  for iph in range(2):
    dp = (iph*Xstep, 0)
    
    oph = OPHS[iph]
    assert path.nelems(oph) == 1
    omv = path.elem(oph,0)
    
    ophr = path.rev(oph) if rev else oph
    omvr = move.rev(omv) if rev else omv

    # Plot the matter shadow:
    plot_trace_matter(c, [oph,], dp, style,color)
    
    # Plot the move:
    clr = color['fill'] 
    rwdf = style['rwd_fill']
    deco = True
    plot_paths(c, [ophr,], dp, [clr,], rwdf,deco, style,color)
    
    # Labeled points ({xxxP}) and label displacements ({xxxD}) on the unreversed moves:
    piniP = move.pini(omv);  
    pfinP = move.pfin(omv);
    pmidP = rn.mix(0.5,piniP, 0.5,pfinP); 
    u, ulen = rn.dir(rn.sub(pfinP, piniP))
    v = ( -u[1], +u[0] )
    
    piniD = rn.mix( -4.1 + 0.4*iph - 0.2*int(rev), u, +0.1,           v )
    pfinD = rn.mix( +1.0 - 0.4*iph,                u, +0.1,           v )
    pmidD = rn.mix( -0.5,                          u, +1.0 - 0.3*iph, v )

    txmvr = r"{\mkern-0.7\thinmuskip\overleftarrow{m}\mkern-0.7\thinmuskip}" if rev else r"m"
    
    txpini = r"$\mathop{p_{\su{INI}}}(%s)$" % txmvr
    txpfin = r"$\mathop{p_{\su{FIN}}}(%s)$" % txmvr
    txpmid = r"$%s$" % txmvr

    mfsize = style['mathfsize']
    plot_math_label(c, txpmid, pmidP, pmidD, dp, mfsize, style)
    if not rev:
      plot_math_label(c, txpini, piniP, piniD, dp, mfsize, style)
      plot_math_label(c, txpfin, pfinP, pfinD, dp, mfsize, style)
    else:                                          
      plot_math_label(c, txpfin, piniP, piniD, dp, mfsize, style)
      plot_math_label(c, txpini, pfinP, pfinD, dp, mfsize, style)
  
  return c
  # ----------------------------------------------------------------------

def plot_figure_cover(subfig, style,color):
  # Plots the figure that shows a path covering one contact and closing another.
  
  assert subfig == "basic", ("invalid subfig %s" % subfig)

  # Get the two paths {P,Q}, and the contacts {ctA,ctB}, where {ctA} is covered by {P} and {Q},
  # {ctB} is closed by {P}:
  mp_cont, mp_fill, mp_link, mp_jump = make_move_parms(style);
  OPHS, CTS = paper_example_B.make_simple_cover(mp_fill, mp_jump)
  
  nph = len(OPHS); assert nph == 2
  nct = len(CTS); assert nct == 2
  
  P = OPHS[0] # Main path.
  Q = OPHS[1] # Secondary path.

  B = path.bbox(OPHS) # Bounding box of all move endpoints.
  B = rn.box_join(B, contact.bbox(CTS)) # Just paranoia.
  
  # Compute the figure's bounding box {Bfig}:
  Bfig = rn.box_expand(B, (2.2,2.2), (1.2,1.2))

  # # Widen {Bfig} symmetrically to standard "math figure" widtdh:
  # Bfig = widen_box_for_math_figure(Bfig, style)
 
  autoscale = False
  c = make_figure_canvas(Bfig, autoscale, style,color)

  dp = (0, 0)

  # Plot the paths:
  rwdf = style['rwd_fill']
  clrP = color['fill'] 
  clrQ = ghostify_colors([color['hifi'],], color['ghost'])[0] 
  deco = True
  plot_paths(c, OPHS, dp, [clrP,clrQ], rwdf,deco, style,color)
  
  # Label the paths and the endpoints:
  mfsize = style['mathfsize'] # Font size for math labels.
  pfsize = style['detafsize'] # Font size for points.

  ptPe0 = path.pini(P)
  plot_dot_on_path(c, ptPe0, dp, style,color)
  Pe0lab_tx = "A"
  Pe0lab_dp = (-1.5, -0.5)
  plot_point_label(c, Pe0lab_tx, ptPe0, Pe0lab_dp, dp, pfsize, style,color)

  ptPe1 = path.pfin(P)
  plot_dot_on_path(c, ptPe1, dp, style,color)
  Pe1lab_tx = "E"
  Pe1lab_dp = (-0.5, +0.7)
  plot_point_label(c, Pe1lab_tx, ptPe1, Pe1lab_dp, dp, pfsize, style,color)
  
  ptphP = rn.mix(0.5,ptPe0, 0.5,ptPe1)
  phPlab_tx = r"$P$"
  phPlab_dp = (+2.0,-3.8)
  plot_math_label(c, phPlab_tx, ptphP, phPlab_dp, dp, mfsize, style)
 
  ptphQ = path.pini(Q)
  phQlab_tx = r"$Q$"
  phQlab_dp = (-1.2,-2.0)
  plot_math_label(c, phQlab_tx, ptphQ, phQlab_dp, dp, mfsize, style)

  # Plot the contacts:
  interbc = False
  shadow = True
  trim = False
  plot_contacts(c, CTS, dp, interbc, shadow,trim, style,color)
  
  # Label points on the path {P} where it covers the midpoint of each contact:
  for ict in range(nct):
    ct = CTS[ict]
    ctends = contact.endpoints(ct)
    pmd = rn.mix(0.5,ctends[0], 0.5,ctends[1]) # Midpoint of contact.
    for jsd in range(2):
      mv = contact.side_move(ct, jsd)
      # Find the move {mv} in either {P} or {Q}:
      for kph in range(2):
        oph = OPHS[kph]
        imv = path.find_move(oph, mv)
        if imv != None: break
      assert imv != None
      omv = path.elem(oph, imv) # Oriented version of {mv}.

      # Find point {pt} on {mv} closest to midpoint of {ct}:
      r = rn.pos_on_line(move.pini(omv), move.pfin(omv), pmd)
      assert 0 < r and r < 1, "move does not cover midpoint of contact"
      pt = rn.mix(1-r,move.pini(omv), r,move.pfin(omv))

      # Label points on the path, moves, and contacts:
      du = (+0.5,+0.2) # Displacement along contacts.
      dv = (-0.2,+0.5) # Displacement perpendicular to contacts.
      if ict == 0 and jsd == 0:
        plab_tx = "B";        plab_dp = rn.mix(1.0,(+0.2, -1.5), +0.3,dv)
        mvlab_tx = r"$m_1'$";  mvlab_dp = rn.mix(1.0,plab_dp, -3.3,du)
        ctlab_tx = r"$c_1$";   ctlab_dp = rn.mix3(1.0,(0.0,0.0), -5.6,du, +2.0,dv)
      elif ict == 0 and jsd == 1:
        plab_tx = None;       plab_dp = rn.mix(1.0,(-0.7, +0.7), +0.2,dv)
        mvlab_tx = r"$m_1''$"; mvlab_dp = rn.mix3(1.0,plab_dp, +0.2,du, +0.6,dv)
        ctlab_tx = None
      elif ict == 1 and jsd == 0:
        plab_tx = "C";        plab_dp = rn.mix(1.0,(+0.2, -1.5), +0.3,dv)
        mvlab_tx = r"$m_2'$";  mvlab_dp = rn.mix(1.0,plab_dp, -4.5,du)
        ctlab_tx = r"$c_2$";   ctlab_dp = rn.mix3(1.0,(0.0,0.0), -5.3,du, +2.0,dv)
      elif ict == 1 and jsd == 1:
        plab_tx = "D";        plab_dp = rn.mix(1.0,(-0.7, +0.7), +0.2,dv)
        mvlab_tx = r"$m_2''$"; mvlab_dp = rn.mix3(1.0,plab_dp, +3.9,du, +0.6,dv)
        ctlab_tx = None

      if plab_tx != None:
        plot_point_label(c, plab_tx, pt, plab_dp, dp, pfsize, style,color)
        plot_dot_on_path(c, pt, dp, style,color)
      if mvlab_tx != None:
        plot_math_label(c, mvlab_tx, pt, mvlab_dp, dp, mfsize, style)
      if ctlab_tx != None:
        ctlab_txx = r"\textcolor{red}{%s}" % ctlab_tx
        plot_math_label(c, ctlab_txx, pt, ctlab_dp, dp, mfsize, style)
  
  return c
  # ----------------------------------------------------------------------

def plot_figure_input(subfig, style,color):
  # Plots the figures with individual rasters, as slected by {subfig}: 
  #
  #   subfig = "rasct" rasters and contacts.
  #   subfig = "links" links.
  #   subfig = "graph" contact graph.
  #
  # Returns the canvas {c} with the figure.
  
  assert subfig == "rasct" or subfig == "links" or subfig == "graph", ("invalid subfig %s" % subfig)

  # Get the contours, fillers, links, contacts, and graph of the "turkey" part:
  mp_cont, mp_fill, mp_link, mp_jump = make_move_parms(style);
  OCRS,OPHS,OLKS,CTS,VGS,EGS = paper_example_B.make_turkey(mp_cont, mp_fill,mp_link,mp_jump) 
  
  assert len(OPHS) > 0
  
  Bfig, cuxlo, cuxhi = get_turkey_figure_bbox(OCRS, OPHS, OLKS, style)
  autoscale = True
  c = make_figure_canvas(Bfig, autoscale, style,color)

  dp = (0,0)
  
  if subfig == "rasct":

    # RASTERS AND CONTACTS
    
    contact.show_times(sys.stderr, OPHS, CTS)

    plot_trace_matter(c, OCRS + OPHS + OLKS, dp, style,color)

    plot_contours(c, OCRS, dp, style, color['cont'])

    # Plot the filling path(s) with strong color:
    clr = color['fill']
    rwdf = style['rwd_fill']
    deco = False
    plot_paths(c, OPHS, dp, [clr,], rwdf,deco, style,color)

    # Plot contacts.
    interbc = False
    shadow = False
    trim = False
    plot_contacts(c, CTS, dp, interbc, shadow,trim, style,color)
    
    # Write TeX file with fabtimes:
    write_tex_params("input", subfig, OPHS, None)

  elif subfig == "links" and len(OLKS) > 0:

    # LINKS

    plot_contours(c, OCRS, dp, style, color['ghost'])

    # Plot the filling path(s) with weak color:
    clr = color['ghost']
    rwdf = style['rwd_fill']
    deco = False
    plot_paths(c, OPHS, dp, [clr,], rwdf,deco, style,color)

    # Plot individual links with various colors: 

    xdir = (1,0)
    ydir = (0,1)
    ystep,yphase = raster.get_spacing_and_phase(OPHS, xdir, ydir)
    bicolor = True
    plot_links(c, OLKS, dp, ystep,yphase, bicolor, style,color) 

    # Plot dots at start and end of links: 
    plot_path_endpoints(c, OLKS, dp, style,color)

    # Write TeX file with fabtimes:
    write_tex_params("input", subfig, OLKS, None)

  elif subfig == "graph":

    # CONTACT GRAPH

    plot_contours(c, OCRS, dp, style, color['ghost'])

    # Plot the filling path(s) with weak color:
    clr = color['ghost']
    rwdf = style['rwd_fill']
    deco = False
    plot_paths(c, OPHS, dp, [clr,], rwdf,deco, style,color)

    plot_contact_graph(c, VGS, EGS, dp, style,color)
    
  else:
    assert False, ("invalid subfig %s" % subfig)

  return c
  # ----------------------------------------------------------------------

def plot_figure_zigzag(subfig, style,color):
  # Plots the figures with scanline tool-path using the rasters {OPHS}: 
  #
  #   subfig = "nonalt" non-alternating
  #   subfig = "alt" alternating
  #
  # Returns the canvas {c} with the figure.

  assert subfig == "nonalt" or subfig == "alt", ("invalid subfig %s" % subfig)

  # Get the contours, fillers, links, contacts, and graph of the "turkey" part:
  mp_cont, mp_fill, mp_link, mp_jump = make_move_parms(style);
  OCRS,OPHS,OLKS,CTS,VGS,EGS = paper_example_B.make_turkey(mp_cont, mp_fill,mp_link,mp_jump) 
  
  alt = (subfig == "alt") # True for alternating path.
  
  dp = (0,0)
  
  Bfig, cuxlo, cuxhi = get_turkey_figure_bbox(OCRS, OPHS, OLKS, style)
  autoscale = True
  c = make_figure_canvas(Bfig, autoscale, style,color)

  if alt:
    # Separate the rasters by scan-line:
    xdir = (1,0)
    ydir = (0,1)
    ystep, yphase = raster.get_spacing_and_phase(OPHS, xdir, ydir)
    SCS = raster.separate_by_scanline(OPHS, xdir, ydir, ystep, yphase)

    # Reverse order and orientation of raster elements of {OPHS} on alternate scanlines:
    PFSnew = [] # Raster paths, reversed if needed.
    rev = False # True if next scan-line is to be reversed.
    for Si in SCS:
      # Now {Si} is a list of indices of rasters in one scanline.
      if rev:
        PFSi = [ path.rev(OPHS[jrs]) for jrs in Si ]
        PFSi.reverse()
      else:
        PFSi = [ OPHS[jrs] for jrs in Si ]
      rev = not rev
      PFSnew += PFSi
    OPHS = PFSnew

  # Assemble the path:
  use_links = alt
  ph = path.concat(OPHS, use_links, mp_jump)
  
  plot_contours(c, OCRS, dp, style, color['ghost'])
  
  # Plot path:
  clr = color['fill']
  rwdf = style['rwd_fill']
  deco = False
  plot_paths(c, [ph,], dp, [clr,], rwdf,deco, style,color) # ??? Should plot only arrows without axes ???

  # Plot the path endpoints:
  plot_path_endpoints(c, [ph,], dp, style,color)

  # Find the coldest contact and its cooling time:
  ncold = 1 # Number of coldest contacts to show.
  CTScold = contact.coldest(ph, CTS, ncold)
  contact.show_times(sys.stderr, [ph,], CTScold)
  
  # Write TeX file with fabtimes and cooling times:
  write_tex_params("zigzag", subfig, [ph,], CTScold)

  # Plot coldest contact(s).
  interbc = False
  shadow = True
  trim = False
  plot_contacts(c, CTScold, dp, interbc, shadow,trim, style,color)
    
  # Label path endpoints:
  tdA = (-1.8, -1.6)
  tdB = (-1.6, +0.5) if alt else (+0.8, +0.8)
  pfsize = style['fullfsize']
  plot_path_endpoint_labels(c, ph, "A", tdA, "B", tdB, dp, pfsize, style,color)

  return c
  # ----------------------------------------------------------------------

def get_cold_paths(OPHS):
  # Retruns a list {CRSS} that describes a typical tool-path produced by
  # {RP3} or {slic3r}. Each element of the list {CRSS} specifies the
  # rasters that comprise a snake sub-path ("continuous raster sequence)
  # of that tool-path. It is a list of pairs {(isc,jrs)}, meaning raster {jrs}
  # of scanline {isc}.
  
  CRSS = ( 
    ( (0,0), (1,0), (2,0), (3,0), (4,0), (5,0), (6,0), (7,0), (8,0), 
      (9,0), (10,0), (11,0), (12,0), (13,0), (14,0),
      (15,0), (16,0), (17,0),
    ),
    ( (14,1), (13,1), (12,1), (11,1), (10,1), (9,1), (8,1), (7,1), (6,1), (5,1), (4,1), (3,1), ),
    ( (0,1), (1,1), (2,1), (3,3), (4,3), (5,3), (6,3), (7,2), (8,2), 
      (9,3), (10,3), (11,3), (12,3), (13,3), (14,3), (15,1),
    ),
    ( (14,2), (13,2), ),
    ( (12,2), (11,2), (10,2), (9,2), ),
    ( (6,2), (5,2), (4,2), (3,2), ),
  )
  return CRSS
  # ----------------------------------------------------------------------
  
def plot_figure_cold(subfig, style,color):
  # Plots the figures with scanline path, alternating or not according to {alt}: 
  #
  #   subfig = "path" the path, colored by CRS, with coldest contact(s)
  #
  # Returns the canvas {c} with the figure.

  assert subfig == "path", ("invalid subfig %s" % subfig)
  
  # Get the contours, fillers, links, contacts, and graph of the "turkey" part:
  mp_cont, mp_fill, mp_link, mp_jump = make_move_parms(style);
  OCRS,OPHS,OLKS,CTS,VGS,EGS = paper_example_B.make_turkey(mp_cont, mp_fill,mp_link,mp_jump) 
  
  dp = (0,0)
  
  Bfig, cuxlo, cuxhi = get_turkey_figure_bbox(OCRS, OPHS, OLKS, style)
  autoscale = True
  c = make_figure_canvas(Bfig, autoscale, style,color)

  # Separate the rasters by scan-line:
  xdir = (1,0)
  ydir = (0,1)
  ystep, yphase = raster.get_spacing_and_phase(OPHS, xdir, ydir)
  SCS = raster.separate_by_scanline(OPHS, xdir, ydir, ystep, yphase)
  
  # Collect the rasters in the guessed {RP3} or {slic3r} order.
  CRSS = get_cold_paths(OPHS)
  
  def choose_direction(p, ij):
    # Returns true if the raster identified by the index pair {ij}
    # should be reversed, based on which end is closer to {p}.
    isc,jrs = ij
    oph = OPHS[SCS[isc][jrs]]
    d0 = rn.dist(p, path.pini(oph))
    d1 = rn.dist(p, path.pfin(oph))
    return d0 > d1
    # ....................................................................

  # Join each Alternate directions between adjacent scanlines.
  p_prev = (15,0) # Last position of the nozzle.
  EFS = []
  for CRS in CRSS:
    # Decide orientation {rev} of first raster in {CRS}
    rev = choose_direction(p_prev, CRS[0])

    # Now collect the raster fill paths specified by {CRS}, alternating directions:
    for isc,jrs in CRS:
      oph = OPHS[SCS[isc][jrs]]
      if rev: oph = path.rev(oph)
      EFS.append(oph)
      p_prev = path.pfin(oph)
      rev = not rev
  
  # Assemble the path:
  use_links = True
  ph = path.concat(EFS, use_links, mp_jump)
      
  # SPlit at jumps for coloring:
  OSKS, JMPS = path.split_at_jumps(ph)
  Ylum = color['Ytrace']
  CLRS = hacks.trace_colors(len(OSKS), Ylum)

  plot_contours(c, OCRS, dp, style, color['ghost'])

  # Plot trace components:
  rwdf = style['rwd_fill']
  deco = False
  plot_paths(c, OSKS, dp, CLRS, rwdf,deco, style,color)

  # Plot the jumps:
  clr_jp = color['jphi']
  plot_jumps(c, JMPS, dp, clr_jp, style,color)
  
  # Plot the path endpoints:
  plot_path_endpoints(c, [ph,], dp, style,color)

  # Find the coldest contact and its cooling time:
  ncold = 2 # Number of coldest contacts to show.
  CTScold = contact.coldest(ph, CTS, ncold)
  contact.show_times(sys.stderr, [ph,], CTScold)

  # Write TeX file with fabtimes and cooling times:
  write_tex_params("cold", subfig, [ph,], CTScold)
  
  # Plot coldest contact(s).
  interbc = False
  shadow = True
  trim = False
  plot_contacts(c, CTScold, dp, interbc, shadow,trim, style,color)
  
  # Label path endpoints:
  tdA = (+0.8, -1.4)
  tdB = (-2.0, -0.4)
  pfsize = style['fullfsize']
  plot_path_endpoint_labels(c, ph, "A", tdA, "B", tdB, dp, pfsize, style,color)
 
  return c
  # ----------------------------------------------------------------------

def get_rivers(OPHS, CTS):
  # Retruns a list {GROUPS} that describes the rivers in the filling raster set {OPHS} 
  # defined by the contacts {CTS}. Each element of the list {GROUPS} specifies the
  # rasters that comprise a river. It is a list of pairs {(isc,jrs)}, meaning raster {jrs}
  # of scanline {isc}.
  
  # !!! Should use {raster_regroup.split_by_group} !!!
  GROUPS = ( 
    ( (0,0), (1,0), (2,0), ),
    ( (3,0), (4,0), (5,0), (6,0), (7,0), (8,0), (9,0), (10,0), (11,0), (12,0), (13,0), (14,0), ),
    ( (3,1), (4,1), (5,1), (6,1), (7,1), (8,1), (9,1), (10,1), (11,1), (12,1), (13,1), (14,1), ),
    ( (15,0), (16,0), (17,0), ),
    ( (14,2), (13,2), ),
    ( (13,3), (14,3), (15,1), ),
    ( (9,2), (10,2), (11,2), (12,2), ),
    ( (9,3), (10,3), (11,3), (12,3), ),
    ( (7,2), (8,2), ),
    ( (3,3), (4,3), (5,3), (6,3), ),
    ( (3,2), (4,2), (5,2), (6,2), ),
    ( (0,1), (1,1), (2,1), ),
  )
  return GROUPS
  # ----------------------------------------------------------------------

def get_sub_rivers(OPHS, CTS):
  # Retruns a list {GROUPS} that describes the sub-rivers in the filling raster set {OPHS} 
  # defined by the contacts {CTS}.
  #
  # Each element of the list {GROUPS} specifies the rasters that comprise
  # a sub-river. It is a list of pairs {(isc,jrs)}, meaning raster {jrs} of
  # scanline {isc}.
  
  # !!! Should use {raster_regroup.split_by_group} !!!
  GROUPS = ( 
    ( (0,0), (1,0), (2,0), ),
    ( (0,1), (1,1), (2,1), ),

    ( (3,0), (4,0), (5,0), (6,0), ), 
    ( (3,1), (4,1), (5,1), (6,1), ),
    ( (3,2), (4,2), (5,2), (6,2), ),
    ( (3,3), (4,3), (5,3), (6,3), ),

    ( (7,0), (8,0), ),
    ( (7,1), (8,1), ),
    ( (7,2), (8,2), ),

    ( (9,0), (10,0), (11,0), (12,0), ),
    ( (9,1), (10,1), (11,1), (12,1), ),
    ( (9,2), (10,2), (11,2), (12,2), ),
    ( (9,3), (10,3), (11,3), (12,3), ),

    ( (13,0), (14,0), ),
    ( (13,1), (14,1), ),
    ( (13,2), (14,2), ),
    ( (13,3), (14,3), ),

    ( (15,0), ),
    ( (15,1), ),

    ( (16,0), (17,0), ),
  )
 
  return GROUPS
  # ----------------------------------------------------------------------
 
def get_essential_cut_lines(OPHS, CTS):
  # Returns a list {LNS} of essential cut-lines.
  # Each element of {LNS} is a pair {(icu, t)} where {icu} is the index of the scan-line just above 
  # the cut-line and {t} is a code (2 for essential, 1 for non-essential, 0 for not relevant). 

  LNS = (
    (0,2),
    (3,2),
    (7,2),
    (9,2),
    (13,2),
    (15,2),
    (16,2),
    (18,2),
  )

  return LNS
  # ----------------------------------------------------------------------

def plot_figure_rivers(subfig, style,color):
  # Plots the figure with rasters grouped into rivers.
  #
  #   subfig == "whole" Shows the rivers, color-coded, and the separating contacts.
  #   subfig == "subs"  Shows the sub-rivers defined by the essential scan-lines.
  #
  # Returns the canvas {c} with the figure.
 
  assert subfig == "whole" or subfig == "subs", ("invalid subfig %s" % subfig)
  
  # Get the contours, fillers, links, contacts, and graph of the "turkey" part:
  mp_cont, mp_fill, mp_link, mp_jump = make_move_parms(style);
  OCRS,OPHS,OLKS,CTS,VGS,EGS = paper_example_B.make_turkey(mp_cont, mp_fill,mp_link,mp_jump) 
 
  Bfig, cuxlo, cuxhi = get_turkey_figure_bbox(OCRS, OPHS, OLKS, style)
  dp = (0,0)
  autoscale = True
  c = make_figure_canvas(Bfig, autoscale, style,color)

  # Separate the rasters by scan-line:
  xdir = (1,0)
  ydir = (0,1)
  ystep, yphase = raster.get_spacing_and_phase(OPHS, xdir, ydir)
  SCS = raster.separate_by_scanline(OPHS, xdir, ydir, ystep, yphase)
  
  # Assign group indices to rivers.
  # Each element of the list {GROUPS} below is a river. 
  # Each river is a list of pairs {(isc,jrs)}, meaning raster {jrs} of scanline {isc}.
  
  if subfig == "whole":
    GROUPS = get_rivers(OPHS, CTS)
  elif subfig == "subs":
    GROUPS = get_sub_rivers(OPHS, CTS)
  else:
    assert False, "invalid subfig"
  
  ngr = len(GROUPS)
  igr = 0
  for RIV in GROUPS:
    # Assign group index {igr} to rasters in {RIV}:
    for isc,jrs in RIV:
      oph = OPHS[SCS[isc][jrs]]
      assert path.nelems(oph) == 1
      path.set_group(oph, igr)
    igr += 1
  assert igr == ngr
  
  plot_contours(c, OCRS, dp, style, color['ghost'])

  # Make a color list, one for each river:
  Ylum = color['Ytrace']
  CLRS = hacks.trace_colors(ngr, Ylum)

  # Plot the filling path(s) with color based on group:
  rwdf = style['rwd_fill']
  deco = False
  for oph in OPHS:
    igr = path.get_group(oph)
    assert igr != None, "group was not assigned"
    clr = CLRS[igr]
    plot_paths(c, [oph,], dp, [clr,], rwdf,deco, style,color)

  if subfig == "whole":
    interbc = True # Only inter-block contacts please
    shadow = True
    trim = False
    plot_contacts(c, CTS, dp, interbc, shadow,trim, style,color)
  elif subfig == "subs":
    LNS = get_essential_cut_lines(OPHS, CTS);
    trim = False
    plot_cut_lines(c, LNS, cuxlo,cuxhi,trim, ystep,yphase, dp, style,color)
  else:
    assert False, "invalid subfig"

  return c
  # ----------------------------------------------------------------------
 
def get_generic_blocks(OPHS, SCS, full, hooks, mp_jump):
  # Returns a list {BCS} of snake blocks built from the raster paths in
  # {OPHS}, as could be the input to {HotPath}.
  #
  # The parameter {OPHS} should be a list of all the filling raster elements.
  # each element oriented from left to right.
  #
  # The parameter {SCS} should be a list of lists of indices. The raster
  # element {jrs} (from 0, left to right) on scanline {isc}, as returned by
  # {raster.separate_by_scanline}.
  # 
  # Each block will be constructed from a subset of rasters from {OPHS}
  # that lie on adjacent scan-lines, with exactly one contact between
  # conscutive rasters. The rasters will be connected by the link paths
  # attached to those rasters, if any. Apart from that constraint, the
  # partition of rasters into blocks is arbitrarily chosen by the
  # procedure
  # 
  # If {full} is true, the choices of each block will be all the snake
  # paths that can be built from those rasters (either 2 or 4, snake
  # paths, depending on the number of scan-lines that it spans)
  #
  # If {full} is false, each block will have only one choice,
  # chosen arbitrarily.
  # 
  # If {hooks} is false or the block has a single raster, 
  # the links from the input rasters are copied to the block choices
  # when applicable.  Otherwise the linsk are incorporated 
  # in the block choices. See {make_snake_block_from_rasters}.
  #
  # No new {Contact} objects are created, but the attachments between
  # the existing contacts and paths of {OPHS} that are used in blocks,
  # as given by {path.get_contacst} and {contact.get_side_paths}, are
  # broken and the contacts are attached to the choices of the blocks
  # instead.
  #
  # The procedure also sets group index (as returned by
  # {path.get_group}) of all original raster paths in {OPHS} to the
  # sequential index of the block in which they were used; or to {None}
  # if they were not used in any block.
  
  # !!! Should use {raster_regroup.split_by_group} !!!
  
  # {SNIXSS} is a list of lists of pairs.  If {SNIXSS[ibc][irs]} is {(isc,jrs)}, 
  # it means that raster {irs} (from bottom) of block {ibc} 
  # is raster {jrs} of scanline {isc}.
  SNIXSS = ( 
    ( (0,0), (1,0), (2,0), ),
    ( (0,1), (1,1), (2,1), ),

    ( (3,0), (4,0), (5,0), ), 
    ( (6,0), (7,0), (8,0), ),
    ( (9,0), (10,0), (11,0), ),
    ( (12,0), (13,0), (14,0), ),
    
    ( (3,1), (4,1), (5,1), (6,1), (7,1), (8,1), ),
    ( (9,1), (10,1), (11,1), (12,1), (13,1), (14,1), ),

    ( (3,2), (4,2), (5,2), (6,2), ),

    ( (7,2), (8,2), ),

    ( (9,2), (10,2), ),
    ( (11,2), (12,2), ),
    ( (13,2), (14,2), ),

    ( (3,3), (4,3), (5,3), (6,3), ),
    ( (9,3), (10,3), (11,3), (12,3), ),

    ( (13,3), (14,3), (15,1), ),

    ( (15,0), (16,0), (17,0), ), 
  )
  
  # Clear the goup index of all raster paths in {OPHS}:
  for oph in OPHS: path.set_group(oph, None)
  
  nsc = len(SCS) # Number of scan-lines.
  nbc = 0 # Number of blocks created.
  
  BCS = [] # List of generic blocks.
  RSrest = []
  for SNIXSj in SNIXSS:
    # SNIXSj is a list of index pairs, specifying the rasters of one block.
    # Create the block {bcj}
    icumax = nsc # Cut-line above last scan-line.
    OPBSj = collect_and_tag_snake_rasters(SNIXSj, SCS, OPHS, icumax, nbc, RSrest)
    bcj = make_snake_block_from_rasters(OPBSj, full, hooks, nbc, mp_jump)
    assert bcj != None
    BCS.append(bcj) 
    nbc += 1

  # All rasters must have been used:
  assert len(RSrest) == 0
  
  sys.stderr.write("got %d blocks\n" % len(BCS))
  for jbc in range(len(BCS)):
    bcj = BCS[jbc]
    if full:
      assert isinstance(bcj, block.Block)
      nchj = block.nchoices(bcj)
      nrsj = path.nelems(block.choice(bcj,0))
      sys.stderr.write("  block %d has %d choices of %d rasters\n" % (jbc,nchj,nrsj))
    else:
      assert isinstance(bcj, path.Path)
      nrsj = path.nelems(bcj)
      sys.stderr.write("  path %d has %d rasters\n" % (jbc,nrsj));

  return BCS
  # ----------------------------------------------------------------------

def collect_and_tag_snake_rasters(SNIXS, SCS, OPHS, icumax, ibc, RSrest):
  # The parameter {OPHS} should be a list of all the filling raster elements.
  # each element oriented from left to right.
  #
  # The parameter {SCS} should be a list of lists of indices. The raster
  # element {jrs} (from 0, left to right) on scanline {isc}, as returned by
  # {raster.separate_by_scanline}.
  # 
  # Parameter {SNIXS} is a list of pairs {(isc,jrs)} that identify the rasters of 
  # one canonical snake block.  Returns the list {OPKS} of those raster paths,
  # in the same orientation as they are in {OPHS}.
  #
  # However, only rasters below the cut-line {icumax} are used in the snake block.
  # Rasters above that cut-line are instead appended to the list {RSrest} 
  # with no change.  If all rasters specified by {SNIXS} are above that cut-line,
  # returns the empty list.
  #
  # Also sets the group index of all the selected raster paths to {ibc}.

  OPKS = [] # Raster elements in this snake.
  for isc, jrs in SNIXS:
    ors = OPHS[SCS[isc][jrs]]
    if isc >= icumax:
      RSrest.append(ors)
    else:
      path.set_group(ors, ibc)
      OPKS.append(ors)
  return OPKS
  # ----------------------------------------------------------------------

def make_snake_block_from_rasters(OPBS, full, hooks, ibc, mp_jump):
  # Builds and returns a snake block from the raster paths in the list
  # {OPBS}. These rasters should be in consecutive scan-lines and all
  # oriented from left to right.
  #
  # The raster paths in {OPBS} should all have group index {ibc}, which
  # should be distinct from the groups of al other original raster paths
  # and of all choices of previously created blocks.
  #
  # If {full} is true, returns a block whose choices are all snake paths
  # that can be built from the rasters on {OPBS}; usually 2 if {OPBS}
  # has a single raster, and 4 if it has two or more. If {full} is
  # false, returns a block with only one choice, specifically the snake
  # path that begins with the bottom-most raster oriented from left to
  # right.
  #
  # However, if the list {OPBS} is empty, returns {None}.
  #
  # If {hooks} is false, or {OPBS} has only one raster, copies the
  # applicable link paths from the raster paths in {OPBS} to the choices
  # of the new block. If {hooks} is true and {OPBS} has two or more
  # rasters, includes the inward-pointing links at the ends of the path
  # in the appropriate choices, to improve adhesion to the contour
  #
  # Also modifies the contact-path attachments of to refer to the block
  # choices instead of those raster paths. Specifically, any contact
  # {ct} that is attached to a raster path {ors} of {OPBS} and to some
  # path not in {OPBS} (original raster or choice of a previosuly
  # created block) will be detached from {ors} and attached to all the
  # choices of the new block. Assumes that the only existing contacts
  # are between rasters, and any contacts with contour elements, if any,
  # will be added later.

  # Collect raster paths {OPHS0,OPHS1} of the snake and (if {full}) of the "mirror" snake:
  OPHS0 = [] # Raster elements in this snake.
  OPHS1 = [] # Reversed raster elemens in this snake (if {full}).
  rev = False # Should reverse the next raster?
  for opb in OPBS:
    assert path.get_group(opb) == ibc
    if rev: opb = path.rev(opb)
    OPHS0.append(opb)
    if full: OPHS1.append(path.rev(opb))
    rev = not rev

  # Did we get any rasters at all?
  if len(OPHS0) == 0: return None
  
  if hooks and len(OPBS) >= 2:
    # Add the inward-pointing links, if any:
    OPHS0 = add_hooks_to_snake_rasters(OPHS0)
    if len(OPHS1) != 0:
      OPHS1 = add_hooks_to_snake_rasters(OPHS1)
  
  # Now join the rasters in each of {OPHS0} and {OPHS1} into a snake block
  use_links = True
  ph0 = path.concat(OPHS0, use_links, mp_jump) # The snake path.
  if not full: 
    bc = block.from_paths([ph0,])
  else:
    ph1 = path.concat(OPHS1, use_links, mp_jump) # The other snake path.
    if len(OPHS0) == 1:
      bc = block.from_paths([ph0, ph1])
    else:
      bc = block.from_paths([ph0, path.rev(ph0), ph1, path.rev(ph1)])

  reattach_snake_block_contacts(OPHS0, bc, ibc)
 
  for ich in range(block.nchoices(bc)): 
    och = block.choice(bc, ich)
    path.set_group(och, ibc)
    path.set_name(och, "B%d.%d" % (ibc,ich), moves = False)

  block.set_name(bc, "B%d" % ibc)

  return bc
  # ----------------------------------------------------------------------
 
def add_hooks_to_snake_rasters(OPHS):
  # Given list {OPHS} of two or more stacked rasters that are to be made
  # into a block, oriented with alternating directions, returns a copy
  # of the list with any "hook" links prepended and/or appended.
  # 
  # A hook is a link from the end of the last raster to the start of the
  # next-to-last one, or from the end of the second raster to the
  # start of the first one.  The hools are slightly shortened at the 
  # internal end.
  
  nrs = len(OPHS) # Number of rasters in the stack
  assert nrs >= 2
  
  olk0 = path.get_connecting_link(OPHS[1], OPHS[0])
  if olk0 != None: 
    olk0 = make_hook_from_link(olk0)
  
  olk1 = None if nrs == 2 else path.get_connecting_link(OPHS[nrs-1], OPHS[nrs-2])
  if olk1 != None: 
    olk1 = path.rev(make_hook_from_link(path.rev(olk1)))
  
  OPHSnew = OPHS.copy()
  if olk0 != None: OPHSnew = [olk0,] + OPHSnew
  if olk1 != None: OPHSnew = OPHSnew + [olk1,]
  
  return OPHSnew
  # ----------------------------------------------------------------------
    
def make_hook_from_link(olk):
  # Converts a link path {olk} into a hook,
  # by shaving of a little from the beginning and simplifying the result:

  # Trim end:
  wd = move.width(path.elem(olk, 0))
  dtrim = 0.50*wd
  olk = path.trim(olk, dtrim, 0)
  assert olk != None
  stol = 0.030*wd
  olk = path.simplify(olk, stol, None)
  return olk
  # ----------------------------------------------------------------------

def reattach_snake_block_contacts(OPHS, bc, ibc):
  # Given a list of raster-fill paths {OPHS} that were used to make one
  # snake block. detaches them from their attached contacts and attaches
  # the latter to the choices of the block {bc}.
  #
  # On input, all the raster paths in {OPHS} should have group index
  # {ibc}, and all other paths (original raster paths or choices of
  # previoulsy created blocks) must have group index {None} or strictly
  # less than {ibc}.  The group index of the choices of {bc} is not used. 
  #
  # Assumes that the only existing contacts are between rasters, and any
  # contacts with contour elements, if any, will be added later.

  debug = False

  # Get all the contacts {CTS} whose sides were used in this block:
  CTS = set()
  for ophi in OPHS:
    for isd in range(2):
      CTS = set.union(CTS, path.get_contacts(ophi, isd))
  
  # Fix the attachments of all those contacts:
  for ct in CTS:
    if debug: contact.show(sys.stderr, "  ct = ", ct, "\n", 0)
    for ksd in range(2):
      # Get the paths that are attached to side {ksd} of {ct}.
      # These are either a single old raster path, possibly used in this block, 
      # or one or more choices of a previously conctructed block.
      SDPSk = contact.get_side_paths(ct, ksd)
      assert len(SDPSk) > 0 # Since the contact must originally have had have both sides on raster paths of {OPHSD}.
      
      # Set {ibck} to the index of the block that contains side {k} of {ct} (which may be {ibc})
      # or to {None} if that side is still not in any block:
      # sys.stderr.write("ibc = %d contact %s side %d SPDSk = %s:\n" % (ibc, contact.get_name(ct),isd, str(SDPSk)))
      igrk = None
      for phkj, drkj, imvkj in SDPSk:
        igrkj = path.get_group(phkj)
        if debug: path.show(sys.stderr, ("    side %d @ " % isd), phkj, (" igr = %s\n" % str(igrkj)), True, 0,0,0)
        # sys.stderr.write("  igrkj = %s igrk = %s\n" % (str(igrkj),str(igrk)))
        if igrkj == ibc or igrkj == None:
          # This side of the contact is in one of original rasters, either used by this
          # block or not yet used by any block:
          assert igrk == None # Since there should be only one such raster.
          assert len(SDPSk) == 1 # Since there should be only one such raster.
          assert imvkj == 0 # Since those original rasters paths have length 1.
        else:
          # This side of the contact is in a choice of a previously constructed block.
          assert igrkj < ibc # Assuming blocks are numbered consecutively by creation.
          assert igrk == None or igrk == igrkj # Since all side paths must be in the same block.
        igrk = igrkj 
        
      if igrk == ibc:
        # Path attachments on side {ksd} of {ct} must be changed from old raster to choices of new block:
        contact.clear_side_paths(ct, ksd)
        for phkj, drkj, imvkj in SDPSk: path.clear_contacts(phkj, ksd)
          
        # Attach side {ksd} of {ct} to the choices of {bc}
        mvk = contact.side_move(ct, ksd)
        for ich in range(block.nchoices(bc)):
          ochi = block.choice(bc, ich)
          jmvki = path.find_move(ochi, mvk)
          if jmvki != None:
            path.add_contact(ophi, ksd, ct)
            contact.add_side_path(ct, ksd, ochi, jmvki)
  return
  # ----------------------------------------------------------------------

def plot_figure_hotfill(subfig, style,color):
  # Plots the figure with canonical paths and canonical bandpath.
  #
  #   subfig "sta{NNNNN}" = state during {hotfill.solve}.
  #   subfig "sol{NNNNN}" = actual solution found by HotFill with {Delta = NNNN/100}.
  #   subfig "tab{NNNNN}" = tableau during {hotfill.solve} with {Delta = NNNN/100}.
  #
  # where {NNNNN} is {1000*Delta} formatted as 5 digits, zero-padded.
  # 
  # Returns the canvas {c} with the figure.
  
  subpre = subfig[0:3]
  assert subpre == "sol" or subpre == "sta" or subpre == "tab", ("invalid subfig %s" % subfig)

  # Get the contours, fillers, links, contacts, and graph of the "turkey" part:
  mp_cont, mp_fill, mp_link, mp_jump = make_move_parms(style);
  OCRS,OPHS,OLKS,CTS,VGS,EGS = paper_example_B.make_turkey(mp_cont, mp_fill,mp_link,mp_jump) 
  
  # Sort and separate by scanline:
  xdir = (1,0)
  ydir = (0,1)
  ystep, yphase = raster.get_spacing_and_phase(OPHS, xdir, ydir)
  OPHS = raster.sort_by_scanline(OPHS, xdir, ydir, ystep, yphase)
  SCS = raster.separate_by_scanline(OPHS, xdir, ydir, ystep, yphase)
  nsc = len(SCS)

  # Assign the cooling time limits:
  assert len(subfig) == 8
  Delta = float(subfig[3:])/1000
  for ct in CTS: contact.set_tcool_limit(ct, Delta)

  # Compute the {hotfill} solution:
  
  # Describe and check the input data:
  maxband = 1000 # No band size limit.
  quick = False # Use quadratic greedy.
  hotfill.describe_and_check_input(sys.stderr, OPHS, mp_jump, maxband, quick)
  fph_best, z_best, i_best, j_best, BPHS_best, BTCV, TUK = hotfill.solve(OPHS, mp_jump, maxband, quick)
  assert j_best == nsc
  
  dp = None
  
  if subpre == "sta" or subpre == "tab":
    z = 0
    i = 8
    j = 11
    T_zij, u_zij, k_zij = TUK[z][i][j]
    fph, BPHS = hotfill.recover_fullpath(BTCV, TUK, u_zij, k_zij, i, mp_jump)
    bph, CTS_lo, TCVS_lo, CTS_hi, TCVS_hi = BTCV[z][i][j]
  elif subpre == "sol":
    z = z_best
    i = i_best
    j = j_best
    fph = fph_best; BPHS = BPHS_best
    bph = None; CTS_lo = None; TCVS_lo = None; CTS_hi = None; BCTS_hi = None;
    
  # Get the indices of the cut-lines that separate the bandpaths of {fph,pbh}:
  LNS = get_hotfill_cut_lines(TUK, z, i, j);

  if subpre == "tab":

    elsz = 1.0         # X and Y size of a tableau cell.
    tbsz = (nsc+1)*elsz  # Total X and Y size of tableau, minus row and column labels. 
    B = ((0, 0), (tbsz, tbsz))  # Bounding box of tabeau.
    # Compute the figure's bounding box {Bfig}, allowing for row and column labels:
    Bfig = rn.box_expand(B, (2.6,2.2), (2.6,2.2))
    autoscale = True
    c = make_figure_canvas(Bfig, autoscale, style,color)
  
    plot_tableau(c, nsc, TUK, BTCV, z, i, j, style,color)
    
  elif subpre == "sta" or subpre == "sol":

    Bfig, cuxlo, cuxhi = get_turkey_figure_bbox(OCRS, OPHS, OLKS, style)
    autoscale = True
    c = make_figure_canvas(Bfig, autoscale, style,color)
  
    plot_contours(c, OCRS, dp, style, color['ghost'])
  
    # Plot the cut-lines that separate the bandpaths of {fph,pbh}:
    trim = True
    plot_cut_lines(c, LNS, cuxlo,cuxhi,trim, ystep,yphase, dp, style,color)

   # Get the link paths {PLKS}, jump moves {JMPS}, and the bandpaths {PPHS} to plot
    clr_ph0 = color['fil0'] # Color for even bandpaths in the fullpath {fph_uki}.
    clr_ph1 = color['fil1'] # Color for odd bandpaths in the fullpath {fph_uki}.
    CLRS_ph = []
    PPHS = []
    PLKS = []
    JMPS = []
    for iph in range(len(BPHS)):
      if iph % 2 == 0:
        # Should be a bandpath of {fph_uki}:
        oph = BPHS[iph]
        assert oph != None and not isinstance(oph, move.Move)
        PPHS.append(oph)
        CLRS_ph.append(clr_ph0 if len(PPHS) % 2 == 0 else clr_ph1)
      else:
        mvph = BPHS[iph]
        if isinstance(mvph, move.Move):
          assert move.is_jump(mvph)
          JMPS.append(mvph)
        else:
          # Assume it is link path:
          PLKS.append(mvph)

    # Plot the link paths:
    bicolor = False
    plot_links(c, PLKS, dp, ystep,yphase, bicolor, style,color)

    # Plot the paths:
    rwdf = style['rwd_fill']
    deco = False
    pfsize = style['fullfsize']

    if subpre == "sta":
      # Add the candidate bandpath {bph}: 
      clr_bph = color['hifi'] # Color for the {(z,i,j)} bandpath {bph}.
      CLRS_ph.append(clr_bph)
      PPHS.append(bph)

    plot_paths(c, PPHS, dp, CLRS_ph, rwdf,deco, style,color)

    # Plot the inter-band jumps:
    clr_jp = color['jphi']
    plot_jumps(c, JMPS, dp, clr_jp, style,color)

    if subpre == "sta":

      # Plot the unused rasters:
      UPHS = [ OPHS[iph] for isc in range(j,nsc) for iph in SCS[isc]]
      clr_ur = color['ghost']
      deco = False
      plot_paths(c, UPHS, dp, [clr_ur,], rwdf,deco, style,color)

      # Plot the contacts between {fph} and {bph}:
      interbc = False
      shadow = True
      trim = False
      plot_contacts(c, CTS_lo, dp, interbc, shadow,trim, style,color)

      # Plot and label the endpoints of {fph} and {bph}:
      plot_path_endpoints(c, [fph,bph,], dp, style,color)
      tdA = (-1.8, -1.6)
      tdB = (-1.8, -1.1)
      plot_path_endpoint_labels(c, fph, "A", tdA, "B", tdB, dp, pfsize, style,color)
      tdC = (-1.8, -0.1)
      tdD = (+0.7, -0.4)
      plot_path_endpoint_labels(c, bph, "C", tdC, "D", tdD, dp, pfsize, style,color)

    elif subpre == "sol":

      # Plot and label the endpoints of {fph}:
      plot_path_endpoints(c, [fph,], dp, style,color)
      tdA = (-1.8, -1.6)
      tdB = (-1.8, 00.0)
      plot_path_endpoint_labels(c, fph, "A", tdA, "B", tdB, dp, pfsize, style,color)
    
      # Find the coldest contact(s) and their cooling times:
      ncold = 2 # Number of coldest contacts to show.
      CTScold = contact.coldest(fph, CTS, ncold)
      contact.show_times(sys.stderr, [fph,], CTScold)

      # Write TeX file with fabtimes and cooling times:
      write_tex_params("hotfill", subfig, [fph,], CTScold)

      # Plot coldest contact(s).
      interbc = False
      shadow = True
      trim = False
      plot_contacts(c, CTScold, dp, interbc, shadow,trim, style,color)
    else:
      assert False
  else:
    assert False

  return c
  # ----------------------------------------------------------------------

def get_hotfill_cut_lines(TUK, z, i, j):
  # Given the tableay {TUK} of {hotfill},
  # Returns a list {LNS} of the cut-lines that 
  # separate the bandpaths of the fullpath {(z,i,j)}.
  
  if j == 0:
    LNS = [ ]
  else:
    T_zij, u_zij, k_zij = TUK[z][i][j]
    LNS = get_hotfill_cut_lines(TUK, u_zij, k_zij, i)
  LNS.append( (j,1), )
  return LNS
  # ----------------------------------------------------------------------

def plot_tableau(c, nsc, TUK, BTCV, z, i, j, style,color):
  # Given the tables {TUK} and {BCTF} used by {hotfill},
  # plots on canvas {c}  the tableau as it was when the procedure was considering
  # the concatenation of {F[u][k][i]} with {B[z][k][j]}.
  
  zX = z  # The {z} argument to {min_fullpath}
  iX = i  # The {i} argument to {min_fullpath}
  jX = j  # The {j} argument to {min_fullpath}.
  
  elsz = 1.0         # X and Y size of a tableau cell.
  tbsz = (nsc+1)*elsz  # Total X and Y size of tableau, minus row and column labels. 

  dp = (0, 0)
  
  elmrg = 0.11 # Indent of element boxes.
  wdess = 0.15 # Width of essential cut-line strokes.
  wddot = 0.50 # Diameter of dots.
  labwd = 5.00 # Estra width for index labels.
  esswd = 0.60 # Overshoot of essential row/cols backgroud stripe.
  loXY = 0 - esswd            # Low coordinate of highlighted row/cols backgroud stripe.
  hiXY = (nsc+1)*elsz + esswd # High coordinate of highlighted row/cols backgroud stripe.
  
  wd = 0.05 # Width of lines.
  clr_unused = pyx.color.rgb( 0.830, 0.820, 0.800 ) # Entries unused because {i >= j}.
  clr_nothot = pyx.color.rgb( 1.000, 0.000, 0.000 ) # Entries that are {None} because of Tcool violation.
  clr_exists = pyx.color.rgb( 0.000, 0.700, 0.000 ) # Entries defined.
  clr_tocomp = pyx.color.rgb( 0.000, 0.000, 0.000 ) # Entries not computed yet.
  clr_phnext = pyx.color.rgb( 0.000, 1.000, 1.000 ) # Dot on entry that is to be computed next.
  clr_phneed = pyx.color.rgb( 1.000, 1.000, 1.000 ) # Dot on entry that is needed for that entry.
  clr_cutcmp = pyx.color.rgb( 0.900, 0.900, 0.900 ) # Highlight of rows/columns that {iX} or {jX}.

  def plot_tableau_row_col_deco(i):
    # Plots the index of row and/or column{i}, if relevant, and draws line across it.

    if i == 0 or i == iX or i == jX or i == nsc:
      # Row {i} is interesting. Column {i} is intersting unless {i==jX}:
      
      pri = (loXY, (nsc - i + 0.5)*elsz)  # Left endpoint of line though row {i}.
      qri = (hiXY, pri[1])                # Right endpoint of line though row {i}.
      pci = ((i + 0.5)*elsz, loXY)        # Bottom endpount of line through column {i}.
      qci = (pci[0], hiXY)                # Top endpoint of line through column {i}.

      # Column {i} is interesting, draw line and column labels:
      hacks.plot_line(c, pci, qci, dp, clr_cutcmp, wdess, None)
      for txY in loXY - 1.3, hiXY + 0.5:
        dxci = 0.4 if i <= 9 else 0.8
        tpci = (pci[0] - dxci, txY)
        hacks.plot_text(c, (r"\texttt{%d}" % i), tpci, dp, style['tbixfsize'], None)

      if i != jX:
        # Row {i} is intersting, draw line and row labels:
        hacks.plot_line(c, pri, qri, dp, clr_cutcmp, wdess, None)
        dxri = 0.1 if i <= 9 else 0.8     # Extra displacement of label on left
        for txX in loXY - 1.0 - dxri, hiXY + 0.3:
          tpri = (txX, pri[1] - 0.4)
          hacks.plot_text(c, (r"\texttt{%d}" % i), tpri, dp, style['tbixfsize'], None)

    return
    # ....................................................................

  def plot_tableau_cell(i, j):
    # Plots the cell {i,j} of the tableau.

    # Plot the box:
    pij = (j*elsz, (nsc-i)*elsz)
    qij = rn.add(pij, (elsz, elsz))
    Bij = (pij, qij)
    Bij = rn.box_expand(Bij, (-elmrg, -elmrg), (-elmrg, -elmrg))
    TUK_zij = TUK[z][i][j]
    if i == nsc or j == 0 or i >= j:
      clr = clr_unused
    elif j > jX or (j == jX and i >= iX):
      clr = clr_tocomp
    elif TUK_zij != None and TUK_zij[0] != None and TUK_zij[0] != +inf:
      clr = clr_exists
    else:
      clr = clr_nothot
    hacks.plot_box(c, Bij, dp, clr)

    # Plot the dot, if any:
    cij = rn.add(pij, (0.5*elsz, 0.5*elsz)) # Center of cell.
    if i == iX and j == jX:
      # Paint a yellow dot to indicate the entry to be computed bext:
      clr = clr_phnext
    elif i < j and j == iX:
      clr = clr_phneed
    else:
      clr = None

    if clr != None:
      hacks.plot_line(c, cij, cij, dp, clr, wddot, None)
    return
    # ....................................................................

  for ipass in range(2):
    # When {ipass} is 0, draws lines for essential rows/colums,
    # When {ipass} is 1, paints the tableau entries.
  
    rnext = 0 # The next essential cutline greater than or equal to {i} is {L[r]}
    for i in range(nsc+1):

      if ipass == 0:
        # Paint the background of rows that are essential cut-lines:   
        plot_tableau_row_col_deco(i)
      else:
        # Paint row {j} of the tableau:
        for j in range(nsc+1):
          plot_tableau_cell(i, j)
  
  return c
  # ----------------------------------------------------------------------

def make_block_from_short_contour(ocr, EPS, ibc, gap, wdf):
  # Given a contour path {ocr}, makes a block with two choices, which
  # are the path starting from the top and the paths starting from the
  # bottom. Each path has its ends trimmed so as to leave the specified
  # {gap}.  The group index of the paths is set to {ibc}.
  #
  # The parameter {EPS} should be a list of pairs {(p,igr)} as extracted
  # by {get_all_move_endpoints} from the filling blocks {BCS_fill}.
  # Assumes that the traces in the filling blocks (including links
  # included in their choices) have moves of width {wdf}.

  nmv = path.nelems(ocr)
  assert nmv >= 3
  assert path.pini(ocr) == path.pfin(ocr)
  mp_cont = move.parameters(path.elem(ocr,0))
  wd_cont = move_parms.width(mp_cont)
  
  P = [None, None] # The two paths.  
  for ich in 0, 1:

    # Rebase contour at a point near the top or bottom, but not quite:
    alf = 0.65*pi + pi*ich
    och = rebase_contour(ocr, alf)

    # Trim ends to make the contour endpoints visible and avoid the button:
    och = path.trim(och, gap/2, gap/2)
    assert och != None

    # Simplify the path, preserving its endpoints:
    stol = 0.030*wdf
    och = path.simplify(och, stol, None)

    path.set_group(och, ibc)
    path.set_name(och, "B%d.%d" % (ibc,ich), moves = True)
    P[ich] = och
  
  bc = block.from_paths(P)
  block.set_name(bc, "B%d" % ibc)
  
  return bc
  # ----------------------------------------------------------------------

def make_blocks_from_tall_contour(ocr, EPS, BCS_fill, gap, wdf):
  # Given a contour path {ocr}, breaks it into sections and creates a block  
  # with two choices from each section {S}, namely {S} and {rev(S)}.  Each path has its ends trimmed
  # so as to leave the specified {gap} between sections. 
  #
  # The parameter {EPS} should be a list of pairs {(p,igr)} as extracted
  # by {get_all_move_endpoints} from the filling blocks {BCS_fill}.
  # Assumes that the traces in the filling blocks (including links
  # included in their choices) have moves of width {wdf}.

  nmv = path.nelems(ocr)
  assert nmv >= 3
  assert path.pini(ocr) == path.pfin(ocr)
  mp_cont = move.parameters(path.elem(ocr,0))
  wd_cont = move_parms.width(mp_cont)
  
  # Rebase contour at a point near the top or bottom, but not quite:
  alf = 0.52*pi
  och = rebase_contour(ocr, alf)

  # Break up the contour path into sections syncronized with adjacent blocks:
  OSCS = chop_contour_path(och, EPS, wdf)
  assert len(OSCS) > 0

  ### # Concatenate sections in pairs to make the blocks:
  ### use_links = True
  ### och = path.concat(OSCS, use_links, None) # Sections should have joining links.

  # Make blocks from the sections:
  BCS_cont = []
  for osc in OSCS:
    bc = block.from_paths((osc, path.rev(osc),))
    BCS_cont.append(bc)
  
  return BCS_cont
  # ----------------------------------------------------------------------
    
def rebase_contour(ocr, alf):
  # Returns a copy of the contour path {ocr}, assumed closed,
  # but starting at the extremal vertex in the direction that makes an
  # angle {alf} with the X-axis.

  nmv = path.nelems(ocr)
  
  # Collect the points of the path, without the closing one:
  q = [ move.pini(path.elem(ocr,imv)) for imv in range(nmv) ]
  ydir = (cos(alf), sin(alf))
  
  # Locate the index {i_max} of the extremal vertex in direction {ydir}:
  y_max = -inf
  i_max = None
  for imv in range(nmv):
    yi = rn.dot(q[imv], ydir)
    if yi > y_max: y_max = yi; i_max = imv
  assert i_max != None and i_max >= 0 and i_max < nmv

  # Copy the moves from {i_max}:
  OMVS = [ path.elem(ocr, (i_max + k) % nmv) for k in range(nmv) ]

  # Assemble the new path:
  och = path.from_moves(OMVS)
  return och
  # ----------------------------------------------------------------------

def chop_contour_path(oph, EPS, wdf):
  # Rebuilds the path {oph}, assumed to be a piece of contour, as a list
  # of {OSCS} of /section/ paths, each adjacent to a single filling block. The group of each section path is copied
  # from that filling block. Each section path may have link
  # paths connecting it to adjacent sections
  #
  # Assumes that {EPS} is a list of pairs {(p,igr)} as extracted by
  # {get_all_move_endpoints} from the filling blocks. Assumes that the
  # traces in the filling blocks (including links included in their
  # choices) have moves of width {wdf}.
  #
  # The procedure uses {atomize_contour_path} to split {oph} into a list
  # of /atomic/ paths/. Then each maximal sequence of atomic paths
  # assigned to the same integer group {igr} is concatenated to form a
  # single section path {osc}, with the same group code {igr}. The
  # result is the list of those sections.
  #
  # The unassigned atomic paths are concatenated too to form link paths
  # for the adjacent section paths. However any unassigned atomic paths
  # before the first assigned one are also assigned to the same group,
  # and ditto at the other end.
  #
  # The procedure fails if it cannot find at least one assigned atom.
 
  OAPS = atomize_contour_path(oph, EPS, wdf)
  nap = len(OAPS)
  assert nap >= 1
  
  # Fix groups of initial and final unassigned paths:
  iap = 0
  while iap < nap and path.get_group(OAPS[iap]) == None: iap += 1
  assert iap >= 0 and iap < nap
  igr0 = path.get_group(OAPS[iap])
  for jap in range(iap): path.set_group(OAPS[jap], igr0)
  
  iap = nap-1
  while iap >= 0 and path.get_group(OAPS[iap]) == None: iap -= 1
  assert iap >= 0 and iap < nap
  igr1 = path.get_group(OAPS[iap])
  for jap in range(iap+1, nap): path.set_group(OAPS[jap], igr1)
  
  # Combine maximal sequences of assigned and unassigned atoms:
  OSCS = []
  OAPSi = [OAPS[0]] # Atomic paths in current section.
  olk = None # Link path before current section, or {None}.
  for iap in range(1, nap+1):
    igr_prev = path.get_group(OAPS[iap-1]) # Group of previous atomic path.
    oapi = None if iap == nap else OAPS[iap]
    igr_this = None if oapi == None else path.get_group(oapi)
    if igr_this != igr_prev:
      use_links = False
      sys.stderr.write("  completed section with %d atomic paths of group %s\n" % (len(OAPSi),igr_prev))
      osc = path.concat(OAPSi, use_links, None)
      path.set_group(osc, igr_prev)
      OAPSi = [ ]
      if igr_prev == None:
        # {osc} is a link path. Save it for later.
        olk = osc
        if len(OSCS) >= 1:
          path.add_link(path.rev(OSCS[-1]), path.rev(olk))
      else:
        # {osc} is a section path.
        if olk != None:
           path.add_link(osc, olk)
        olk = None
        OSCS.append(osc)
    # Either way, append this atomic path to the list of the current section:
    if oapi != None: OAPSi.append(oapi)
  assert len(OAPSi) == 0 

  sys.stderr.write("got %d sections of contour\n" % len(OSCS))
  GRS = [ path.get_group(osc) for osc in OSCS ]
  show_group_stats(GRS)

  return OSCS
  # ----------------------------------------------------------------------
  
def atomize_contour_path(oph, EPS, wdf):
  # Breaks the path {oph}, assumed to be part of the slice's contour,
  # into a list of /atomic/ paths, most of them syncronized with moves
  # of adjacent filling blocks.
  #
  # Assumes that {EPS} is a list of pairs {(p,igr)} as extracted by
  # {get_all_move_endpoints} from the filling blocks. Assumes that the
  # traces in the filling blocks (including links included in their
  # choices) have moves of width {wdf}.
  # 
  # If the two endpoints of an atomic path are projection of endpoints from {EPS} with the same group
  # code {igr} , the atomic path is temporarily assigned to that same group;
  # otherwise it is /unassigned/ (assigned to group {None}).

  OAPS = [ ]
  nmv = path.nelems(oph)
  if nmv == 0: return OAPS
  
  wdc = move.width(path.elem(oph,0)) # Width of contour traces, assumed uniform.
  dcmax = 1.5*(wdc + wdf)/2          # Max distance from midline to filling endpoint.
  
  # Compute a box {B} that will contains all "close enough" points of {EPS}:
  B = path.bbox([oph,])
  B = rn.box_expand(B, (dcmax,dcmax), (dcmax,dcmax))
  
  # Project onto the midline of {oph} the points of {EPS} that are
  # within distance {dcmax} of the path, obtaining the list {CPS}.
  # Each element of {CPS} is a quadruple {(q, imv, t, igr)} where {q} is a point of 
  # {EPS} projected onto the midline of {oph}, {imv} is the index of a move of
  # {oph} that contains {q}, {t} is total length of the path from {pini(oph)} to {p},
  # and {igr} is the group index of the block that has the endpoint that was projected.
  
  # Collect the nearby filler endpoints and project them onto contour:
  CPS = set()
  for p, igr in EPS:
    if rn.box_inside(p, B):
      q, imv, t, d = path.find_nearest_midline_point(oph, p)
      if d < dcmax:
        CPS.add((q, imv, t, igr))
  CPS = list(CPS)
  list.sort(CPS, key = lambda x: x[2])
  ncp = len(CPS)
  sys.stderr.write("found %d nearby filler endpoints projected onto contour\n" % ncp)
  
  # Discard coincident projections, retaining those with non-{None} group index:
  CPS_new = []
  q_prev = None
  imv_prev = None
  t_prev = None
  igr_prev = None
  nrp = 0 # Count of repeated points.
  for q, imv, t, igr in CPS:
    if q == q_prev:
      assert abs(t - t_prev) < 1.0e-8
      # Ignore this point but retain the non-{None} group:
      if igr != None: igr_prev = igr
      nrp += 1
    else:
      # Store {q_prev} and update:
      if q_prev != None:
        CPS_new.append((q_prev, imv_prev, t_prev, igr_prev))
      q_prev = q; imv_prev = imv; t_prev = t; igr_prev = igr
  # Flush the last one:
  if q_prev != None:
    CPS_new.append((q_prev, imv_prev, t_prev, igr_prev))
  sys.stderr.write("rejected %d repeated projected points\n" % nrp)
  
  CPS = CPS_new
  ncp = len(CPS)
  sys.stderr.write("retained %d projected points\n" % ncp)
  
  GRS = [ igr for q, imv, t, igr in CPS ]
  show_group_stats(GRS)
  
  # Now break {oph} at the points of {CPS} to obtain the atomic paths:
  
  q_prev = path.pini(oph) # End point of previous atomic path (even if trivial).
  imv_prev = 0 # Move of {oph} that includes {q_prev}.
  igr_prev = None # Group of previous projected endpoint.
  nap_0 = 0 # Number of atomic paths with zero length.
  nap_1 = 0 # Number of atomic paths with exactly one non-zero move.
  for icp in range(1,ncp+1):
    if icp < ncp:
      # Get netx projected endpoint:
      qi, imvi, ti, igri = CPS[icp]
    else:
      # Sentinel entry to flush last section.
      qi = path.pfin(oph); imvi = nmv-1; ti = +inf; igri = None
    oapi = path.extract_section(oph, imv_prev,q_prev, imvi,qi)
    
    # Assign to a group {igri} if both endpoints are projected from 
    # that group, otherwise assign to {None}
    sys.stderr.write("  igr_prev = %s igri = %s\n" % (igr_prev,igri))
    if igri != None and igri == igr_prev:
      # Path may be a projection of a filling move.
      path.set_group(oapi, igri)
    else:
      # Path is not a projection of a filling move; don't simplify it.
      path.set_group(oapi, None)
    
    assert path.pini(oapi) != path.pfin(oapi) # Because we eliminated repeated points.
    # Add atomic path to {OAPS}
    OAPS.append(oapi)
    if path.nelems(oapi) == 1: nap_1 += 1
    
    igr_prev = igri
    imv_prev = imvi
    q_prev = qi

  sys.stderr.write("built %d atomic paths including %d with 1 move\n" % (len(OAPS),nap_1))
  
  path.show_list(sys.stderr, "    ", OAPS, None, True,True)

  return OAPS
  # ----------------------------------------------------------------------
  
def show_group_stats(GRS):
  # Given a list {GRS} of group indices (possibly {None}) prints
  # some statistics from it.
  nix = len(GRS) # Count of indices given.  
  nng = 0 # Count of group indices {None}.
  igr_max = -1 # max non-None group
  for igr in GRS:
    if igr == None:
      nng += 1
    else:
      assert igr >= 0
      igr_max = max (igr_max, igr)
  GRCT = [ 0 ]*(igr_max + 1)
  for igr in GRS:
    if igr != None: GRCT[igr] += 1
  sys.stderr.write("seeing %d group indices in 0..%d plus %d with no group\n" % (nix-nng,igr_max,nng))
  if nng < nix:
    sys.stderr.write("counts per group: ")
    for igr in range(len(GRCT)):
      if GRCT[igr] != 0: sys.stderr.write(" %d:%d" % (igr,GRCT[igr]))
    sys.stderr.write("\n")
  
  return
  # ----------------------------------------------------------------------

def get_contour_blocks_and_contacts(OCRS, BCS_fill, mp_fill):
  # Given the list {OCRS} of contours of the {turkey} test piece, 
  # breaks the contours into sections, and makes blocks out of the
  # sections.  Returns a list {BCS} of those blocks.
  # 
  # Also returs a list {CTS} of the contacts between the choices of the
  # blocks of {BCS_fill} and these new blocks in {BCS}. Creates the
  # appropriate path-contact pointers.
  #
  # Also returns a list {OCRS_rest} with the contours of {OCRS}
  # that were not turned into blocks, for some reason.
  #
  # Sets the group index of the new blocks sequentially starting with
  # index {nbc = len(BCS_fill)}.
  #
  # Each move of the returned blocks will tend to be the projection on
  # the contour of a move of {BCS_fill} that is close enough to the
  # contour to make contact with it. This is done to maximize the number
  # of contacts between {BCS} and {BCS_fill}.
  
  wdf = move_parms.width(mp_fill);
  gap = 0.4*wdf # Length of gap to leave where contours are cut.
  nbc = len(BCS_fill)
  
  # Get the endpoints of all moves in {BCS_fil}, with group indices:
  EPS = get_all_move_endpoints(BCS_fill)
  
  BCS = []
  OCRS_rest = []
  for ocr in OCRS:
    Bi = path.bbox([ocr,])
    szy = Bi[1][1] - Bi[0][1]
    if szy > 5*wdf:
      # Tall contour - split into sections, make each section into a block:
      BCSi = make_blocks_from_tall_contour(ocr, EPS, BCS_fill, gap, wdf)
    else:
      # Short contour: split at top and bottom, make two choices.
      bci = make_block_from_short_contour(ocr, EPS, nbc, gap, wdf)
      BCSi = [bci,]
    nbc += len(BCSi)
    BCS += BCSi
    
  # Set the group indices:
  ngr = len(BCS_fill)
  for ibc in range(len(BCS)):
    bc = BCS[ibc]
    for ich in range(block.nchoices(bc)):
      oph = block.choice(bc,ich)
      path.set_group(oph, ngr)
    ngr += 1

  CTS = []
  szmin = 0.10*wdf
  rszmin = 0.0
  tol = 0.40*wdf  # Tolerance for overlaps, tilts, etc.

  ydir = None # We don't care about the order of the sides.
  for ibc in range(len(BCS)):
    bci = BCS[ibc]
    gri = path.get_group(block.choice(bci,0))
    for jbc in range(len(BCS_fill)):
      bcj = BCS_fill[jbc]
      grj = path.get_group(block.choice(bcj,0))
      CTSij = contact.from_blocks(bci, bcj, szmin, rszmin, tol, ydir)
      for ict in range(len(CTSij)):
        name = "C(%d:%d)%d" % (ibc,jbc,ict)
        contact.set_name(CTSij[ict], name)
      CTS += CTSij
      
  sys.stderr.write("found %d contour-filling contacts\n" % len(CTS))
  
  return BCS, CTS, OCRS_rest
  # ----------------------------------------------------------------------

def get_all_move_endpoints(BCS):
  # Given a list of filling blocks {BCS}, enumerates all the endpoints
  # of all the moves in all choices of all blocks of {BCS}, including
  # link paths that were included in those choices but not link paths
  # between blocks.
  # 
  # Returns a list of pairs {(p,igr)} where {p} is an move endpoint and
  # {igr} is the group index of the block to which the move belongs.
  EPS = set()
  for bc in BCS:
    for ich in range(block.nchoices(bc)):
      och = block.choice(bc, ich)
      igr = path.get_group(och)
      EPS.add((path.pini(och), igr))
      for imv in range(path.nelems(och)):
        omvi = path.elem(och, imv)
        EPS.add((move.pfin(omvi), igr))
      for ochr in och, path.rev(och):
        for olk in path.get_links(ochr):
          assert path.pfin(olk) == path.pini(ochr)
          for jmv in range(path.nelems(olk)):
            omvj = path.elem(olk, jmv)
            EPS.add((move.pini(omvj), igr))
  EPS = list(EPS)

  # Statistics for debugging:
  sys.stderr.write("collected %d filler endpoints total\n" % len(EPS))
  GRS = [ igr for p, igr in EPS ]
  show_group_stats(GRS)

  return EPS
  # ----------------------------------------------------------------------

def plot_figure_blocks(subfig, style,color):
  # A figure showing rasters partitioned into blocks or pre-blocks ("If").
  # 
  #   {subfig} == "genblocks",  plots the figure with the generic blocks, including links and contacts. 
  #   {subfig} == "genchoices", plots the four alternatives of a selected block from the generic blocks figure. 
  #   {subfig} == "gengraph",   plots the contact graphs of the generic blocks. 
  #   {subfig} == "ctrblocks",  plots the generic blocks with contours made into blocks. 
  #   {subfig} == "ctrchoices", plots the choices of a block of "ctrblocks". 
  #
  # The "ctrblocks" option has hooks in the filler blocks, instead of link paths.
  #
  # Returns the canvas {c} with the figure.

  # Get the contours, fillers, links, contacts, and graph of the "turkey" part:
  mp_cont, mp_fill, mp_link, mp_jump = make_move_parms(style);
  OCRS,OPHS,OLKS,CTS,VGS,EGS = paper_example_B.make_turkey(mp_cont, mp_fill,mp_link,mp_jump) 

  # Separate the rasters by scan-line:
  xdir = (1,0)
  ydir = (0,1)
  ystep, yphase = raster.get_spacing_and_phase(OPHS, xdir, ydir)
  SCS = raster.separate_by_scanline(OPHS, xdir, ydir, ystep, yphase)
  nsc = len(SCS) # Number of scan-lines.
  
  # Create the blocks to be plotted, as a list {BCSS} of lists of
  # blocks. For the sub-figures with canonical blocks, each element of
  # {BCSS} is a list of the blocks of one of the bands defined by the
  # cut-lines returned by {get_hotfill_cut_lines}. For the sub-figures
  # with generic blocks, each element of {BCSS} is a list with a single
  # block.

  # Get the blocks, canonical or generic, possibly with links and contacts or cut-lines,
  # comprising all rasters; and pick up a list {HBCS} of blocks to be highlighted or expanded;

  dp = (0,0)

  if subfig == "genblocks" or subfig == "genchoices" or subfig == "gengraph":
    # Create generic blocks, sub-rivers but not defined by cut-lines:
    full = True
    hooks = False
    BCS = get_generic_blocks(OPHS, SCS, full, hooks, mp_jump)
    # Select a block to highlight in the generic blocks figure:
    IBCHS = [ 13, ]
  elif subfig == "ctrblocks" or subfig == "ctrchoices":
    # Blocks are the same as in "genblocks" but contours are
    # chopped and added to blocks, with extra contacts:
    
    full = True
    hooks = True
    BCS_fill = get_generic_blocks(OPHS, SCS, full, hooks, mp_jump)
    # Get endpoints of link moves of {BCS_fill}:
    BCS_ctr, CTS_ctr, OCRS = get_contour_blocks_and_contacts(OCRS, BCS_fill, mp_fill)
    BCS = BCS_fill + BCS_ctr
    CTS = CTS + CTS_ctr
    # Select a block to highlight in the generic blocks figure:
    IBCHS = [ 13, 18 ]
  else:
    assert False, ("invalid subfig = %s" % subfig)

  nbc = len(BCS) # Number of blocks.
  BCHS = [ BCS[ibc] for ibc in IBCHS ]

  # Paranoia:
  for ibc in range(len(BCS)):
    bc = BCS[ibc]
    for ich in range(block.nchoices(bc)):
      och = block.choice(bc,ich)
      assert path.get_group(och) == ibc
  
  # Determine the figure's bounding box {Bfig} and the cut-line X range {cuxlo,cuxhi}:
  if subfig == "genchoices" or subfig == "ctrchoices":

    # Bounding box is for the blocks {bcH0,bcH1} only:
    choiceGap = 3.0 # X and Y gap between choices.
    DPS, Bfig = get_block_choice_disps_figbox(BCHS, choiceGap)

    wdf = style['wd_fill']
    mrg = (1.0*wdf, 1.0*wdf) # To account for sausage overflow and links.
    Bfig = rn.box_expand(Bfig, mrg, mrg)
    cuxlo = None; cuxhi = None

  else:

    # Bounding box is for the whole slice including contours.
    # Assume that links are contained in that box.
    B = path.bbox(OCRS+OPHS)
    Bfig, cuxlo, cuxhi = get_turkey_figure_bbox(OCRS, OPHS, OLKS, style)

  autoscale = True
  c = make_figure_canvas(Bfig, autoscale, style,color)

  if subfig == "genblocks" or subfig == "gengraph" or subfig == "ctrblocks":
  
    plot_figure_blocks_whole(c, subfig, OCRS, OPHS, CTS, BCS, BCHS, dp, cuxlo,cuxhi, nsc,ystep,yphase, style,color)
    
  elif subfig == "genchoices" or subfig == "ctrchoices":
  
    plot_figure_block_choices(c, subfig, BCHS, dp, DPS, ystep,yphase, style,color)
    
  else:
    assert False, ("invalid subfig %s" % subfig)

  return c
  # ----------------------------------------------------------------------

def get_block_choice_disps_figbox(BCS, gap):
  # Given a list {BXS} of blocks, returns a list of lists of displacements {DPS}, 
  # and a plot bounding box {Bfig}.
  #
  # The results are such that choice {ich} of block {ibc} should be
  # plotted with a displacement {DPS[ibc][ich]}
  
  nbc = len(BCS)
  BXS = [ None ] * nbc
  DPS = [ None ] * nbc
  
  szFig = [0, 0] # Upper corner of figure, with lower corner at (0,0).
  for ipass in 0, 1:
    # {ipass=0} compute the block bounding boxes {BXS} and global size {szFig} (incl {gap} around):
    # {ipass=1} computes displacements {DPS}
    if ipass == 1: dy = 0
    for ibc in range(nbc):
      bc = BCS[ibc]
      if bc != None:
        nch = block.nchoices(bc)
        if ipass == 0:
          B = block.bbox([bc,], True, False)
          BXS[ibc] = B
          DPS[ibc] = [ None ]*nch
        elif ipass == 1:
          B = BXS[ibc]
        else:
          assert False
        szB = rn.box_size(B)
        step = rn.add(szB, (gap,gap)) # Displacement between plots of choices.
        szbc = (szB[0] + (nch-1)*step[0], szB[1]) # Size of box of whole block.
        if ipass == 0:
          szFig[0] = max(szFig[0], szbc[0])
          szFig[1] = szFig[1] + (0 if ibc == 0 else gap) + szbc[1]
        elif ipass == 1:
          dx = (szFig[0] - szbc[0])/2
          for ich in range(nch):
            DPS[ibc][ich] = rn.sub((dx, dy), B[0])
            dx += step[0]
          dy += step[1]
        else:
          assert False
  
  Bfig = ((0,0), tuple(szFig))
  
  return DPS, Bfig
  # ----------------------------------------------------------------------

def plot_figure_blocks_whole(c, subfig, OCRS, OPHS, CTS, BCS, BCHS, dp, cuxlo,cuxhi, nsc,ystep,yphase, style,color):
  # Plots version {subfig} of the "blocks" figure, where {subfig is not 3.
  # The parameters {OCRS,OPHS,CTS} as as returned by {paper_example_B.make_turkey}.
  # The {BCS} parameter is a flat list of full-choice blocks buit from those rasters;
  # {BCHS} is a list of the blocks from {BCS} that are to be highlighted; {nsc} is the number of 
  # scan-lines, and {ystep,yphase} definde their positions;
  # {cuxlo,cuxhi} is the X-range of the cut-lines (excluding their labels).

  plot_contours(c, OCRS, dp, style, color['ghost'])
  
  nbc = len(BCS)

  # Plot the link paths for the generic blocks figure:
  if subfig == "genblocks" or subfig == "ctrblocks":
    # Plot the links beneath the blocks:
    bicolor = False
    for bc in BCS:
      plot_block_links(c, [bc,], dp, ystep,yphase, bicolor, style,color)
  elif subfig == "gengraph":
    # Neither cut-lines nor links:
    pass
  else:
    assert False

  # Plot the blocks:
  if subfig == "genblocks" or subfig == "ctrblocks" or subfig == "gengraph":

    # Plot the blocks, with uniform color, highlighting the selected ones:
    BCHS_set = set(BCHS)
    rwdf = style['rwd_fill']
    deco = False
    for ibc in range(nbc):
      bc = BCS[ibc]
      assert isinstance(bc, block.Block)
      clr = color['hifi'] if bc in BCHS_set else color['fill']
      if subfig == "gengraph": clr = ghostify_colors([clr,], color['ghost'])[0]
      nch = block.nchoices(bc)
      for ich in range(nch):
        oph = block.choice(bc, ich)
        assert path.get_group(oph) == ibc
        plot_paths(c, [oph,], dp, [clr,], rwdf,deco, style,color)
       
    # Label the highlighted blocks:
    for ibch in range(len(BCHS)):
      bch = BCHS[ibch]
      pbc = block.barycenter(bch)
      sys.stderr.write("highlight block %s" % block.get_name(bch))
      sys.stderr.write(" barycenter = %s\n" % str(pbc))
      bclab_dp = ( (+3.5,-1.5),  (-1.0,-0.5), (00.0,00.0) ) [ibch]
      if len(BCHS) == 1:
        bclab_tx = r"$B$"
      else:
        bclab_tx = r"$B_{%d}$" % ibch
      mfsize = style['fumafsize']
      plot_math_label(c, bclab_tx, pbc, bclab_dp, dp, mfsize, style)

  if subfig == "genblocks" or subfig == "ctrblocks":

    # Plot the inter-block contacts:
    interbc = True # Only inter-block contacts.
    shadow = False
    trim = (subfig == "genblocks")
    plot_contacts(c, CTS, dp, interbc, shadow,trim, style,color)

  elif subfig == "gengraph":

    # Compute the vertices of the contact graph as the barycenters of blocks.
    VGS = []
    for bc in BCS:
      pbar = block.barycenter(bc)
      VGS.append(pbar)
    # Compute the edges of the graph:
    EGS = []
    for ct in CTS:
      igr = contact_side_groups(ct)
      assert igr[0] != None and igr[1] != None, "contact side not in any block"
      if igr[0] != igr[1]:
        EGS.append(igr)
    # Plot the graph:
    plot_contact_graph(c, VGS, EGS, dp, style,color) 

  return
  # ----------------------------------------------------------------------
        
def plot_figure_block_choices(c, subfig, BCS, dp, DPS, ystep,yphase, style,color):
  # Plots version {subfig} of the "blocks" figure, where {subfig} is
  # "genchoices" or "ctrchoices", showing the choices of a specific
  # block with any attached links. The {BCS} paramerer a list of the blocks
  # whose choices are to be shown; {DPS[ibc][ich]} is the displacement
  # of the plot of choice {ich} of block {BCS[ibc]}; and {ystep,yphase}
  # defines the position of scan-lines.
    
  assert subfig == "genchoices" or  subfig == "ctrchoices"

  for ibc in range(len(BCS)):
    bc = BCS[ibc]
    nch = block.nchoices(bc)

    # Plot the block choices and links:
    clr = color['hifi']
    rwdf = style['rwd_fill']
    deco = True
    for ich in range(nch):
      ophi = block.choice(bc, ich)
      dpi = rn.add(dp, DPS[ibc][ich])

      # Plot the matter shadow of the choice and links:
      plot_trace_matter(c, [ophi,], dpi, style,color)
      for olk in path.get_links(ophi) + path.get_links(path.rev(ophi)):
        plot_trace_matter(c, [olk,], dpi, style,color)

      # Plot the choice and links:
      bicolor = False
      plot_path_links(c, [ophi,], dpi, ystep,yphase, bicolor, style,color)
      plot_paths(c, [ophi,], dpi, [clr,], rwdf,deco, style,color)
      # ??? Should show the attached contacts too ???
  return
  # ----------------------------------------------------------------------
  
# ######################################################################
# PLOTTING STUFF

def plot_paths(c, OPHS, dp, CLRS, rwd,deco, style,color):
  # Plot the traces and jumps of the paths {OPHS}, /without/ the matter traces.
  # The trace widths are reduced by {rwd} from what their {MoveParms} records say.
  # If {deco} is true, plots axes, arrows, and endpoints of all moves.
  # If {deco} is false, plots arrows and endpoints only of jumps.
  wd_axes = style['wd_axes']
  axes = deco
  dots = deco
  arrows = deco
  matter = False
  path.plot_standard(c, OPHS, dp, None, CLRS, rwd,wd_axes, axes, dots, arrows, matter)
  return
  # ----------------------------------------------------------------------

def plot_path_links(c, OPHS, dp, ystep,yphase, bicolor, style,color):
  # Plots onto {c} the link paths of the paths in {OPHS}. 
  for oph in OPHS:
    for ophr in oph, path.rev(oph):
      OLKS = path.get_links(ophr)
      plot_links(c, OLKS, dp, ystep,yphase, bicolor, style,color)
  return
  # ----------------------------------------------------------------------
      
def plot_links(c, OLKS, dp, ystep,yphase, bicolor, style,color):
  # Plots the paths in {OLKS} with the proper link style. 
  #
  # If {bicolor} is true, the color is chosen depending on the parity of
  # the scan-line index of the lowest endpoint of the link.
  # 
  rwd_link = style['rwd_link']
  wd_axes = style['wd_axes']
  axes = False
  dots = False
  arrows = False
  matter = False
  for ilk in range(len(OLKS)):
    olki = OLKS[ilk]
    # Make sure the link is directed up:
    if path.pini(olki)[1] > path.pfin(olki)[1]: olki = path.rev(olki)
    # Determine the scan-line index:
    isc = raster.point_scanline_index(path.pini(olki), (1,0), (0,1), ystep, yphase)
    clr = color['links'] if not bicolor else color['link0'] if isc % 2 == 0 else color['link1']
    path.plot_standard(c, [olki,], dp, None, [clr,], rwd_link, wd_axes, axes, dots, arrows, matter)
  return
  # ----------------------------------------------------------------------

def plot_block_links(c, BCS, dp, ystep,yphase, bicolor, style,color):
  # Plots onto {c} the link paths of the choices of the blocks in {BCS}.  The color is chosen
  # depending on the parity of the scan-line index of the lowest endpoint of the link.
  
  for bc in BCS:
    nch = block.nchoices(bc)
    for ich in range(nch):
      oph = block.choice(bc, ich)
      plot_path_links(c, [oph,], dp, ystep,yphase, bicolor, style,color)
  return
  # ----------------------------------------------------------------------
    
def plot_dot_on_path(c, p, dp, style,color):
  clr = color['dots']
  wd_edots = style['wd_edots']
  hacks.plot_line(c, p, p, dp, clr, wd_edots, None)
  return
  # ----------------------------------------------------------------------

def plot_path_endpoints(c, OPHS, dp, style,color):
  # Plots the endpoints of the paths in {OPHS}.
  for oph in OPHS:
    for p in path.endpoints(oph):
      plot_dot_on_path(c, p, dp, style,color)
  return
  # ----------------------------------------------------------------------
  
def plot_contacts(c, CTS, dp, interbc, shadow,trim, style,color):
  # Plots the contacts in {CTS} with the style used in the paper, with color {color['ctac']}.
  # If {interbc} is true, plots only contacts between paths that have different group indices.
  # If {shadow} is true, plots a white shadow unde the contacts.
  # If {trim}, shrinks the contacts to avoid interference with link paths.
   
  wd_shad = style['wd_shad']
  wd_ctac = style['wd_ctac']
  wd_link = style['rwd_link']*style['wd_fill']
  dashpat = (1.0*wd_ctac, 2.0*wd_ctac)
  ctlen_min = 2*wd_ctac # Min length of plotted contact line.
  if trim:
    ext_max = - (0.5*wd_link + wd_ctac) # To avoid links in some plots.
  else:
    ext_max = 0
  sz_tics = 0
  arrows_ct = False
  ncp = 0 # Interblock contacts plotted.
  for ct in CTS:
    igr = contact_side_groups(ct)
    plotit = True if not interbc else igr[0] != igr[1]
    if trim: 
      # Adjust the shrinking so as to not eliminate short contacts:
      ep = contact.endpoints(ct)
      ctlen = rn.dist(ep[0],ep[1])
      ext_ct = min(ext_max, ctlen - ctlen_min)
    else:
      ext_ct = 0
    if plotit:
      if shadow: 
        contact.plot_single(c, ct, dp, color['white'], None, ext_ct, wd=wd_ctac+wd_shad, sz_tic=sz_tics, arrow=arrows_ct)
      contact.plot_single(c, ct, dp, color['ctac'], dashpat, ext_ct, wd=wd_ctac, sz_tic=sz_tics, arrow=arrows_ct)
      ncp += 1
      
  sys.stderr.write("plotted %d out of %d contacts\n" % (ncp,len(CTS)))
  return
  # ----------------------------------------------------------------------

def plot_cut_lines(c, LNS, cuxlo,cuxhi,trim, ystep,yphase, dp, style,color):
  # The parameter {LNS} must be a list of pairs {(icu, t)}, as in the 
  # output of {get_essential_cut_lines}.
  # Plots the essential ({t=2}) and non-essential ({t=1}) cut-lines listed in {LNS}.
  # Ignores cut-lines with {t=0}.
  #
  # If {trim} is True, plots onlt the stubs of the lines on each side.
  #
  wd_cuts = style['wd_cuts'] 
  for icu, t in LNS:
    if t != 0:
      clr = color['cuess'] if t == 2 else color['cunon']
      dashpat = [ 0.100, 0.250 ] if t == 1 else None
      y = ystep*(icu-0.5) + yphase
      p = (cuxlo, y)
      q = (cuxhi, y)
      if trim:
        stublen = 1.25
        pf = rn.add(p, (stublen,0))
        hacks.plot_line(c, p, pf, dp, clr, wd_cuts, dashpat)
        qf = rn.sub(q, (stublen,0))
        hacks.plot_line(c, q, qf, dp, clr, wd_cuts, dashpat)
      else:
        hacks.plot_line(c, p, q, dp, clr, wd_cuts, dashpat)

      # Label the cut-line on both sides:
      fsize = style['cutlfsize']
      tx = r"\texttt{%d}" % icu
      tdp = (-1.7,-0.3) if icu > 9 else (-1.1,-0.3)
      pt = rn.add(p, tdp)
      hacks.plot_text(c, tx, pt, dp, fsize, None)
      tdq = (1,-0.3)
      qt = rn.add(q, tdq)
      hacks.plot_text(c, tx, qt, dp, fsize, None)

  return
  # ----------------------------------------------------------------------
    
def plot_trace_matter(c, OPHS, dp, style,color):
  # Plot the matter traces of all the paths in {OPHS}:
  clr = color['matter']
  jmp = False     # Plot traces only.
  wd_matter = 0   # Extra width of material footprint.
  rwd_matter = style['rwd_matter']
  dashed = False
  wd_dots = 0
  sz_arrows = 0
  for ph in OPHS: 
    path.plot_layer(c, ph, dp, jmp, clr, rwd_matter, wd_matter, dashed, wd_dots, sz_arrows)
  return
  # ----------------------------------------------------------------------
     
def plot_jumps(c, JMPS, dp, clr, style,color):
  # Plot the jumps in {JMPS} with color {clr}:
  wd_axis = style['wd_axes']
  dashed = True
  wd_dots = 2.5*wd_axis
  sz_arrow = 6*wd_axis
  for jmp in JMPS: 
    move.show(sys.stderr, "  jmp = ", jmp, "\n", 0)
    move.plot_layer(c, jmp, dp, clr, wd_axis, dashed, wd_dots, sz_arrow)
  return
  # ----------------------------------------------------------------------
 
def plot_contours(c, OCRS, dp, style, clr):
  # Plot the contour paths {OCRS} with color {clr}:

  rwd_cont = style['rwd_cont']
  wd_axes = style['wd_axes']
  axes = False
  dots = False
  arrows = False
  matter = False
  path.plot_standard(c, OCRS, dp, None, [clr,], rwd_cont, wd_axes, axes, dots, arrows, matter)
  return
  # ----------------------------------------------------------------------

def plot_math_label(c, tx, qP, qD, dp, mfsize, style):
  # Plots the text {tx} at point {dp+qP+qD} with font size {mfsize}.
  # If {dp} is {None}, assumes {(0,0)}
  
  q = rn.add(qP, qD)
  hacks.plot_text(c, tx, q, dp, mfsize, None)
  return 
  # ----------------------------------------------------------------------
    
def plot_contact_graph(c, VGS, EGS, dp, style,color):
  # Plots the contact graph given the vertices {VGS} (a list of points) and the edges {EGS}
  # (a list of pairs of indices into {VGS}).
  
  # Plot edges:
  wd_edges = style['wd_edges']
  clr_edges = color['edges']
  for e in EGS:
    p = VGS[e[0]]
    q = VGS[e[1]]
    hacks.plot_line(c, p, q, dp, clr_edges, wd_edges, None)

  # Plot vertices:
  wd_verts = style['wd_verts'] # Diameter of vertices.
  clr_verts = color['verts']
  for p in VGS:
    hacks.plot_line(c, p, p, dp, clr_verts, wd_verts, None)
  return
  # ----------------------------------------------------------------------

def plot_point_label(c, lab, p, td, dp, pfsize, style,color):
  # Places label {lab} next with the lower left corner at the point {p}
  # displaced by {td} and {dp} (if not {None}), with font size {pfsize}.  Also plots a white shadow all around, with width
  # proportional to the font size.
  
  p = rn.add(p, td)
  tad = 0.002*pfsize
  for kx in range(5):
    for ky in range(5):
      dx = kx - 2; dy = ky - 2;
      if dx != 0 or dy != 0:
        pk = rn.add(p, (dx*tad, dy*tad))
        txk = r"\textbf{\textsf{%s}}" % lab
        hacks.plot_text(c, txk, pk, dp, pfsize, 'white')
  tx = r"\textbf{\textsf{%s}}" % lab
  hacks.plot_text(c, tx, p, dp, pfsize, None)
  return
  # ----------------------------------------------------------------------

def plot_path_endpoint_labels(c, oph, labA, tdA, labB, tdB, dp, pfsize, style,color):
  # Places labels {labA} and {labB} next to the start and end of 
  # path {oph}, displaced by {tdA} and {tdB}, respectively, with fonts size {pfsize}.
  plot_point_label(c, labA, path.pini(oph), tdA, dp, pfsize, style,color)
  plot_point_label(c, labB, path.pfin(oph), tdB, dp, pfsize, style,color)
  return
  # ----------------------------------------------------------------------

# ######################################################################
# UTILITY FUNCTIONS

def contact_side_groups(ct):
  # Returns the group indices of the paths that contain the two sides of {ct}.
  # Uses {contact.get_side_paths} and {path.get_group}. 
  # Each side must have at least one attached path, and all groups attached 
  # to the same side must have the same group index.
  
  debug = False
  if debug: contact.show(sys.stderr, "  ct = ", ct, "\n", 0)
  igr = [-1, -1] # Group indices on the two sides of contact:
  for isd in range(2):
    for phij, drij, imvij in contact.get_side_paths(ct, isd):
      igrij = path.get_group(phij)
      assert igrij == None or igrij >= 0
      if debug: path.show(sys.stderr, ("    side %d @ " % isd), phij, (" igr = %s\n" % str(igrij)), True, 0,0,0)
      if igr[isd] == -1:
        igr[isd] = igrij
      else:
        assert igr[isd] == igrij, ("side %d of contact on different groups" % isd)
    assert igr[isd] != -1, "contact has no side paths"
  return tuple(igr)
  # ----------------------------------------------------------------------

def tex_float_fmt(v, maxprec):
  # Formats float value {v} with at most {maxprec} decimal fraction digits.
  vv = v
  prec = 0
  while prec < maxprec and floor(vv) != vv:
    vv = vv*10; prec += 1
  s = "%.*f" % (prec, v)
  return s
  # ----------------------------------------------------------------------

def write_text_parms(fname, style):
  # Writes to file {fname} the defintion of TeX macros
  # with the printer and job parameters used in the figures.
  
  sys.stderr.write("writing file %s ...\n" % fname)
  wr = open(fname, "w")
  
  acc =       style['acc']       # Acceleration/deceleration for moves and jumps (mm/s^2).
  csp_trace = style['csp_trace'] # Cruise speed for traces (mm/s).
  csp_jump =  style['csp_jump']  # Cruise speed for jumps (mm/s).
  udp =       style['udp']       # Trace/jump transition penalty (s).
  
  wd_cont = style['wd_cont']  # Contour trace width (mm).
  wd_fill = style['wd_fill']  # Fill trace width (mm).
  
  mname = "FigParms" # Prefix for parameter names
  
  wr.write((r"\newcommand{\%sContourWidth}{%s}" + "\n") % (mname, tex_float_fmt(wd_cont,2)))
  wr.write((r"\newcommand{\%sFillingWidth}{%s}" + "\n") % (mname, tex_float_fmt(wd_fill,2)))
  wr.write((r"\newcommand{\%sAcceleration}{%s}" + "\n") % (mname, tex_float_fmt(acc,2)))
  wr.write((r"\newcommand{\%sTraceSpeed}{%s}" + "\n")   % (mname, tex_float_fmt(csp_trace,2)))
  wr.write((r"\newcommand{\%sJumpSpeed}{%s}" + "\n")    % (mname, tex_float_fmt(csp_jump,2)))
  wr.write((r"\newcommand{\%sJumpPenalty}{%s}" + "\n")  % (mname, tex_float_fmt(udp,2)))

  wr.close()
  return
  # ----------------------------------------------------------------------

def make_move_parms(style):
  # Returns the {MoveParms} records to use.
  
  acc =       style['acc']       # Acceleration/deceleration for moves and jumps (mm/s^2).
  csp_trace = style['csp_trace'] # Cruise speed for traces (mm/s).
  csp_jump =  style['csp_jump']  # Cruise speed for jumps (mm/s).
  udp =       style['udp']       # Trace/jump transition penalty (s).
  
  wd_cont = style['wd_cont']  # Contour trace width (mm).
  wd_fill = style['wd_fill']  # Fill trace width (mm).
  
  mp_cont = move_parms.make(wd_cont, acc, csp_trace, 0.0)
  mp_fill = move_parms.make(wd_fill, acc, csp_trace, 0.0)
  mp_link = mp_fill
  mp_jump = move_parms.make(0.0,     acc, csp_jump,  udp)

  sys.stderr.write("printer parameters:\n")
  move_parms.show(sys.stderr, "contours = { ", mp_cont, " }\n")
  move_parms.show(sys.stderr, "filling =  { ", mp_fill, " }\n")
  move_parms.show(sys.stderr, "links   =  { ", mp_link, " }\n")
  move_parms.show(sys.stderr, "jumps =    { ", mp_jump, " }\n")

  return mp_cont, mp_fill, mp_link, mp_jump
  # ----------------------------------------------------------------------

def make_style_dict():
  # Returns a Python dict with the plot dimension parameters ({wd_fill},
  # {wd_axes}, etc) to be used in figures of the paper.
  
  # They should be the same as used in the tests section of the paper,
  # exccept that the trace width will be 1.0 (instead of 0.4) to 
  # make the figures more readable.
  
  wd_cont =  0.75 # Contour trace width (mm).
  wd_fill =  1.00 # Fill trace width (mm).

  rwd_path = 0.80 # Plot trace sausages with relative width {rwd_path} on isolated path plots.
  rwd_fill = 0.50 # Plot fill trace sausages with relative width {rwd_fill} on turkey plots.

  wd_fisa = rwd_fill*wd_fill # Width of plotted sausages in fillin.

  wd_axes = 0.10*wd_fill
  wd_ctac = 0.20*wd_fill
  wd_cuts = 0.10*wd_fill

  style = {
    'detafsize':  30,             # Font size for point labels in detail and notation figures.
    'mathfsize':  36,             # Font size for math labels in notation figures.
    'fullfsize':  48,             # Font size for point labels in full-size slice figures.
    'fumafsize':  58,             # Font size for math labels in notation figures.
    'cutlfsize':  36,             # Font size for labels of cut-lines on full-size slice figure.
    'tbixfsize':  40,             # Font size for row and column indices on tableau figure.
    
    'cuext':                2.0,   # Extra left and right overhang of cut-lines, excl. label (mm).
    'culab':                2.5,   # Estimated width of cut-line labels, including space (mm).

    'acc':                 3000,   # Acceleration/deceleration for moves and jumps (mm/s^2).
    'csp_trace':             40,   # Cruise speed for traces (mm/s).
    'csp_jump':             130,   # Cruise speed for jumps (mm/s).
    'udp':                 0.05,   # Trace/jump transition penalty (s).
  
    'wd_cont':          wd_cont,   # Contour trace width (mm).
    'wd_fill':          wd_fill,   # Fill trace width (mm).

    'wd_axes':          wd_axes,   # Width of jumps and axis lines.
    'wd_ctac':          wd_ctac,   # Width of contact lines.
    'wd_shad':          wd_ctac,   # Width of white shadow under contacts.

    'wd_cuts':          wd_cuts,   # Width of cut-lines.
    'wd_edots':    0.75*wd_fisa,   # Diameter of black dots at end of paths.
    
    'wd_edges':    0.20*wd_fill,   # Width of graph edges.
    'wd_verts':    1.20*wd_fisa,   # Diameter of graph vertices.
    
    'rwd_path':        rwd_path,   # Relative width of fill trace sausages in isolated path plots.
    'rwd_fill':        rwd_fill,   # Relative width of fill trace sausages in turkey plots.
    'rwd_cont':        rwd_fill,   # Relative width of contour trace sausages in turkey plots.
    'rwd_link':   0.67*rwd_fill,   # Relative width of link trace sausages in turkey plots.
    'rwd_matter':          1.13,   # Relative width of material footprint in plots.

    'wd_mfig':             34.0,   # Standard width of figures with math formulas.
  }
  return style
  # ----------------------------------------------------------------------
  
def make_color_dict():
  # Returns a Python dict with the standard color scheme for paper figures.
  color = {
    'white':  pyx.color.rgb( 1.000, 1.000, 1.000 ),  # Color for invisible frame and label shadow.
    'black':  pyx.color.rgb( 0.000, 0.000, 0.000 ),  # Default color.
    'matter': pyx.color.rgb( 0.900, 0.870, 0.850 ),  # Color for estimated material footprints.
    'fill':   pyx.color.rgb( 0.000, 0.900, 0.050 ),  # Color for relevant fill traces.
    'fil0':   pyx.color.rgb( 0.200, 0.900, 0.050 ),  # Color for fill traces in even special paths.
    'fil1':   pyx.color.rgb( 0.000, 0.600, 0.250 ),  # Color for fill traces in odd special paths.
    'hifi':   pyx.color.rgb( 0.000, 0.700, 1.000 ),  # Color for highlighted fill traces.
    'cont':   pyx.color.rgb( 0.600, 0.700, 0.800 ),  # Color for contours, when somewhat relevant.
    'ghost':  pyx.color.rgb( 0.850, 0.850, 0.850 ),  # Color for non-relevant traces.
    'ctac':   pyx.color.rgb( 1.000, 0.200, 0.000 ),  # Color for contacts.
    'dots':   pyx.color.rgb( 0.000, 0.000, 0.000 ),  # Color of dots (maybe not used).
    'links':  pyx.color.rgb( 0.450, 0.050, 1.000 ),  # Color of links without distinction.
    'link0':  pyx.color.rgb( 0.950, 0.750, 0.000 ),  # Color of links above even scan-lines.
    'link1':  pyx.color.rgb( 0.450, 0.000, 1.000 ),  # Color of links above odd scan-line.
    'edges':  pyx.color.rgb( 0.000, 0.000, 0.000 ),  # Color or edges of contact graph
    'verts':  pyx.color.rgb( 0.000, 0.000, 0.000 ),  # Color of vertices of contact graph.
    'jump':   pyx.color.rgb( 0.000, 0.000, 0.000 ),  # Color of jumps (maybe not used).
    'jphi':   pyx.color.rgb( 0.450, 0.000, 1.000 ),  # Color of highlighted jumps.
    'cuess':  pyx.color.rgb( 1.000, 0.200, 1.000 ),  # Color of essential cut-lines.
    'cunon':  pyx.color.rgb( 0.400, 0.400, 1.000 ),  # Color of non-essential cut-lines.
    
    'Ytrace': 0.600  # Luminosity of colors in trace/path/block palette.
  }
  return color
  # ----------------------------------------------------------------------
  
def make_figure_canvas(Bfig, autoscale, style,color):
  # Creates the canvas for a figure whose contents has bounding box {Bfig}.
  
  # Extra extra margin:
  mrg = (0.2,0.2)
  Bplot = rn.box_expand(Bfig, mrg, mrg)

  dp = None
  scale = None if autoscale else 0.5
  c, szx, szy = hacks.make_canvas(Bplot, dp, scale, False, False, 1, 1)

  # # Plot an invisible frame to keep the figure sizes uniform:
  # wd_frame = 0.5*style['wd_axes']
  # clr_frame = color['white'] # Invisible.
  # # clr_frame = color['black'] # Just to check.
  # hacks.plot_frame(c, clr_frame, wd_frame, dp, Bfig, wd_frame/2)
  
  return c
  # ----------------------------------------------------------------------
  
def get_turkey_figure_bbox(OCRS, OPHS, OLKS, style):
  # Given the lists of contours {OCRS}, raster paths {OPHS}, and link
  # paths {OLKS} of the "turkey" slice, returns a bounding box suitable
  # for {make_canvas}.
  #
  # The box includes all the move endpoints in those paths, plus some extra space on the sides for 
  # the cut-lines that may be inserted in it, and some extra space to account for matter and sausage width.
  #
  # Also returns the abscissas {cuxlo} and {cuxhi} of the endpoints of cut-lines (excluding the labels).

  # Get the bounding box {B} of all contents:
  B = path.bbox(OPHS)
  Bsize = rn.sub(B[1], B[0])
  sys.stderr.write("bounding box (excl. contours) = %6.2f x %6.2f mm\n" % (Bsize[0], Bsize[1]))
  if len(OCRS) != 0: B = rn.box_join(B, path.bbox(OCRS))
  if len(OLKS) != 0: B = rn.box_join(B, path.bbox(OLKS))
  # Assume that the contacts are automatically included in the box.

  # Compute X-range of cut-lines:
  cuxlo = B[0][0] - style['cuext']
  cuxhi = B[1][0] + style['cuext']

  # Expand box for cut-lines and their labels:
  dxcut = style['cuext'] + style['culab']
  Bfig = B
  Bfig = rn.box_expand(Bfig, (dxcut,0.3), (dxcut,0.3))

  # Extra margin for sausage widths:
  wdf = style['wd_fill']
  mrg = (wdf, wdf)
  Bfig = rn.box_expand(Bfig, mrg, mrg)

  return Bfig, cuxlo, cuxhi
  # ----------------------------------------------------------------------

def ghostify_colors(CLRS, cgh):
  # Mixes a lot of {cgh} into each color of {CLRS}.
  CLRS_new = []
  f = 0.3
  for clr in CLRS:
    R = f*clr.r + (1-f)*cgh.r
    G = f*clr.g + (1-f)*cgh.g
    B = f*clr.b + (1-f)*cgh.b
    clr_new = pyx.color.rgb(R,G,B)
    CLRS_new.append(clr_new)
  return CLRS_new
  # ----------------------------------------------------------------------
  
def widen_box_for_math_figure(B, style):
  wd_mfig = style['wd_mfig']
  xmrg = (wd_mfig - (B[1][0] - B[0][0]))/2
  assert xmrg >= 0, "math-containing box is wider than standard"
  B = rn.box_expand(B, (xmrg,0), (xmrg,0))
  return B
  # ----------------------------------------------------------------------
 
def write_tex_params(fig, subfig, OPHS, CTS):
  # Writes a file "out/paper_figures_B_TST_{fig}_{subfig}.tex" 
  # with TeX macros with the total fabtimes of the paths in {OPHS}.
  # 
  # If {CTS} is not empty or {None}, then {OPHS} mustbe a single path
  # {[ophs]}, and it writes also TeX macros with the cooling times of
  # the contacts in {CTS}, must be closed by {oph}.
  
  fname = "tests/out/paper_figures_B_TST_p_" + fig + "_" + subfig + ".tex"
  sys.stderr.write("writing %s ...\n" % fname)
  wr = open(fname, "w")
  
  # Replace digits 0-9 by letters A-J in {subfig} for use in TeX macro names:
  d2a = str.maketrans('0123456789', 'ABCDEFGHIJ')
  fig_tex = fig.translate(d2a)
  subfig_tex = subfig.translate(d2a)
  mname = r"F%sS%s" % (fig_tex, subfig_tex)
  
  # Fabrication times:
  tfab_tot = 0 # Total fabtime of paths in {OPHS}
  tfab_trc = 0 # Total fabtime of traces in {OPHS}
  for oph in OPHS:
    tfab_tot += path.fabtime(oph)
    OSQS, OJMS = path.split_at_jumps(oph)
    for osq in OSQS:
      tfab_trc += path.fabtime(osq)
  wr.write((r"\newcommand{\%sFabTime}{%s}" + "\n") % (mname, tex_float_fmt(tfab_tot,2)))
  wr.write((r"\newcommand{\%sExtTime}{%s}" + "\n") % (mname, tex_float_fmt(tfab_trc,2)))
  wr.write((r"\newcommand{\%sAirTime}{%s}" + "\n") % (mname, tex_float_fmt(tfab_tot - tfab_trc,2)))
    
  if CTS != None and len(CTS) > 0:
    # Compute and show the cooling times:
    assert len(OPHS) == 1
    oph = OPHS[0]
    for ict in range(len(CTS)):
      ct = CTS[ict]
      tcool = contact.tcool(oph, ct)
      if tcool == None:
        ctna = contact.get_name(ct)
        wr.write(" contact %s is not closed\n", ctna)
      else:
        ict_tex = (str(ict)).translate(d2a)
        wr.write((r"\newcommand{\%sContact%stcool}{%s}" + "\n") % (mname, ict_tex, tex_float_fmt(tcool,2)))
  wr.close()
  return
  # ----------------------------------------------------------------------
