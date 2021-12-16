import path
import rn
import move
import hacks
import move_parms
import contact
import pyx
import sys
import math as m
import txt_write

# --------------------------------------------------------------- #

def plot_plots(ph, OPHS, BPHS, CTScold, tag, alg, outfolder):
  fname = outfolder + "/" + alg + "-" + tag + '-'

  style = make_style_dict(OPHS)
  color = make_color_dict()

  parms_fname = outfolder + "/" + alg + "_" + tag + "_printer_parms.tex"
  write_text_parms(parms_fname, style)

  tex_fname = outfolder + "/" + alg + "_" + tag + "_times.tex"
  # Write TeX file with fabtimes and cooling times:
  write_tex_params(tex_fname, "cold", "path", [ph,], CTScold)

  c, Bfig = plot_figure(style, color, OPHS, ph, CTScold)
  
  if c != None:
    hacks.write_plot(c, fname + 'fp', 6)

  if Bfig != None and 'hotfill' in alg:
    c = plot_figure_2(style, color, BPHS, CTScold, Bfig)

    if c != None:
      hacks.write_plot(c, fname + 'bp', 6)

  return

# --------------------------------------------------------------- #

def make_style_dict(BPHS):
  # Returns a Python dict with the plot dimension parameters ({wd_fill},
  # {wd_axes}, etc) to be used in figures of the paper.
  
  # They should be the same as used in the tests section of the paper,
  # exccept that the trace width will be 1.0 (instead of 0.4) to 
  # make the figures more readable.
  
  wd_cont =  0.75 # Contour trace width (mm).
  mp_fill = move.parameters(path.elem(BPHS[0],0))
  wd_fill = move_parms.width(mp_fill)

  rwd_path = 0.8 # Plot trace sausages with relative width {rwd_path} on isolated path plots.
  rwd_fill = 0.5 # Plot fill trace sausages with relative width {rwd_fill} on turkey plots.

  wd_fisa = rwd_fill*wd_fill # Width of plotted sausages in fillin.

  wd_axes = 0.05*wd_fill
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

# --------------------------------------------------------------- #

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
    'jphi':   pyx.color.rgb( 0.450, 0.450, 0.450 ),  # Color of highlighted jumps.
    'cuess':  pyx.color.rgb( 1.000, 0.200, 1.000 ),  # Color of essential cut-lines.
    'cunon':  pyx.color.rgb( 0.400, 0.400, 1.000 ),  # Color of non-essential cut-lines.
    
    'Ytrace': 0.600  # Luminosity of colors in trace/path/block palette.
  }
  return color
  # ----------------------------------------------------------------------
 
# --------------------------------------------------------------- #

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

# --------------------------------------------------------------- #

def write_txt_file(fname, OPHS):
  wr = open(fname, "w")
  Z = 10.0
  angle = 0
  shift = (0,0)
  for oph in OPHS: path.set_group(oph,0)
  txt_write.write(wr, [], OPHS, Z, angle, shift)
  wr.close()

  # ----------------------------------------------------------------------

# --------------------------------------------------------------- #

def plot_figure(style, color, OPHS, ph, CTScold):  
  dp = (0,0)
  
  Bfig, cuxlo, cuxhi = get_figure_bbox(OPHS, style)
  autoscale = True
  c = make_figure_canvas(Bfig, autoscale, style,color)

  # Plot coldest contact(s).
  interbc = False
  shadow = True
  trim = False
  plot_contacts(c, CTScold, dp, interbc, shadow, trim, style, color)

  # SPlit at jumps for coloring:
  OSKS, JMPS = path.split_at_jumps(ph)
  Ylum = color['Ytrace']
  CLRS = hacks.trace_colors(len(OSKS), Ylum)

  # Plot trace components:
  rwdf = style['rwd_fill']
  deco = False
  plot_paths(c, OSKS, dp, CLRS, rwdf,deco, style,color)

  # Plot the jumps:
  clr_jp = color['jump']
  plot_jumps(c, JMPS, dp, clr_jp, style, color)
  
  # Plot the path endpoints:
  plot_path_endpoints(c, [ph,], dp, style,color)

  contact.show_times(sys.stderr, [ph,], CTScold)
    
  # Label path endpoints:
  tdA = (+0.8, -1.4)
  tdB = (-2.0, -0.4)
  pfsize = style['fullfsize']
  #plot_path_endpoint_labels(c, ph, "A", tdA, "B", tdB, dp, pfsize, style,color)
  
  if c == None:
    Bfig = None

  return c, Bfig
  # ----------------------------------------------------------------------

# --------------------------------------------------------------- #

def get_figure_bbox(OPHS, style):
  # Get the bounding box {B} of all contents:
  B = path.bbox(OPHS)
  Bsize = rn.sub(B[1], B[0])
  sys.stderr.write("bounding box (excl. contours) = %6.2f x %6.2f mm\n" % (Bsize[0], Bsize[1]))
  # Assume that the contacts are automatically included in the box.

  # Compute X-range of cut-lines:
  cuxlo = B[0][0] 
  cuxhi = B[1][0]

  # Expand box for cut-lines and their labels:
  Bfig = B

  # Extra margin for sausage widths:
  wdf = style['wd_fill']
  mrg = (wdf, wdf)
  Bfig = rn.box_expand(Bfig, mrg, mrg)

  return Bfig, cuxlo, cuxhi
  # ----------------------------------------------------------------------

# --------------------------------------------------------------- #

def make_figure_canvas(Bfig, autoscale, style, color):
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

# --------------------------------------------------------------- #

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

# --------------------------------------------------------------- #

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

# --------------------------------------------------------------- #

def plot_paths(c, OPHS, dp, CLRS, rwd,deco, style, color):
  # Plot the traces and jumps of the paths {OPHS}, /without/ the matter traces.
  # The trace widths are reduced by {rwd} from what their {MoveParms} records say.
  # If {deco} is true, plots axes, arrows, and endpoints of all moves.
  # If {deco} is false, plots arrows and endpoints only of jumps.
  wd_axes = style['wd_axes']
  axes = False
  dots = deco
  arrows = deco
  matter = False
  path.plot_standard(c, OPHS, dp, None, CLRS, rwd, wd_axes, axes, dots, arrows, matter)
  return
  # ----------------------------------------------------------------------

# --------------------------------------------------------------- #

def plot_jumps(c, JMPS, dp, clr, style,color):
  # Plot the jumps in {JMPS} with color {clr}:
  wd_axis = 2.5*style['wd_axes']
  dashed = True
  wd_dots  = 3.0*wd_axis
  sz_arrow = 5.0*wd_axis
  
  for jmp in JMPS:
    move.show(sys.stderr, "  jmp = ", jmp, "\n", 0)
    move.plot_layer(c, jmp, dp, clr, wd_axis, dashed, wd_dots, sz_arrow)

  return
  # ----------------------------------------------------------------------

# --------------------------------------------------------------- #

def plot_path_endpoints(c, OPHS, dp, style,color):
  # Plots the endpoints of the paths in {OPHS}.
  for oph in OPHS:
    for p in path.endpoints(oph):
      plot_dot_on_path(c, p, dp, style,color)
  return
  # ----------------------------------------------------------------------
  
# --------------------------------------------------------------- #

def plot_dot_on_path(c, p, dp, style,color):
  clr = color['dots']
  wd_edots = style['wd_edots']
  hacks.plot_line(c, p, p, dp, clr, wd_edots, None)
  return
  # ----------------------------------------------------------------------

# --------------------------------------------------------------- #

def write_tex_params(fname, fig, subfig, OPHS, CTS):
  # Writes a file "out/paper_figures_B_TST_{fig}_{subfig}.tex" 
  # with TeX macros with the total fabtimes of the paths in {OPHS}.
  # 
  # If {CTS} is not empty or {None}, then {OPHS} mustbe a single path
  # {[ophs]}, and it writes also TeX macros with the cooling times of
  # the contacts in {CTS}, must be closed by {oph}.
  
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

# --------------------------------------------------------------- #

def tex_float_fmt(v, maxprec):
  # Formats float value {v} with at most {maxprec} decimal fraction digits.
  vv = v
  prec = 0
  while prec < maxprec and m.floor(vv) != vv:
    vv = vv*10; prec += 1
  s = "%.*f" % (prec, v)
  return s
  # ----------------------------------------------------------------------

# --------------------------------------------------------------- #

def plot_figure_2(style, color, BPHS, CTScold, Bfig):
  dp = (0,0)
  
  autoscale = True
  c = make_figure_canvas(Bfig, autoscale, style,color)

  interbc = False
  shadow = False
  trim = False
  plot_contacts(c, CTScold, dp, interbc, shadow, trim, style, color)

  PPHS = []
  for iph in range(len(BPHS)):
    if iph % 2 == 0:
      PPHS.append(BPHS[iph])

  nbp = len(PPHS) # Number of bandpaths.
  CLRS = hacks.trace_colors(nbp, None)

  # Plot trace components:
  rwdf = style['rwd_fill']
  deco = True
  plot_paths(c, PPHS, dp, CLRS, rwdf, deco, style, color)
  
  return c

# --------------------------------------------------------------- #