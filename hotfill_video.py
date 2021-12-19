import hotfill
import paper_example_B
import raster
import path 
import move_parms
import contact
import time
import sys
import main_plots
import hacks
import move
import pyx

##
# Some arbitrary dynamics parameters:
acc = 3000      # Acceleration/deceleration for moves and jumps (mm/s^2).
csp_trace = 40  # Cruise speed for traces (mm/s).
csp_jump = 130  # Cruise speed for jumps (mm/s).
udp = 0.05      # Trace/jump transition penalty (s).
outfolder = ""

# --------------------------------------------------------------- #
# --------------------------------------------------------------- #
# --------------------------------------------------------------- #
# --------------------------------------------------------------- #
# --------------------------------------------------------------- #
# --------------------------------------------------------------- #

def plot_OPHS(c, OPHS, rwd, dp):
  clr = pyx.color.rgb(0.850, 0.850, 0.850)
  
  for oph in OPHS:
    omv = path.elem(oph, 0)
    move.plot_layer(c, omv, dp, clr, rwd, dashed=False, wd_dots=0, sz_arrow=None)

  return

# --------------------------------------------------------------- #

def plot_lines(c, OMVS, rwd, dp, CLRS, style, color):
  index_clr = 0
  clr = pyx.color.rgb(0.250, 0.360, 0.850)

  JMPS = []
  
  for omv in OMVS:
    if move.is_jump(omv):
      index_clr = index_clr + 1
      JMPS.append(omv)
      
    else:
      clr = CLRS[index_clr]
      move.plot_layer(c, omv, dp, clr, rwd, dashed=False, wd_dots=0, sz_arrow=None)
  
  main_plots.plot_jumps(c, JMPS, dp, color['jump'], style, color)

  return

# --------------------------------------------------------------- #

def plot_fp(OPHS, ph, tag):
  index_image = 0
  
  color = main_plots.make_color_dict()
  style = main_plots.make_style_dict(OPHS) 

  OSKS, JMPS = path.split_at_jumps(ph)
  CLRS = hacks.trace_colors(len(JMPS) + 1, None)

  OMVS = ph.OMVS
  rwd = style['rwd_fill']
  wd_axis = style['wd_axes']

  for imv in range(len(OMVS) + 1):
    omv = OMVS[0:imv]

    fname = "%s/%s_%03d" % (outfolder, tag, index_image)
    dp = (0, 0)
    
    Bfig, cuxlo, cuxhi = main_plots.get_figure_bbox(OPHS, style)
    autoscale = True
    c = main_plots.make_figure_canvas(Bfig, autoscale, style, color)
    
    plot_OPHS(c, OPHS, rwd, dp)

    if len(omv) > 0:
      plot_lines(c, omv, rwd, dp, CLRS, style, color)
    
    hacks.write_plot(c, fname, 3)
    index_image = index_image + 1

  return

# --------------------------------------------------------------- #

def do_test_build(OPHS, tag, Delta, maxband, outfile):

  sys.stderr.write("--- testing {build} tag = %s ---\n" % tag)
   
  # Set the contact time limits:
  CTS = set()
  for oph in OPHS:
    for isd in 0,1:
      for ct in path.get_contacts(oph,isd):
        CTS.add(ct)
        contact.set_tcool_limit(ct, Delta)
  CTS = list(CTS)
  
  wd_jump = 0.00; mp_jump = move_parms.make(wd_jump, acc, csp_jump,  udp)

  # Compute input raster fabtime:
  ydir = (0,1)  # Y direction.
  ytol = 1.0e-5 # Tolerance for Y coordinate variation.
  minlen = 0.1  # Minimum raster length.
  NrastA, TrastA, Nlink, Tlink, Njump, Tjump = raster.analyze_fabtime(OPHS, ydir, ytol, minlen)

  start_time = time.time()
  fph, z, i, j, BPHS, BTCV, TUK = hotfill.solve(OPHS, mp_jump, maxband, False)
  end_time = time.time()

  execution_time = end_time - start_time

  if fph == None:
    fname = "tests/out/" + outfile + ".txt"
    sys.stderr.write("writing %s ...\n" % fname)
    wr = open(fname, "w")
    wr.write("CpuTime:  %7.5f\n" % (execution_time))
    wr.close()
    sys.stderr.write("!! {hotfill.solve} returned {None}\n")
  
  else:
    # Compute input raster fabtime:
    #NrastB, TrastB, Nlink, Tlink, Njump, Tjump = raster.analyze_fabtime([fph,], ydir, ytol, minlen)

    # ??? Should validate the fullpath ???
    #ncold = 1 # Number of coldest contacts to show.
    #CTScold = contact.coldest(fph, CTS, ncold)
    #contact.show_times(sys.stderr, [fph,], CTScold)
    #contact.write_times(outfile, [fph,], CTScold, execution_time, NrastA, TrastA, Nlink, Tlink, Njump, Tjump)

    plot_fp(OPHS, fph, tag)
    #main_plots.plot_plots(fph, OPHS, BPHS, CTScold, tag, 'hotfill', outfolder)
  return

# --------------------------------------------------------------- #

def test_paper_example_B(Delta, maxband, quick):
  tag = "mb%03d_d%05d" % (maxband, int(Delta*100))

  wd_cont = 0.75; mp_cont = move_parms.make(wd_cont, acc, csp_trace, 0.0)
  wd_fill = 1.00; mp_fill = move_parms.make(wd_fill, acc, csp_trace, 0.0)
  wd_link = 0.50; mp_link = move_parms.make(wd_link, acc, csp_trace, 0.0)
  wd_jump = 0.00; mp_jump = move_parms.make(wd_jump, acc, csp_jump,  udp)

  wd_axes = 0.15*min(wd_fill,wd_cont)
  OCRS, OPHS, OLKS, CTS, VGS, EGS = paper_example_B.make_turkey(mp_cont, mp_fill, mp_link, mp_jump)
  
  return OPHS, tag
  # ----------------------------------------------------------------------
  
# --------------------------------------------------------------- #

Delta = 2.3
maxband = 5
outfile  = 'out'
outfolder = './tests/out/'
OPHS, tag = test_paper_example_B(Delta, maxband, False)
do_test_build(OPHS, tag, Delta, maxband, outfile)