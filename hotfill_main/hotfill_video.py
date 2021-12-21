import hotfill
import paper_example_B
import path 
import move_parms
import contact
import os
import sys
import main_plots
import hacks
import move
import pyx
import matplotlib.pyplot as plt

# --------------------------------------------------------------- #

# Some arbitrary dynamics parameters:
acc = 3000      # Acceleration/deceleration for moves and jumps (mm/s^2).
csp_trace = 40  # Cruise speed for traces (mm/s).
csp_jump = 130  # Cruise speed for jumps (mm/s).
udp = 0.05      # Trace/jump transition penalty (s).

# --------------------------------------------------------------- #

def plot_OPHS(c, dp, rwd, OPHS):
  clr = pyx.color.rgb(0.850, 0.850, 0.850)
  
  for oph in OPHS:
    omv = path.elem(oph, 0)
    move.plot_layer(c, omv, dp, clr, rwd, dashed=False, wd_dots=0, sz_arrow=None)

  return

# --------------------------------------------------------------- #

def get_open_contacts(CTS, OMVS):
  cts = []
  
  for ct in CTS:
    mvs = ct.side
    mv0, d = move.unpack(mvs[0])
    mv1, d = move.unpack(mvs[1])

    if mv0 in OMVS and mv1 not in OMVS:
      cts.append(ct)
    elif mv0 not in OMVS and mv1 in OMVS:
      cts.append(ct)
  
  return cts

# --------------------------------------------------------------- #

def plot_contacts(c, dp, style, color, CTS, OMVS):
  wd_ctac = style['wd_ctac']
  dashpat = (1.0*wd_ctac, 2.0*wd_ctac)
  sz_tics = 0
  arrows_ct = False
  ext_ct = 0
  
  cts = get_open_contacts(CTS, OMVS)

  for ct in cts:
    contact.plot_single(c, ct, dp, color['ctac'], dashpat, ext_ct, wd=wd_ctac, sz_tic=sz_tics, arrow=arrows_ct)

  return

# --------------------------------------------------------------- #

def plot_contacts_time(c, style, CTS, OMVS, ph):
  wd_ctac = style['wd_ctac']
  cts = get_open_contacts(CTS, OMVS)
  radius = 0.35

  clr_tot  = pyx.color.rgb(0.000, 0.000, 1.000)
  clr_cool = pyx.color.rgb(1.000, 0.250, 0.000)

  for ct in cts:
    x = (ct.pts[0][0] + ct.pts[1][0])/2
    y = ct.pts[0][1] + (wd_ctac + 1.1*radius)   

    circle_ph = pyx.path.circle(x, y, radius)
    c.stroke(circle_ph, [pyx.style.linewidth.Thick, clr_cool, pyx.deco.filled([clr_cool])])

    t_limit = contact.tcool_limit(ct)
    t_fin = ph.cumtex[len(OMVS)]
    t_cool0, t_cool1 = contact.path_tcovs(ph, ct)

    t_cool_min = min(t_cool0, t_cool1)

    t_cool = t_fin-t_cool_min
    t_cool = t_cool/t_limit

    sector_ph = pyx.path.path(pyx.path.moveto(x, y), pyx.path.arc(x, y, radius, 90, 90-(360*t_cool)), pyx.path.closepath())
    c.stroke(sector_ph, [pyx.style.linewidth.Thick, clr_tot, pyx.deco.filled([clr_tot])])

  return

# --------------------------------------------------------------- #

def plot_lines(c, dp, rwd, style, color, CLRS, OMVS):
  rwd = style['rwd_fill']
  index_clr = 0
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

def plot_fp(OPHS, ph, CTS, tag, outfolder):
  index_image = 0
  
  color = main_plots.make_color_dict()
  style = main_plots.make_style_dict(OPHS) 

  OSKS, JMPS = path.split_at_jumps(ph)
  CLRS = hacks.trace_colors(len(JMPS) + 1, None)

  OMVS = ph.OMVS
  OMVS_ph = [ move.unpack(omv)[0] for omv in OMVS]

  rwd = style['rwd_fill']
  Bfig, cuxlo, cuxhi = main_plots.get_figure_bbox(OPHS, style)
  autoscale = True
  dp = (0, 0)

  for imv in range(len(OMVS) + 1):
    omv = OMVS[0:imv]
    omv_ph = OMVS_ph[0:imv]

    fname = "%s/%s_%03d" % (outfolder, tag, index_image)
  
    c = main_plots.make_figure_canvas(Bfig, autoscale, style, color)
    
    plot_OPHS(c, dp, rwd, OPHS)

    if len(omv) > 0:
      plot_contacts(c, dp, style, color, CTS, omv_ph)
      plot_lines(c, dp, rwd, style, color, CLRS, omv)
      plot_contacts_time(c, style, CTS, omv_ph, ph)
    
    hacks.write_plot(c, fname, 3)
    index_image = index_image + 1

  #plot_heat_map(OMVS, CTS, ph, tag, outfolder, Bfig, autoscale, style, color, dp, rwd, CLRS)

  return

# --------------------------------------------------------------- #

def do_test_build(OPHS, tag, Delta, maxband, outfolder):

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

  fph, z, i, j, BPHS, BTCV, TUK = hotfill.solve(OPHS, mp_jump, maxband, False)

  if fph != None:
    plot_fp(OPHS, fph, CTS, tag, outfolder)
    alg = 'hotfill'

    CTScold = contact.coldest(fph, CTS, 1)

    main_plots.plot_plots(fph, OPHS, BPHS, CTScold, tag, alg, outfolder)
    os.remove(outfolder + "/" + alg + "_" + tag + "_printer_parms.tex")
    os.remove(outfolder + "/" + alg + "_" + tag + "_times.tex")
  
  return

# --------------------------------------------------------------- #

def test_paper_example_B(Delta, maxband):
  tag = "mb%03d_d%05d" % (maxband, int(Delta*100))

  wd_cont = 0.75; mp_cont = move_parms.make(wd_cont, acc, csp_trace, 0.0)
  wd_fill = 1.00; mp_fill = move_parms.make(wd_fill, acc, csp_trace, 0.0)
  wd_link = 0.50; mp_link = move_parms.make(wd_link, acc, csp_trace, 0.0)
  wd_jump = 0.00; mp_jump = move_parms.make(wd_jump, acc, csp_jump,  udp)

  OCRS, OPHS, OLKS, CTS, VGS, EGS = paper_example_B.make_turkey(mp_cont, mp_fill, mp_link, mp_jump)
  
  return OPHS, tag
  # ----------------------------------------------------------------------
  
# --------------------------------------------------------------- #

def main(Delta, maxband):
  outfolder = './tests/out/mb%03d_d%05d' % (maxband, Delta*100)
  if not os.path.exists(outfolder):
    os.makedirs(outfolder)

  OPHS, tag = test_paper_example_B(Delta, maxband)
  do_test_build(OPHS, tag, Delta, maxband, outfolder)

  return

# --------------------------------------------------------------- #

#main(Delta=1.7, maxband=5)
#main(Delta=2.3, maxband=5)
#main(Delta=3.5, maxband=5)
main(Delta=50, maxband=50)
