import path
from math import sqrt, sin, cos, atan2, log, exp, floor, ceil, inf, nan, pi
import move_parms
import rn
import contact
import move
import sys


# --------------------------------------------------------------- #

def read(rd, islice, angle, shift, mp_cont, mp_fill, mp_link, mp_jump):
  if 'rp3' in rd.name:
    C, PTS_raster, Z = read_gcode_rp3(rd, islice, mp_fill, angle, shift)
    tol = 0.005
  elif 'slic3r' in rd.name:
    C, PTS_raster, Z = read_gcode_slic3r(rd, islice, mp_fill, angle, shift)
    tol = 0.05
  else:
    return None, None

  TRS_raster = create_lines(PTS_raster, mp_fill, mp_link, tol)
  find_contacts(TRS_raster, mp_fill, tol)

  OPHS = []

  for TRS in TRS_raster:
    oph = path.concat(TRS, False, mp_jump)
    OPHS.append(oph)      
  
  return OPHS, Z

# --------------------------------------------------------------- #

def read_gcode_rp3(rd, islice, mp_fill, angle, shift):
  C = list()
  R = list()
  Z = 0

  PTS_read = list()

  read_slice = False
  is_raster = False
  new_list = False

  xdir = (cos(angle), sin(angle))   # Dir vector parallel to rasters, "left" to "right".
  ydir = (-xdir[1], xdir[0])        # Dir vector perpendicular to rasters, "up".

  def unrotate_and_shift(p_raw):
    # Given a point as read from the file,
    # separated by the character {sep}, returns the point whith those coordinates (as pair of {float}s) 
    # after rotating by {-alpha} and shifting by {shift}
    #
    xr = rn.dot(p_raw, xdir) + shift[0]
    yr = rn.dot(p_raw, ydir) + shift[1]
    p = (xr,yr)
    return p
    # ....................................................................

  for line in rd:
    if '(Layer' in line:
      layer = line.replace('(Layer', '')
      layer = int(layer.replace(')', ''))

      if layer == islice:
        read_slice = True
      elif layer > islice:
        break

    elif read_slice:
      if 'Z' in line:
        line = line.split(' ')
        for iz in line:
          if 'Z' in iz:
            Z = iz.replace('Z', '')

      elif '(Offset ' in line:
        is_raster = False
        
        if 'Offset 1)' not in line: 
          new_list = True 
      
      elif '(Raster ' in line:
        if is_raster == False:
          if len(PTS_read) > 0:
            C.append(PTS_read)
          PTS_read = list()
          is_raster = True

        if 'Raster 1)' not in line: 
          new_list = True 

      elif line[0:3] == 'G1 ' and 'X' in line and 'Y' in line:
        if new_list == True:
          if not is_raster:
            if len(PTS_read) > 0:
              C.append(PTS_read)
          else: 
            if len(PTS_read) > 0:
              R.append(PTS_read)

          PTS_read = list()
          new_list = False

        XY = line.replace('  ', ' ').split(' ')
              
        for i in XY:
          if 'X' in i:
            x = float(i.replace('X',''))
          elif 'Y' in i:
            y = float(i.replace('Y',''))
        
        p_raw = (x, y) # One endpoint.
        p = unrotate_and_shift(p_raw)

        PTS_read.append(p)

  if len(PTS_read) > 0:
    R.append(PTS_read)
  
  return C, R, Z

# --------------------------------------------------------------- #

def read_gcode_slic3r(rd, islice, mp_fill, angle, shift):
  C = list()
  R = list()
  Z = islice*0.250
  
  PTS_read = list()
  new_list = False

  read_slice = False
  layer_z = 'Z' + str(Z)

  xdir = (cos(angle), sin(angle))   # Dir vector parallel to rasters, "left" to "right".
  ydir = (-xdir[1], xdir[0])        # Dir vector perpendicular to rasters, "up".

  def unrotate_and_shift(p_raw):
    # Given a point as read from the file,
    # separated by the character {sep}, returns the point whith those coordinates (as pair of {float}s) 
    # after rotating by {-alpha} and shifting by {shift}
    #
    xr = rn.dot(p_raw, xdir) + shift[0]
    yr = rn.dot(p_raw, ydir) + shift[1]
    p = (xr,yr)
    return p
    # ....................................................................
  
  for line in rd:
    if layer_z in line and 'lift nozzle' not in line:
      read_slice = True
    
    elif 'Z' in line:
      read_slice = False
      
    elif read_slice:
      if 'X' in line and 'Y' in line and 'F7800' in line:
        new_list = True 

      if line[0:3] == 'G1 ' and 'X' in line and 'Y' in line:
        XY = line.replace('  ', ' ').split(' ')
        
        for i in XY:
          if 'X' in i:
            x = float(i.replace('X',''))
          elif 'Y' in i:
            y = float(i.replace('Y',''))

        p_raw = (x, y) # One endpoint.
        p = unrotate_and_shift(p_raw)
        PTS_read.append(p)

        if new_list == True:
          aux = None
    
          if len(PTS_read) > 0:
            aux = PTS_read[-1]   
            PTS_read.pop()   
      
            if len(PTS_read) > 0:
              if rn.dist(PTS_read[0], PTS_read[-1]) < 0.09:
                C.append(PTS_read)
              else:
                R.append(PTS_read)

          PTS_read = list()
    
          if aux != None:
            PTS_read.append(aux)
      
          new_list = False

  if len(PTS_read) > 0:
    if rn.dist(PTS_read[0], PTS_read[-1]) < 0.09:
      C.append(PTS_read)
    else:
      R.append(PTS_read)

  return C, R, Z

# --------------------------------------------------------------- #

def create_lines(PTS_raster, mp_fill, mp_link, tol):
  TRS_raster = []

  for PTS in PTS_raster:
    TRS = []
    p = None

    for index_point in range(len(PTS)):
      q = PTS[index_point]

      if p != None:
        if (abs(p[1] - q[1]) <= tol):
          mp = mp_fill
          q = (q[0], p[1])
          tag = 'TRACE'
      
        else:
          mp = mp_link
          tag = 'LINK'
        
        mv = move.make(p, q, mp)
        ph = path.from_moves((mv,))
        path.set_name(ph, tag, False)
        TRS.append(ph)

      p = q
    
    TRS_raster.append(TRS)

  return TRS_raster

# --------------------------------------------------------------- #

def find_contacts(TRS, mp_fill, tol):
  CTS = []

  ystep = move_parms.width(mp_fill) # Expected spacing of rasters.

  for itrc in range(len(TRS)):
    traces = TRS[itrc]

    for i in range(len(traces)):
      oph0 = traces[i]

      if 'TRACE' in path.get_name(oph0):
        for j in range(i + 1, len(traces)):
          oph1 = traces[j]
          
          if 'TRACE' in path.get_name(oph1):
            ph0, ph1 = check_contacts(oph0, oph1, ystep, tol)

            if ph0 != None and ph1 != None:
              ct = create_and_attach_contact(ph0, ph1)
              CTS.append(ct)
        
        for jtrc in range(itrc + 1, len(TRS)):
          jtraces = TRS[jtrc]
          
          for j in range(len(jtraces)):
            oph1 = jtraces[j]
            
            if 'TRACE' in path.get_name(oph1):
              ph0, ph1 = check_contacts(oph0, oph1, ystep, tol)

              if ph0 != None and ph1 != None:
                ct = create_and_attach_contact(ph0, ph1)
                CTS.append(ct)         
            
  return CTS

# --------------------------------------------------------------- #

def check_contacts(oph0, oph1, ystep, tol):
  imv0, dr0 = move.unpack(path.elem(oph0, 0))
  imv1, dr1 = move.unpack(path.elem(oph1, 0))

  pts0 = move.endpoints(imv0)
  pts1 = move.endpoints(imv1)

  if (abs(pts0[0][1] - pts1[0][1]) - ystep) <= tol:
    if pts0[0][1] < pts1[0][1]:
      return oph0, oph1
    else:
      return oph1, oph0

  return None, None

# --------------------------------------------------------------- #

def create_and_attach_contact(oph0, oph1):
  mv0, dr0 = move.unpack(path.elem(oph0, 0))
  mv1, dr1 = move.unpack(path.elem(oph1, 0))
  mvdir = (1, 0)
  
  wd0 = move.width(mv0)
  wd1 = move.width(mv1)
  tol = 0.20*min(wd0, wd1)  # Tolerance for overlaps, tilts, etc.
  ct = contact.from_moves(mv0, mv1, mvdir, 0.1, 0.05, tol)

  if ct == None: 
    sys.stderr.write("!! contact creation failed\n")
    move.show(sys.stderr, "  mv0 = ", mv0, "\n", 0)
    move.show(sys.stderr, "  mv1 = ", mv1, "\n", 0)
  else:
    contact.add_side_path(ct, 0, oph0, 0)
    contact.add_side_path(ct, 1, oph1, 0)

    path.add_contact(oph0, 0, ct)
    path.add_contact(oph1, 1, ct)

  return ct

# --------------------------------------------------------------- #

def solve(OPHS, mp_jump):
  fph = path.concat(OPHS, False, mp_jump)
  return fph