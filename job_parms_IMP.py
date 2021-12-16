# Implementaton of module {job_parms}.
# Last edited on 2021-11-09 01:55:33 by stolfi

import job_parms
from math import sqrt, sin, cos, floor, ceil, inf, nan, pi

# ??? Parameters subject to changes ???

def units():
  units = {}
  # PRINTER:
  units['acceleration']           = "mm/s^2"
  units['max_extrusion_speed']    = "mm/s"
  units['max_jump_speed']         = "mm/s"
  units['max_z_speed']            = "mm/s"
  units['max_filament_speed']     = "mm/s"
  units['extrusion_on_off_time']  = "s" 
  # JOB:
  units['job_contour_speed']      = "mm/s"
  units['job_filling_speed']      = "mm/s"
  units['job_jump_speed']         = "mm/s"
  units['job_z_speed']            = "mm/s"
  units['slice_thickness']        = "mm"
  units['solid_raster_width']     = "mm"
  units['contour_trace_width']    = "mm"
  units['nozzle_temperature']     = "C"
  units['filament_diameter']      = "mm"

  return units


def typical_js():
  parms = {}

  # PRINTER PARAMETERS

  # Timing parameters. These are limits for the job-specific values below.
  parms['acceleration']           = 3000  # Acceleration/deceleration at ends of moves (mm/s^2).
  parms['max_extrusion_speed']    = 20    # Max cruise speed when extruding (mm/s).
  parms['max_jump_speed']         = 40    # Max cruise speed when jumping (mm/s).
  parms['max_z_speed']            = 40    # Max cruise speed when moving vertically (mm/s).
  parms['max_filament_speed']     = 2     # Max speed of filament feeding (mm/s).
  parms['extrusion_on_off_time']  = 0.500 # Nozzle/filament up/down time (s). 
  
  # JOB PARAMETERS

  # Nominal geometry design parameters.
  # These are nominal dimensions used in computing the 
  # geometry of paths. They are also used to compute the
  # volume of material to extrude during traces.
  parms['contour_trace_width']    = 0.25 # Spacing between adjacent traces (mm).
  parms['solid_raster_width']     = 0.50 # Spacing between adjacent traces (mm).
  parms['slice_thickness']        = 0.25 # Slice thickness (mm).

  # Timing parameters to use in a specific job:
  parms['job_contour_speed']      = 10 # Cruise speed to use when extruding contours (mm/s).
  parms['job_filling_speed']      = 20 # Cruise speed to use when extruding fillings (mm/s).
  parms['job_jump_speed']         = 40 # Cruise speed to use when jumping (mm/s).
  parms['job_z_speed']            = 40 # Cruise speed to use when moving vertically (mm/s).

  # Parameters for G-code generation:
  parms['nozzle_temperature']     = 210  # Nozzle temperature (C).
  parms['filament_diameter']      = 1.75 # Diameter of feedstock filament (mm).

  return parms
  # ----------------------------------------------------------------------

def typical_elis():
  parms = {}

  # PRINTER PARAMETERS

  # Timing parameters. These are limits for the job-specific values below.
  parms['acceleration']           = 3000  # Acceleration/deceleration at ends of moves (mm/s^2).
  parms['max_extrusion_speed']    = 40    # Max cruise speed when extruding (mm/s).
  parms['max_jump_speed']         = 100   # Max cruise speed when jumping (mm/s).
  parms['max_z_speed']            = 40    # Max cruise speed when moving vertically (mm/s).
  parms['max_filament_speed']     = 2     # Max speed of filament feeding (mm/s).
  parms['extrusion_on_off_time']  = 0.050 # Nozzle/filament up/down time (s). 
  
  # JOB PARAMETERS

  # Nominal geometry design parameters.
  # These are nominal dimensions used in computing the 
  # geometry of paths. They are also used to compute the
  # volume of material to extrude during traces.
  parms['contour_trace_width']    = 0.40 # Spacing between adjacent traces (mm).
  parms['solid_raster_width']     = 0.40 # Spacing between adjacent traces (mm).
  parms['slice_thickness']        = 0.25 # Slice thickness (mm).

  # Timing parameters to use in a specific job:
  parms['job_contour_speed']      = 40  # Cruise speed to use when extruding contours (mm/s).
  parms['job_filling_speed']      = 40  # Cruise speed to use when extruding fillings (mm/s).
  parms['job_jump_speed']         = 100 # Cruise speed to use when jumping (mm/s).
  parms['job_z_speed']            = 40  # Cruise speed to use when moving vertically (mm/s).

  # Parameters for G-code generation:
  parms['nozzle_temperature']     = 210  # Nozzle temperature (C).
  parms['filament_diameter']      = 1.75 # Diameter of feedstock filament (mm).

  return parms
  # ----------------------------------------------------------------------

def typical_paper():
  parms = {}

  # PRINTER PARAMETERS

  # Timing parameters. These are limits for the job-specific values below.
  parms['acceleration']           = 3000  # Acceleration/deceleration at ends of moves (mm/s^2).
  parms['max_extrusion_speed']    = 40    # Max cruise speed when extruding (mm/s).
  parms['max_jump_speed']         = 130   # Max cruise speed when jumping (mm/s).
  parms['max_z_speed']            = 40    # Max cruise speed when moving vertically (mm/s).
  parms['max_filament_speed']     = 2     # Max speed of filament feeding (mm/s).
  parms['extrusion_on_off_time']  = 0.050 # Nozzle/filament up/down time (s). 
  
  # JOB PARAMETERS

  # Nominal geometry design parameters.
  # These are nominal dimensions used in computing the 
  # geometry of paths. They are also used to compute the
  # volume of material to extrude during traces.
  parms['contour_trace_width']    = 0.40 # Spacing between adjacent traces (mm).
  parms['solid_raster_width']     = 0.40 # Spacing between adjacent traces (mm).
  parms['slice_thickness']        = 0.25 # Slice thickness (mm).

  # Timing parameters to use in a specific job:
  parms['job_contour_speed']      = 40  # Cruise speed to use when extruding contours (mm/s).
  parms['job_filling_speed']      = 40  # Cruise speed to use when extruding fillings (mm/s).
  parms['job_jump_speed']         = 100 # Cruise speed to use when jumping (mm/s).
  parms['job_z_speed']            = 40  # Cruise speed to use when moving vertically (mm/s).

  # Parameters for G-code generation:
  parms['nozzle_temperature']     = 210  # Nozzle temperature (C).
  parms['filament_diameter']      = 1.75 # Diameter of feedstock filament (mm).

  return parms
  # ----------------------------------------------------------------------

def slow():
  parms = typical_js()

  # PRINTER:
  parms['acceleration']           = 30    # Acceleration/deceleration at ends (mm/s^2).
  parms['extrusion_on_off_time']  = 0.200 # Nozzle/filament up/down time (s).

  # JOB:
  parms['contour_trace_width']    = 0.25  # Spacing between adjacent traces (mm).
  parms['solid_raster_width']     = 0.50  # Spacing between adjacent traces (mm).
  parms['job_contour_speed']      = 5     # Cruise speedto use when extruding contours (mm/s).
  parms['job_filling_speed']      = 10    # Cruise speed to use when extruding fillings (mm/s).
  parms['job_jump_speed']         = 20    # Cruise speed to use when jumping (mm/s).
  parms['job_z_speed']            = 5     # Cruise speed to use when moving vertically (mm/s).
  return parms
  # ----------------------------------------------------------------------

def very_slow():
  parms = typical_js()

  # PRINTER:
  parms['acceleration']           = +inf  # Acceleration/deceleration at ends (mm/s^2).
  parms['extrusion_on_off_time']  = 1     # Nozzle/filament up/down time (s).

  # JOB:
  parms['job_contour_speed']      = 1     # Cruise speedto use when extruding contours (mm/s).
  parms['job_filling_speed']      = 1     # Cruise speed to use when extruding fillings (mm/s).
  parms['job_jump_speed']         = 2     # Cruise speed to use when jumping (mm/s).
  parms['job_z_speed']            = 1     # Cruise speed to use when moving vertically (mm/s).
  return parms
  # ----------------------------------------------------------------------

def write(wr, pref, parms, suff):
  uns = units()
  for key, val in parms.items():
    if pref != None:
      wr.write(pref)
    wr.write("%-24s" % (key + ":"))
    if type(val) is float:
      wr.write("%8.3f" % val)
    elif type(val) is int:
      wr.write("%8d" % val)
    elif type(val) is bool:
      wr.write("%8s" % (("F", "T")[int(val)]))
    elif type(val) is str:
      wr.write("%8s" % val)
    else:
      wr.write("%8s" % str(val))
    un = uns[key]
    if un != None:
      wr.write(" %s" % un)
    if suff != None: 
      w.write(suff)
    wr.write("\n")
  wr.flush()
  return
  # ----------------------------------------------------------------------
