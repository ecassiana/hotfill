from PIL.Image import register_decoder
import gcode_read_IMP

def solve(OPHS, mp_jump):
  return gcode_read_IMP.solve(OPHS, mp_jump)

def read(rd, islice, angle, shift, mp_cont, mp_fill, mp_link, mp_jump):
  return gcode_read_IMP.read(rd, islice, angle, shift, mp_cont, mp_fill, mp_link, mp_jump)