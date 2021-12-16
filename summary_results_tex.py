import os
import shutil
import pdfkit
import pandas as pd

# --------------------------------------------------------------- #

def write_begin_table(wr, hbx, hby):
  wr.write(r'\newcommand{\hbx}[1]{\hbox to' + str(hbx) +'mm{\hss #1}}')
  wr.write('\n')
  wr.write(r'\newcommand{\hby}[1]{\hbox to' + str(hby) + 'mm{\hss #1}}')
  wr.write('\n')
  wr.write(r'\newcommand{\mc}[2]{\multicolumn{#1}{c|}{#2}}')
  wr.write('\n')
  wr.write(r'\newcommand{\mcc}[2]{\multicolumn{#1}{c||}{#2}}')
  wr.write('\n')
  wr.write(r'\newcommand{\phz}{\phantom{00}}')
  wr.write('\n')
  wr.write('\n')
  
  wr.write(r'\bigskip')
  wr.write('\n')
  wr.write(r'\footnotesize\rm')
  wr.write('\n')
  wr.write(r'\setlength{\tabcolsep}{2pt}')
  wr.write('\n')
  wr.write('\n')

  return

# --------------------------------------------------------------- #

def write_end_table(wr):
  wr.write('\n')
  wr.write(r' ')
  wr.write(r'\hline')
  wr.write('\n')
  wr.write(r'\end{tabular}')
  wr.write('\n')
  return

# --------------------------------------------------------------- #
# --------------------------------------------------------------- #
# --------------------------------------------------------------- #

def write_header_exp_tcool(wr, delta_keys):  
  line = len(delta_keys[4:])*'r|'
  wr.write(r'\begin{tabular}{|rl||r|r|r|r||' + line + r'}')
  wr.write('\n')

  wr.write(r' ')
  wr.write(r'\hline')
  wr.write('\n')

  wr.write(r'  ')
  wr.write(r'& & & & & &')
  wr.write(r'  ')
  wr.write(r'\mc{' + str(len(delta_keys[4:])) + r'}{{\Hotfill} with various values of $\Delta$} \\')
  wr.write('\n')

  wr.write(r'  ')
  wr.write(r'\cline{7-' + str(6 + len(delta_keys[4:])) + r'}')
  wr.write('\n')

  wr.write(r'  ')
  wr.write(r' & ')
  wr.write(r'Dataset')
  wr.write(r' & ')
  wr.write(r'\hbx{Slic3r}')
  wr.write(r' & ')
  wr.write(r'\hbx{RP3}')
  wr.write(r' & ')
  wr.write(r'\hbx{SCN}')
  wr.write(r' & ')
  wr.write(r'\hbx{SCA}')
  wr.write('\n')
  wr.write(r'  ')

  index_line = 0
    
  for delta in delta_keys[4:]:
    wr.write(r' & ')
    wr.write(r'\hby{' + str(round(float(delta), 1)) + r'}')
    index_line = index_line + 1

    if index_line == 4 and delta != delta_keys[-1]:
      index_line = 0
      wr.write('\n')
      wr.write(r'  ')

  wr.write(r' \\')
  wr.write('\n')
  wr.write(r' ')
  wr.write(r'\hline')
  wr.write('\n')
  wr.write(r' ')
  wr.write(r'\hline')

  wr.write('\n')
  wr.write('\n')
  
  return

# --------------------------------------------------------------- #

def write_exp_tcool(wr, delta_keys, max_band, result_data, df_dict):
  delta_keys  = delta_keys[:3] + delta_keys[4:]
  data_slic3r = result_data['SLIC3R']
  data_rp3    = result_data['RP3']
  data_scn    = result_data['SCANLINE']
  data_sca    = result_data['ALTSCANLINE']
  data_result = result_data[max_band]

  for index in range(1, len(df_dict.keys())):
    key_part = df_dict[index]['part']
    key_angle = df_dict[index]['angle']

    wr.write(r'  ')
    wr.write(str(index))
    wr.write(r' & ')

    wr.write(r'\sfo{' + key_part + r'}{' + key_angle + r'}')
    wr.write(r' & ')

    wr.write(str(data_slic3r[key_part][key_angle]['SLIC3R']['CoolTime']))
    wr.write(r' & ')

    wr.write(str(data_rp3[key_part][key_angle]['RP3']['CoolTime']))
    wr.write(r' & ')

    wr.write(str(data_scn[key_part][key_angle]['SCANLINE']['CoolTime']))
    wr.write(r' & ')

    wr.write(str(data_sca[key_part][key_angle]['ALTSCANLINE']['CoolTime']))
    
    for key_delta in delta_keys[3:]:
      wr.write(r' & ')

      if data_result[key_part][key_angle][key_delta]['CoolTime'] != '-':
        wr.write(str(data_result[key_part][key_angle][key_delta]['CoolTime']))
      else:
        wr.write(r'---')
    
    wr.write(r' \\')
    wr.write('\n')

  return

# --------------------------------------------------------------- #

def create_exp_tcool(folder_out, delta_keys, max_band, result_data, df_dict):
  filename = folder_out + 'exp-tcool.tex'
  wr = open(filename, 'w')
  write_begin_table(wr, 8, 6.5)
  write_header_exp_tcool(wr, delta_keys)
  write_exp_tcool(wr, delta_keys, max_band, result_data, df_dict)
  write_end_table(wr)
  return

# --------------------------------------------------------------- #
# --------------------------------------------------------------- #
# --------------------------------------------------------------- #

def write_header_exp_tfab(wr, delta_keys):
  line = 'r|'*len(delta_keys[4:])
  wr.write(r'\begin{tabular}{|rl||r|r|r|r|r||' + line + r'}')
  wr.write('\n')

  wr.write(r' ')
  wr.write(r'\hline')
  wr.write('\n')
  
  wr.write(r'  ')
  wr.write(r'& & & & & & &')
  wr.write('\n')
  wr.write(r'  ')
  wr.write(r'\mc{' + str(len(delta_keys[4:])) + r'}{{\Hotfill} for various values of $\Delta$} \\')
  wr.write('\n')
  wr.write(r'  ')
  wr.write(r'\cline{8-' + str(7 + len(delta_keys[4:])) + r'}')
  wr.write('\n')
  
  wr.write(r'  ')
  wr.write(r' & ')
  wr.write(r'Dataset')
  wr.write(r' & ')
  wr.write(r'$\Trast$')
  wr.write(r' & ')
  wr.write(r'Slic3r')
  wr.write(r' & ')
  wr.write(r'RP3')
  wr.write(r' & ')
  wr.write(r'SCN')
  wr.write(r' & ')
  wr.write(r'SCA')
  wr.write('\n')
  wr.write(r'  ')

  for delta in delta_keys[4:]:
    wr.write(r' & ')
    wr.write(r'\mc1{' + str(delta) + r'}')
  
  wr.write(r' \\')
  wr.write('\n')
  wr.write(r' ')
  wr.write(r'\hline')
  wr.write('\n')
  wr.write(r' ')
  wr.write(r'\hline')

  wr.write('\n')
  wr.write('\n')

  return

# --------------------------------------------------------------- #

def write_exp_tfab(wr, max_band, delta_keys, result_data, df_dict):
  data_slic3r = result_data['SLIC3R']
  data_rp3    = result_data['RP3']
  data_scn   = result_data['SCANLINE']
  data_sca  = result_data['ALTSCANLINE']
  data_result = result_data[max_band]

  for index in range(1, len(df_dict.keys())):
    key_part = df_dict[index]['part']
    key_angle = df_dict[index]['angle']
    tras = data_scn[key_part][key_angle]['SCANLINE']['RastTime']

    write_exp_tfab_line1(wr, index, key_part, key_angle, delta_keys, data_slic3r, data_rp3, data_scn, data_sca, data_result, tras)
    write_exp_tfab_line2(wr, key_part, key_angle, delta_keys, data_slic3r, data_rp3, data_scn, data_sca, data_result)
    
    wr.write('\n')
  
  return

# --------------------------------------------------------------- #

def write_exp_tfab_line1(wr, index, key_part, key_angle, delta_keys, data_slic3r, data_rp3, data_scn, data_sca, data_result, tras):
  wr.write(r'  ')
  wr.write(r'\multirow{2}{*}{' + str(index) + r'}')
  wr.write(r' & ')
  wr.write('\n')
  
  wr.write(r'  ')
  wr.write(r'\multirow{2}{*}{\sfo{' + key_part + r'}{' + key_angle + r'}}')
  wr.write(r' & ')
  wr.write('\n')

  wr.write(r'  ')
  wr.write(tras)
  wr.write(r' & ')
  wr.write('\n')
  
  wr.write(r'  ')
  wr.write(data_slic3r[key_part][key_angle]['SLIC3R']['FabTime'])
  wr.write(r' & ')
  wr.write('\n')

  wr.write(r'  ')
  wr.write(data_rp3[key_part][key_angle]['RP3']['FabTime'])
  wr.write(r' & ')
  wr.write('\n')

  wr.write(r'  ')
  wr.write(data_scn[key_part][key_angle]['SCANLINE']['FabTime'])
  wr.write(r' & ')
  wr.write('\n')

  wr.write(r'  ')
  wr.write(data_sca[key_part][key_angle]['ALTSCANLINE']['FabTime'])

  for key_delta in delta_keys[4:]:
    wr.write(r' & ')
    wr.write('\n')
    wr.write(r'  ')

    if data_result[key_part][key_angle][key_delta]['FabTime'] != '-':          
      tfab = float(data_result[key_part][key_angle][key_delta]['FabTime'])
      wr.write(str(tfab))
    
    else:
      wr.write(r'\multirow{2}{*}{---}')
  
  wr.write(r' \\')
  wr.write('\n')
  
  return

# --------------------------------------------------------------- #

def write_exp_tfab_line_per(wr, tfab_ref, tfab):
  dif = (tfab-tfab_ref)/tfab_ref
  dif = 100*dif

  wr.write(r'$\hbx{')
  wr.write('%+.0f' % dif)
  wr.write(r'\%}$')

  return

# --------------------------------------------------------------- #

def write_exp_tfab_line2(wr, key_part, key_angle, delta_keys, data_slic3r, data_rp3, data_scn, data_sca, data_result):
  tfab_ref = float(min(data_slic3r[key_part][key_angle]['SLIC3R']['FabTime'], data_rp3[key_part][key_angle]['RP3']['FabTime']))

  wr.write(r'  ')
  wr.write(r'& & & & ')

  wr.write('\n')
  wr.write(r'  ')
  wr.write(r'& ')
  tfab = float(data_scn[key_part][key_angle]['SCANLINE']['FabTime'])
  write_exp_tfab_line_per(wr, tfab_ref, tfab)

  wr.write('\n')
  wr.write(r'  ')
  wr.write(r'& ')
  tfab = float(data_sca[key_part][key_angle]['ALTSCANLINE']['FabTime'])
  write_exp_tfab_line_per(wr, tfab_ref, tfab)
  
  for key_delta in delta_keys[4:]:
    wr.write('\n')
    wr.write(r'  ')
    wr.write(r'& ')

    if data_result[key_part][key_angle][key_delta]['FabTime'] != '-': 
      tfab = float(data_result[key_part][key_angle][key_delta]['FabTime'])
      write_exp_tfab_line_per(wr, tfab_ref, tfab)
    
    else:
      wr.write(r' ')
  
  wr.write(r' \\ ')
  wr.write(r'\hline')
  wr.write('\n')
  
  return

# --------------------------------------------------------------- #

def create_exp_tfab(folder_out, delta_keys, max_band, result_data, df_dict):
  filename = folder_out + 'exp-tfab.tex'
  wr = open(filename, 'w')
  write_begin_table(wr, 9, 6.5)
  write_header_exp_tfab(wr, delta_keys)
  write_exp_tfab(wr, max_band, delta_keys, result_data, df_dict)
  write_end_table(wr)
  return

# --------------------------------------------------------------- #
# --------------------------------------------------------------- #
# --------------------------------------------------------------- #

def write_header_exp_tcpu(wr, max_bands):
  line = 'r|'*len(max_bands)
  wr.write(r'  \begin{tabular}{|rl||r|r||' + line + r'}')
  wr.write('\n')

  wr.write(r' ')
  wr.write(r'\hline')
  wr.write('\n')

  wr.write(r'  ')
  wr.write(r'& & & &')
  wr.write(r'\mc{' + str(len(max_bands)) + r'}{$\maxband$}')
  wr.write(r' \\')
  wr.write('\n')

  wr.write(r'  ')
  wr.write(r'\cline{5-' + str(4 + len(max_bands)) + r'}')
  wr.write('\n')

  wr.write(r'  ')
  wr.write(r' & ')
  wr.write(r'Dataset')
  wr.write(r' & ')
  wr.write(r'\hbx{$n$}')
  wr.write(r' & ')
  wr.write(r'\hbx{$m$}')

  for max_band in max_bands:
    wr.write('\n')
    wr.write(r' ')
    wr.write(r' & ')
    wr.write(r'\hbx{' + max_band + '}')
  
  wr.write(r' \\')
  wr.write('\n')
  wr.write(r' ')
  wr.write(r'\hline')
  wr.write('\n')
  wr.write(r' ')
  wr.write(r'\hline')
  
  wr.write('\n')
  wr.write('\n')

  return

# --------------------------------------------------------------- #

def write_exp_tcpu(wr, max_bands, delta, result_data, df_dict, dict_scanline):
  for index in range(1, len(df_dict.keys())):
    key_part = df_dict[index]['part']
    key_angle = df_dict[index]['angle']
    nrast = result_data[max_bands[0]][key_part][key_angle][delta]['Nrast']

    wr.write(r'  ')
    wr.write(str(index))
    wr.write(r' & ')

    wr.write(r'\sfo{' + key_part + r'}{' + key_angle + r'}')
    wr.write(r' & ')

    wr.write(nrast.replace(' ', ''))
    wr.write(r' & ')
    wr.write('\n')

    wr.write(r'  ')
    wr.write(str(dict_scanline[key_part][key_angle]['nsc']))

    for key_band in max_bands:
      wr.write(r' & ')
      wr.write('\n')
      wr.write(r'  ')
      wr.write(result_data[key_band][key_part][key_angle][delta]['CpuTime'])
  
    wr.write(r' \\')
    wr.write('\n')
  
  wr.write('\n')
        
  return

# --------------------------------------------------------------- #

def create_exp_tcpu(folder_out, delta, max_bands, result_data, df_dict, dict_scanline):
  filename = folder_out + 'exp-tcpu.tex'
  wr = open(filename, 'w')
  write_begin_table(wr, 8, 6.5)
  write_header_exp_tcpu(wr, max_bands)
  write_exp_tcpu(wr, max_bands, delta, result_data, df_dict, dict_scanline)
  write_end_table(wr)
  return

# --------------------------------------------------------------- #
# --------------------------------------------------------------- #
# --------------------------------------------------------------- #