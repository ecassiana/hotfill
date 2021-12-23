import os
import rn

# --------------------------------------------------------------- #

folder_in = './tests/in/'
file_out = './tests/results/tab-dataset-summary.tex'
details_file = './tests/details.txt'

# --------------------------------------------------------------- #

def get_info(folder_in, file_in):
  min_l = None
  max_l = None
  tot_l = 0
  n_ras = 0

  with open(folder_in + file_in, 'r') as f:
    for line in f:
      if line[0] == 'R':
        r = line.split(',')
        p = (r[1], r[2]) ; p = (float(p[0]), float(p[1]))
        q = (r[3], r[4]) ; q = (float(q[0]), float(q[1]))
        dist = rn.dist(p, q)
        
        if min_l == None or min_l > dist:
          min_l = dist
        
        if max_l == None or max_l < dist:
          max_l = dist
        
        tot_l = tot_l + dist
        n_ras = n_ras + 1
        
  return min_l, max_l, tot_l, n_ras

# --------------------------------------------------------------- #

def write_header(file_out):
  with open(file_out, 'w') as wr:
    wr.write(r'\footnotesize\rm')
    wr.write('\n')
    wr.write(r'\setlength{\tabcolsep}{2mm}')
    wr.write('\n')
    wr.write('\n')

    wr.write(r'\begin{tabular}{|rl|c|r|r|r|c|c|rrr|}')
    wr.write('\n')
    wr.write(r' ')
    wr.write(r'\hline')
    wr.write('\n')

    wr.write(r'  ')
    wr.write(r'\multirow{2}{*}{}')
    wr.write(r'& ')
    wr.write('\n')
    wr.write(r'  ')
    wr.write(r'\multirow{2}{*}{Dataset}')
    wr.write(r'& ')
    wr.write('\n')
    wr.write(r'  ')
    wr.write(r'\multirow{2}{*}{Slice}')
    wr.write(r'& ')
    wr.write('\n')
    wr.write(r'  ')
    wr.write(r'\multirow{2}{*}{$\theta$}')
    wr.write(r'& ')
    wr.write('\n')
    wr.write(r'  ')
    wr.write(r'\multirow{2}{*}{$X$}')
    wr.write(r'& ')
    wr.write('\n')
    wr.write(r'  ')
    wr.write(r'\multirow{2}{*}{$Y$}')
    wr.write(r'& ')
    wr.write('\n')
    wr.write(r'  ')
    wr.write(r'\multirow{2}{*}{$n$}')
    wr.write(r'& ')
    wr.write('\n')
    wr.write(r'  ')
    wr.write(r'\multirow{2}{*}{$m$}')
    wr.write(r'& ')
    wr.write('\n')
    wr.write(r'  ')
    wr.write(r'\multicolumn{3}{c|}{Raster length (mm)}')
    wr.write(r' \\')
    wr.write('\n')
    wr.write(r'  ')
    wr.write(r'\cline{9-11}')
    wr.write('\n')
    wr.write(r'  ')
    wr.write(r' & & & & & & & & ')
    wr.write('\n')
    wr.write(r'  ')
    wr.write(r'$\Lmax$')
    wr.write(r'& ')
    wr.write(r'$\Lavg$')
    wr.write(r'& ')
    wr.write(r'$\Ltot$')
    wr.write(r' \\')
    wr.write('\n')

    wr.write(r' ')
    wr.write(r'\hline')
    wr.write('\n')
    wr.write('\n')

  return

# --------------------------------------------------------------- #

def write_table(file_out, details_dict, df_dict, dict_scanline):
  with open(file_out, 'a') as wr:
    for index in range(1, len(df_dict.keys())):
      key_part  = df_dict[index]['part']
      key_angle = df_dict[index]['angle']

      wr.write(r'  ')
      wr.write(str(index))
      wr.write(r' & ')
      wr.write('\n')

      wr.write(r'  ')
      wr.write(r'\multirow{1}{*}{\sfo{' + key_part + r'}{' + key_angle +  r'}}')
      wr.write(r' & ')
      wr.write('\n')

      wr.write(r'  ')
      wr.write(r'\multirow{1}{*}{' +  details_dict[key_part]['islice'] + r'}')
      wr.write(r' & ')
      wr.write('\n')

      wr.write(r'  ')
      wr.write(r'$' + key_angle + r'$')
      wr.write(r' & ')
      wr.write(str(dict_scanline[str(key_part)][str(key_angle)]['x']))
      wr.write(r' & ')
      wr.write(str(dict_scanline[str(key_part)][str(key_angle)]['y']))
      wr.write(r' & ')
      wr.write(str(details_dict[key_part][key_angle]['n_ras']))
      wr.write(r' & ')
      wr.write('\n')

      wr.write(r'  ')
      wr.write(r'$' + str(dict_scanline[str(key_part)][str(key_angle)]['nsc']) + r'$')
      wr.write(r' & ')
      wr.write(str(details_dict[key_part][key_angle]['max_l']))
      wr.write(r' & ')
      wr.write(str(details_dict[key_part][key_angle]['avg_l']))
      wr.write(r' & ')
      wr.write(str(details_dict[key_part][key_angle]['tot_l']))
      wr.write(r' \\')

      wr.write('\n')
      wr.write('\n')

  return

# --------------------------------------------------------------- #

def write_end_table(file_out):
  with open(file_out, 'a') as wr:
    wr.write(r' ')
    wr.write(r'\hline')
    wr.write('\n')
    wr.write(r'\end{tabular}')
    wr.write('\n')
  return

# --------------------------------------------------------------- #
def create_input_table(df_dict, dict_scanline):
  details_dict = dict()

  with open(details_file, 'r') as f:
    for line in f:
      line = line.replace('\n', '')
      line = line.split(';')
      details_dict[line[1]] = dict()
      details_dict[line[1]]['part'] = line[0]
      details_dict[line[1]]['website'] = line[2]
      details_dict[line[1]]['date'] = line[3]

  # --------------------------------------------------------------- #

  for file_in in os.listdir(folder_in):
    line = file_in.replace('.txt', '')
    line = line.split('_')
    part_key = line[0]
    islice = line[1]
    angle = line[2]

    details_dict[part_key]['islice'] = islice
    details_dict[part_key][angle] = dict()

    min_l, max_l, tot_l, n_ras = get_info(folder_in, file_in)

    details_dict[part_key][angle]['min_l'] = round(min_l, 1)
    details_dict[part_key][angle]['max_l'] = round(max_l, 1)
    details_dict[part_key][angle]['tot_l'] = round(tot_l, 1)
    details_dict[part_key][angle]['avg_l'] = round(tot_l/n_ras, 1)
    details_dict[part_key][angle]['n_ras'] = n_ras

  # --------------------------------------------------------------- #

  write_header(file_out)
  write_table(file_out, details_dict, df_dict, dict_scanline)
  write_end_table(file_out)

  return
