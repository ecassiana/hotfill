import pandas as pd

# --------------------------------------------------------------- #

def create_line_tconn(result_data, key_part, key_angle, key_delta, tfab_rp3, tfab_slic3r):
  if result_data['FabTime'] == '-':
    line = [key_part + '_' + key_angle, key_delta, '-', '-', '-', '-']
  
  else:
    trast = float(result_data['RastTime'])
    tfab_hp = float(result_data['FabTime'])

    line = [key_part + '_' + key_angle, key_delta, str(trast), str(tfab_rp3), str(tfab_slic3r), str(tfab_hp)]

  return line

# --------------------------------------------------------------- #

def create_tconn_dataframe(result_data, rp3, slic3r):
  headers_name = ['Model', 'Delta', 'Trast', 'Tfab_rp3', 'Tfab_slic3r', 'Tfab_hp']
  df = pd.DataFrame(columns = headers_name)

  for key_part in result_data.keys():
    for key_angle in result_data[key_part].keys():
      tfab_rp3 = float(rp3[key_part][key_angle]['RP3']['FabTime'])
      tfab_slic3r = float(slic3r[key_part][key_angle]['SLIC3R']['FabTime'])

      for key_delta in result_data[key_part][key_angle].keys():
        line = create_line_tconn(result_data[key_part][key_angle][key_delta], key_part, key_angle, key_delta, tfab_rp3, tfab_slic3r)
        df = df.append(pd.Series(line, index = headers_name), ignore_index = True)

  return df

# --------------------------------------------------------------- #

def create_plot_data(folder_out, result_data, rp3, slic3r, not_filter):
  for key in result_data.keys():
    df = create_tconn_dataframe(result_data[key], rp3, slic3r)
    if isinstance(df, pd.DataFrame):
      if not_filter:
        filename = folder_out + '/CSV/plot_data_complete_maxband' + key
        df.to_csv(filename + '.csv', index = False, header = False, sep = ' ')

      df_filter = df.loc[(df.Delta != "inf") & (df.Delta != "SLIC3R") & (df.Delta != "RP3") & (df.Delta != "ALTSCANLINE") & (df.Delta != "SCANLINE")]
      filename = folder_out + '/CSV/plot_data_filter_maxband' + key
      df_filter.to_csv(filename + '.csv', index = False, header = False, sep = ' ')

  return

# --------------------------------------------------------------- #