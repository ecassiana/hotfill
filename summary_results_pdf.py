import os
import shutil
from pandas.core.accessor import register_index_accessor
import pdfkit
import pandas as pd

# --------------------------------------------------------------- #

kpi_keys = ['FabTime', 'CoolTime', 'AirTime', 'CpuTime']

# --------------------------------------------------------------- #

def create_line(slic3r, rp3, scanline, altscanline, data_key_part, key_part, key_angle, delta_keys):
  # Create and returns a line to be used in the report
  # containing the part name and the information about it.

  line = [key_part, key_angle]

  data_key_part['SLIC3R'] = slic3r[key_part][key_angle]['SLIC3R']
  data_key_part['RP3'] = rp3[key_part][key_angle]['RP3']
  data_key_part['SCANLINE'] = scanline[key_part][key_angle]['SCANLINE']
  data_key_part['ALTSCANLINE'] = altscanline[key_part][key_angle]['ALTSCANLINE']

  tfab_ref = min(float(data_key_part['RP3']['FabTime']), float(data_key_part['SLIC3R']['FabTime']))

  tras = scanline[key_part][key_angle]['SCANLINE']['RastTime']
  line.append(tras)

  sort_field = -999.0

  for key_delta in delta_keys:
    if key_delta in delta_keys[:4]:
      line.append(data_key_part[key_delta]['FabTime'])
      line.append(data_key_part[key_delta]['CoolTime'])
      line.append(data_key_part[key_delta]['AirTime'])

    else:
      for key_value in kpi_keys:
        line.append(data_key_part[key_delta][key_value])
      
      if data_key_part[key_delta]['FabTime'] != '-':
        tfab = float(data_key_part[key_delta]['FabTime'])
        dif = (tfab-tfab_ref)/tfab_ref
        dif = round(100*dif, 1)
        line.append(str(dif))
        if (key_delta == '8.0'):
            sort_field = str(dif)
      else:
        line.append('-')

  line.append(sort_field)
  return line

# --------------------------------------------------------------- #

def create_dataframe_table(slic3r, rp3, scanline, altscanline, result_data, delta_keys):
  # Create and return a dataframe containing the results of the 
  # hotpath_TST executions.

  delta_keys_aux = get_deltas_aux(result_data)

  if len(delta_keys_aux) < 5:
    return None

  headers_name = [
    'Model', 'Deg', 'Tras', 
    'Slic3r', '', '', 
    'RP3', '', '', 
    'Scan', '', '', 
    'Alt', '', ''
  ]

  for delta in delta_keys_aux[4:]:
    headers_name = headers_name + ['Hot ' + str(delta), '', '', '', '']
  
  headers_name = headers_name + ['Sort']

  sub_headers_name_a = ['Tfab', 'Tcool', 'Tair']
  sub_headers_name_b = ['Tfab', 'Tcool', 'Tair', 'Tcpu', 'Dif(%)']
  sub_headers_name = ['Model', 'Deg', 'Tras'] + 4*sub_headers_name_a + len(delta_keys_aux[4:])*sub_headers_name_b + ['-1.0']

  df = pd.DataFrame(columns = headers_name)
  df = df.append(pd.Series(sub_headers_name, index = headers_name), ignore_index = True)
 
  for key_part in result_data.keys():
    for key_angle in result_data[key_part].keys():
      line = create_line(slic3r, rp3, scanline, altscanline, result_data[key_part][key_angle], key_part, key_angle, delta_keys_aux)
      df = df.append(pd.Series(line, index = headers_name), ignore_index = True)
  
  return df

# --------------------------------------------------------------- #

def create_pdf_report(filename, df):
  # If the dataframe isn't empty, create a report and remove the
  # temporary files.

  df.to_csv(filename + '.csv', index = False, header = True, sep = ';')
    
  html_string = df.to_html(classes = 'table_style')
  html_string = html_string.replace('<tr>', '<tr style="text-align: right;">')
  
  with open(filename + '.html', 'w') as f:
    f.write(html_string)
 
  options = {
         'dpi': 240,
         'page-size': 'A2',
         'margin-top': '0.25in',
         'margin-right': '0.25in',
         'margin-bottom': '0.25in',
         'margin-left': '0.25in',
         'encoding': "UTF-8",
         'orientation': 'landscape',
         'custom-header' : [
            ('Accept-Encoding', 'gzip')
         ],
         'no-outline': None
  }

  pdfkit.from_file(filename + '.html', filename + '.pdf', options = options)

  os.remove(filename + '.csv')
  os.remove(filename + '.html')

  return

# --------------------------------------------------------------- #

def create_dict_order(df):
  df_dict = dict()
  index_dict = 0

  for i, row in df.iterrows():
    df_dict[index_dict] = dict()
    df_dict[index_dict]['part'] = row['Model']
    df_dict[index_dict]['angle'] = row['Deg']
    index_dict = index_dict + 1

  return df_dict

# --------------------------------------------------------------- #

def create_tables(folder_out, slic3r, rp3, scanline, altscanline, result_data):
  df_dict = None
  delta_keys = get_deltas(result_data)

  for key in result_data.keys():
    df = create_dataframe_table(slic3r, rp3, scanline, altscanline, result_data[key], delta_keys)
    if isinstance(df, pd.DataFrame):
      filename = folder_out + '/PDF/summary_results_maxband' + key
      df['Sort'] = df['Sort'].astype('float')
      df.sort_values('Sort', inplace=True)
      df.drop('Sort', axis=1, inplace=True)
      if df_dict == None and key == '20':
        df_dict = create_dict_order(df)
      create_pdf_report(filename, df)
  
  return delta_keys, df_dict

# --------------------------------------------------------------- #

def get_deltas(result_data_original): 
  result_data = result_data_original[list(result_data_original.keys())[0]]
  result_data = result_data[list(result_data.keys())[0]]
  result_data = result_data[list(result_data.keys())[0]]

  delta_keys_str = list(result_data.keys())

  deltas_keys_float = [ float(delta) for delta in delta_keys_str ]
  deltas_keys_float.sort()

  delta_keys = [ str(delta) for delta in deltas_keys_float ]

  delta_keys = ['SLIC3R', 'RP3', 'SCANLINE', 'ALTSCANLINE'] + delta_keys
  
  return delta_keys

# --------------------------------------------------------------- #

def get_deltas_aux(result_data_original): 
  result_data = result_data_original[list(result_data_original.keys())[0]]
  delta_keys_aux = list(result_data['0'].keys())

  delta_keys_str = []

  for delta in delta_keys_aux:
    if '.' in delta:
      delta_keys_str.append(delta)


  deltas_keys_float = [ float(delta) for delta in delta_keys_str ]
  deltas_keys_float.sort()

  delta_keys = [ str(delta) for delta in deltas_keys_float ]

  delta_keys = ['SLIC3R', 'RP3', 'SCANLINE', 'ALTSCANLINE'] + delta_keys
  
  return delta_keys

# --------------------------------------------------------------- #