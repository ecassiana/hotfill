import os
import shutil
import pdfkit
import pandas as pd

# --------------------------------------------------------------- #

def create_folders(folder):
  # Check if the {folder} to save the report files exists. 
  # If doesn't exist, create the {folder}.
  # Then, check if the {folder} contains a subfolder to save
  # the figures of the solutions. 
  # If doesn't exist, create the subfolder.

  if not os.path.exists(folder):
    os.makedirs(folder)
  
  if not os.path.exists( folder + '/solutions'):
    os.makedirs(folder + '/solutions')

  if not os.path.exists( folder + '/PDF'):
    os.makedirs(folder + '/PDF')
    
  if not os.path.exists( folder + '/CSV'):
    os.makedirs(folder + '/CSV')
  
  return

# --------------------------------------------------------------- #

def found_keys(folder, test_list):
  # Search for the details about the hotpath_TST execution and returns
  # {key_delta}: delta used in the test.
  # {key_best}: indicate if was the best solution test.
  # {key_part}: name of the part ued in the test. 
  # {key_angle}: angle of the slice used in the test.

  split_list = folder.split('_')
  key_part = split_list[0].split('/')[-1]
  key_angle = split_list[2]
  
  if split_list[3] in test_list:
    key_delta = split_list[3]
    key_maxband = split_list[3]
  else:
    key_delta = split_list[3].replace('delta', '')
    key_maxband = split_list[4].replace('maxband', '')
    
  return key_part, key_angle, key_delta, key_maxband

# --------------------------------------------------------------- #

def result_data_item():
  # Returns a dictionary that will be used to save a function results.
  #

  return {
    'FabTime':  '-',
    'ExtTime':  '-',
    'AirTime':  '-',
    'CpuTime':  '-',
    'CoolTime': '-',
    'Nrast':    '-',
    'RastTime': '-',
    'Nlink':    '-',
    'LnkTime':  '-',
    'Njump':    '-',
    'JumpTime': '-'
  }

# --------------------------------------------------------------- #

def create_folder_name(folder_in):
  # Finds the name and create, if doesn't exist, a folder to save 
  # the solutions figures.
  # Returns {folder_out_solution}, a string that represents the name
  # of the folder.

  folder_out_solution = './tests/results/solutions/'
  if 'delta' in folder_in:
    folder = folder_in.split('delta')
    folder[-1] = 'delta' + folder[-1]
  else:
    folder = folder_in.split('_')
  folder_out_solution = folder_out_solution + folder[-1]

  if not os.path.exists(folder_out_solution):
    os.makedirs(folder_out_solution)
  
  return folder_out_solution

# --------------------------------------------------------------- #

def read_txt_result(filename, result_data):
  # Open and read a file containing the hotpath_TST results.
  # Update the {result_data} dictionary with the information
  # about the results from th execution.

  with open(filename, 'r') as txt_result:
    for line in txt_result:
      if ':' in line:
        key, value = line.replace('\n', '').split(':')
        if key in result_data.keys():
          if '.' in value:
            if key == 'CpuTime':  
              value = round(float(value), 2)
              value = str("{:0.0f}".format(value))
            else:  
              value = round(float(value), 2)
              value = str("{:0.1f}".format(value))

          result_data[key] = value
  
  return

# --------------------------------------------------------------- #

def copy_file(folder_in, folder_out, filename):
  # Copy the {filename} from {folder_in} to {folder_out}.

  folder_in = folder_in + '/' + filename

  if not os.path.exists(folder_out):
    os.makedirs(folder_out)

  folder_out = folder_out + '/' + filename  

  shutil.copy2(folder_in, folder_out)

  return

# --------------------------------------------------------------- #

def read_result(folder_in):
  # Search the folder for the result file and the solution figure. 
  # Get the information about the results in the result file and
  # copy the solution figure to the solution folder.

  result_data = result_data_item()
  folder_out = create_folder_name(folder_in)
  
  for filename in os.listdir(folder_in):
    if 'out.txt' in filename:
      filename = folder_in + '/' + filename
      read_txt_result(filename, result_data)

    elif '.png' in filename:
      copy_file(folder_in, folder_out, filename)

  return result_data

# --------------------------------------------------------------- #

def create_dict(data_parts, test_list):
  slic3r = data_parts['SLIC3R']
  rp3 = data_parts['RP3']
  scanline = data_parts['SCANLINE']
  altscanline = data_parts['ALTSCANLINE']
  
  result_data = dict()

  for key in data_parts.keys():
    if key not in test_list:
      result_data[key] = data_parts[key]

  return slic3r, rp3, scanline, altscanline, result_data

# --------------------------------------------------------------- #

def read(folder_in, folder_out, deltas_aux):
  data_parts = dict()
  create_folders(folder_out)

  lst = os.listdir(folder_in)
  lst.sort()

  for folder in lst:
    if '.txt' not in folder:
      folder = folder_in + folder
      key_part, key_angle, key_delta, key_maxband = found_keys(folder, deltas_aux)

      if key_maxband not in data_parts.keys():
        data_parts[key_maxband] = dict()

      if key_part not in data_parts[key_maxband].keys():
        data_parts[key_maxband][key_part] = dict()

      if key_angle not in data_parts[key_maxband][key_part].keys():
        data_parts[key_maxband][key_part][key_angle] = dict()
      
      result_data = read_result(folder)

      data_parts[key_maxband][key_part][key_angle][key_delta] = result_data

  return data_parts

# --------------------------------------------------------------- #