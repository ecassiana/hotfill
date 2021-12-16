import os
import pandas as pd
import summary_results_read
import summary_results_pdf
import summary_results_plot
import summary_results_tex
import summary_input
import summary_scanlines

# --------------------------------------------------------------- #

deltas_aux = ['SLIC3R', 'RP3', 'SCANLINE', 'ALTSCANLINE']
deltas_cool = ['1.4', '2.0', '2.8', '4.0', '5.7', '8.0', '11.3', '16.0', '22.6', '32.0', '45.3', '64.0']
deltas_fab  = ['1.4', '2.0', '2.8', '4.0', '5.7', '8.0', '11.3', '16.0', '32.0', '64.0']
max_bands = ['10', '20', '30', '40', '50', '60', '70', '80']

# --------------------------------------------------------------- #

folder_in = './tests/outfiltered/'
folder_out = './tests/results/'

# --------------------------------------------------------------- #

data_parts = summary_results_read.read(folder_in, folder_out, deltas_aux)
dict_scanline = summary_scanlines.main(data_parts)

slic3r, rp3, scanline, altscanline, result_data = summary_results_read.create_dict(data_parts, deltas_aux)

delta_keys, df_dict = summary_results_pdf.create_tables(folder_out, slic3r, rp3, scanline, altscanline, result_data)

summary_input.create_input_table(df_dict, dict_scanline)

summary_results_plot.create_plot_data(folder_out, result_data, rp3, slic3r, not_filter = False)

summary_results_tex.create_exp_tcool(folder_out, deltas_aux + deltas_cool, '20', data_parts, df_dict)

summary_results_tex.create_exp_tfab(folder_out, deltas_aux + deltas_fab, '20', data_parts, df_dict)

summary_results_tex.create_exp_tcpu(folder_out, '8.0', max_bands, result_data, df_dict, dict_scanline)