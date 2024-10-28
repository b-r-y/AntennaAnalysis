''' Class with the structure of the DATA '''
import numpy as np

import antenna_analysis
import antenna_analysis.data_load_store

# Initialize an object from the class
DATA = antenna_analysis.Patterns()

# Example usage of the analytical pattern generation with defaults - half lambda dipole
DATA.generate_analytical_pattern(input_parameters={"efficiency": 50}) # pattern index 0 with -3 dB efficiency
DATA.normalize_pattern(0) # pattern index 1 with 0 dB efficiency via normalization to total power of 1
DATA.offset_efficiency(0, offset=3)  # pattern index 2 with 0 dB efficiency via offset gain of + 3dB
DATA.offset_efficiency(0, offset=-1.5)  # pattern index 1 with -4.5 dB efficiency via offset gain of -1.5 dB

# Print the list of loaded Patterns
PATT_LIST = DATA.list_patterns()
print('Patterns List:')
print("\n".join(PATT_LIST))
print('\n')

print(10*np.log10(DATA.get_tot_field_power(pat_ind=0, field='Gabs')))
print(10*np.log10(DATA.get_tot_field_power(pat_ind=1, field='Gabs')))
print(10*np.log10(DATA.get_tot_field_power(pat_ind=2, field='Gabs')))
print(10*np.log10(DATA.get_tot_field_power(pat_ind=3, field='Gabs')))

FIELD_TO_ANALYZE = 'Gabs'
PATTERNS_TO_MANIPULATE = [0, 1, 2, 3]

PLOT_FIELD, PLOT_PHASE, LEGENDS = DATA.prepare_for_plotting(PATTERNS_TO_MANIPULATE, \
                                                            FIELD_TO_ANALYZE)
# Creates the plotting object
PLOTS = antenna_analysis.PatPlot(PLOT_FIELD, LEGENDS)

# # Plot some linear stuff
PLOTS.set_plot_params({'dynamic_range': ['auto', 10, 'log'],\
                       'plot_range': [[0,180],[0,360]],\
                       'plot_type': 'linear',\
                       'output_path_pad': '2d_linear' \
                       })
PLOTS.make_standard_plot(data_select = 'XY')

# Complete the script
print('\n')
print('Done')