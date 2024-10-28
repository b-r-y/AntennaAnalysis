"""Example usage of the pyramidal horn."""
import numpy as np

import antenna_analysis
from antenna_analysis.imports.pattern_creator import pattern_creator
from antenna_analysis.data_load_store import get_analytical_parameters
from antenna_analysis.imports.pyramidal_horn import PyramidalHorn

# Initialize an object from the class
DATA = antenna_analysis.Patterns()

input_parameters = {'ant_type': 'pyramidal horn',
                    'resolution': [1, 1],
                    'frequency': 24.5E9,
                    'efficiency': 70,
                    'target_gain': 15,
                    'a': 0.0107,
                    'b': 0.0043,
                    # 'a1': 0.05619927,
                    # 'b1': 0.043181136,
                    # 'rho1': 0.071445,
                    # 'rho2': 0.07945898,
                    }

DATA.generate_analytical_pattern(input_parameters)

input_parameters = get_analytical_parameters(input_parameters)
horn: PyramidalHorn = pattern_creator.get('Pyramidal Horn', **input_parameters)
D_E, D_H, D_P = horn.calculate_approximate_directivity()
print(f'Expected directivity:\n'
      f'DE: {10*np.log10(D_E)} dB\n'
      f'DH: {10*np.log10(D_H)} dB\n'
      f'DP: {10*np.log10(D_P)} dB\n'
      )
print(horn.export_mechanical_design_parameters())

# Print the list of loaded Patterns
PATT_LIST = DATA.list_patterns()
print('Patterns List:')
print("\n".join(PATT_LIST))
print('\n')

eff, _, max_gain, _ = DATA.get_field_stats(pat_ind=0, field='Gabs')

print(10*np.log10(eff))
print(10*np.log10(max_gain))

FIELD_TO_ANALYZE = 'Gabs'
PATTERNS_TO_MANIPULATE = [0]

PLOT_FIELD, PLOT_PHASE, LEGENDS = \
      DATA.prepare_for_plotting(PATTERNS_TO_MANIPULATE, FIELD_TO_ANALYZE)
# Creates the plotting object
PLOTS = antenna_analysis.PatPlot(PLOT_FIELD, LEGENDS)

# # Plot some linear stuff
PLOTS.set_plot_params({'dynamic_range': ['auto', 40, 'log'],
                       'plot_range': [[0, 180], [0, 360]],
                       'plot_type': 'polar',
                       'output_path_pad': '2d_linear'
                       })
PLOTS.make_standard_plot(data_select='XZ')

# Complete the script
print('\n')
print('Done')
