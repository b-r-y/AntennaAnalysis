"""Main example file."""
import antenna_analysis

# Initialize an object from the class
DATA = antenna_analysis.Patterns()

# Example usage of the analytical pattern generation with defaults -
# half lambda dipole
DATA.generate_analytical_pattern()  # this will be pattern index 0

# Example usage of the analytical pattern generation with defaults -
# one lambda dipole
# this will be pattern index 1
DATA.generate_analytical_pattern({'length': 1.5,'resolution': [5,5]})

# Example usage of the analytical pattern generation with defaults -
# pyramidal horn
# this will be pattern index 1
DATA.generate_analytical_pattern({'ant_type': 'pyramidal horn',
                                  'resolution': [1, 1],
                                  'frequency': 1.5E9,
                                  'target_gain': 15})

# Example loading a single CST file from the dipole example
DATA.load_data("CST_File", "example\\example_antennas\\dipole\\Export\\Farfield\\farfield (f=1000) [1].txt")
# this will be pattern index 2

# Example loading all files in a CST export folder from the dipole example
DATA.load_data("CST_Folder", "example\\example_antennas\\dipole\\Export\\Farfield")
# these will become indexes 3, 4 and 5. Index 4 (1GHz) will be identical to index 2 loaded earlier.

###### Before using the next example you must run the parametric sweep in CST ###########

# Example loading all files in a CST parametric sweep folder from the dipole example

# DATA.load_data("CST_Par_Sweep", "example\\dipole")

# The sweep defined in the example dipole file is
# of lambda 1/2, 1 and 2 for the dipole length, each producing three
# Patterns on export. This makes
# a total of 9 additional Patterns appended after index 5.

# Example loading from JSON file
#DATA.load_data('json', "example\\example_json.json")

# Example storing to JSON
antenna_analysis.data_load_store.save_to_json(DATA,
                                              "example\\example_antennas\\example_json.json")

# Print the list of loaded Patterns
PATT_LIST = DATA.list_patterns()
print('Patterns List:')
print("\n".join(PATT_LIST))
print('\n')

# Selecting to fetch DATA from the index 0.
#FIELD,PHASE,CONTEXT = DATA.fetch_field_with_context(0,field='Gabs',field_format='MA')

FIELD_TO_ANALYZE = 'Gabs'
PATTERNS_TO_MANIPULATE = [0,1,2,10,15,20]

PLOT_FIELD, PLOT_PHASE, LEGENDS = DATA.prepare_for_plotting(PATTERNS_TO_MANIPULATE, \
                                                            FIELD_TO_ANALYZE)
# Creates the plotting object
PLOTS = antenna_analysis.PatPlot(PLOT_FIELD, LEGENDS)

PLOTS.set_plot_params({'dynamic_range': ['auto', 40, 'log'],\
                      'plot_range': [[0,180],[0,360]],\
                      'plot_type': 'linear',\
                      'output_path_pad': 'example_contour' \
                      })
PLOTS.make_countour_plot(0, tick_spacing = [45,60])
PLOTS.make_countour_plot(1, linear_plot_offset = [90,0])
#PLOTS.make_countour_plot(2)

PLOTS.set_plot_params({\
                      'output_path_pad': 'example_3D' \
                      })
PLOTS.make_3d_plot(0)
PLOTS.make_3d_plot(1, n_bins=15)

# # Plot some linear stuff
PLOTS.set_plot_params({'dynamic_range': ['auto', 10000, 'lin'],\
                       'plot_range': [[0,180],[0,360]],\
                       'plot_type': 'linear',\
                       'output_path_pad': '2d_linear' \
                       })
PLOTS.make_standard_plot(data_select = 'XZ',linear_plot_offset = 180)
PLOTS.make_standard_plot(data_select = ['theta', [90]],xtick_spacing = 45, \
                             linear_plot_offset = 180)
PLOTS.make_standard_plot(data_select = ['phi', [30, 90, 180]],linear_plot_offset = 180)

# # Plot same thing but as polar
PLOTS.set_plot_params({'output_path_pad': '2d_polar',\
                       'plot_type' : 'polar'\
                       })
PLOTS.make_standard_plot(data_select = 'XZ',linear_plot_offset = 180)
PLOTS.make_standard_plot(data_select = 'YZ',xtick_spacing = 45)
PLOTS.make_standard_plot(data_select = ['phi', [30, 90, 180]])

# Setting the statistics to omit the poles within 10 degrees - makes some sense for dipoles
FIELD_TO_ANALYZE = 'Gabs'
PATTERNS_TO_MANIPULATE = [0,1] # take just the dipoles here to make it simple.

PLOT_FIELD, PLOT_PHASE, LEGENDS = DATA.prepare_for_plotting(PATTERNS_TO_MANIPULATE, \
                                                            FIELD_TO_ANALYZE)
# Creates a new plotting object
PLOTS = antenna_analysis.PatPlot(PLOT_FIELD, LEGENDS)

PLOTS.set_plot_params({'output_path_pad': 'statistics',\
                       'dynamic_range' : ['MaxFull3D', 40, 'log'],\
                      'plot_range' : [[5,175],[0,360]],\
                      })
PLOTS.make_statistics_plot( stat_plot_type = 'CDF')
PLOTS.make_statistics_plot( stat_plot_type = 'CCDF')

# Example of some correlation calculations
CORR_TABLE = DATA.compute_correlation(PATTERNS_TO_MANIPULATE, corr_type = 'complex')
print('Complex correlation --> ',CORR_TABLE[0,:])

CORR_TABLE = DATA.compute_correlation(PATTERNS_TO_MANIPULATE, corr_type = 'envelope')
print('Envelope correlation --> ',CORR_TABLE[0,:])

# Example for group analysis - finds the 500-1500 MHZ patterns and compures stats.
GROUP_NAMES, MEMBERIDX = DATA.get_filtered_groups({'frequency': [500E6,1500E6]})
print(GROUP_NAMES)
print(MEMBERIDX)
GR_TOT_PWR, GR_MIN, GR_MAX, GR_AVG = DATA.get_group_stats(MEMBERIDX,FIELD_TO_ANALYZE)
print('Total Power --> ', GR_TOT_PWR)
print('Field MIN --> ', GR_MIN)
print('Field MAX --> ', GR_MAX)
print('Field AVG --> ', GR_AVG)

## Export to Exel
antenna_analysis.data_load_store.save_to_excel(DATA, [0,1])

# Complete the script
print('\n')
print('Done')
