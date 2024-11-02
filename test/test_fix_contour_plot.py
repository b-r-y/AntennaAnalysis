import sys

sys.path.insert(0, "C:\\Users\\boyan\\git\\AntennaAnalysis")

import numpy as np

import antenna_analysis
import antenna_analysis.data_load_store

from antenna_analysis.imports.pattern_creator import pattern_creator
from antenna_analysis.includes.field_operations import generate_brightness_temperature
from antenna_analysis.includes.utilities import close_sphere

# Initialize an object from the class
DATA = antenna_analysis.Patterns()

input_parameters = {
    "ant_type": "pyramidal horn",
    "resolution": [5, 5],
    "frequency": 26e9,
    "efficiency": 100,
    "target_gain": 30,
    "a": 0.0107,
    "b": 0.0043,
    # 'a1': 0.05619927,
    # 'b1': 0.043181136,
    # 'rho1': 0.071445,
    # 'rho2': 0.07945898,
}

# Example usage of the analytical pattern generation with defaults - half lambda dipole
DATA.generate_analytical_pattern(
    input_parameters
)  # pattern index 0 with -3 dB efficiency
DATA.set_permanent_rotation_offset(0, [0, 90, 0])

input_parameters = antenna_analysis.data_load_store.get_analytical_parameters(
    input_parameters
)
horn = pattern_creator.get("Pyramidal Horn", **input_parameters)
D_E, D_H, D_P = horn.calculate_approximate_directivity()
print(
    f"Expected directivity:\n"
    f"DE: {10*np.log10(D_E)} dB\n"
    f"DH: {10*np.log10(D_H)} dB\n"
    f"DP: {10*np.log10(D_P)} dB\n"
)
print(horn.export_mechanical_design_parameters())

# Print the list of loaded Patterns
PATT_LIST = DATA.list_patterns()
print("Patterns List:")
print("\n".join(PATT_LIST))
print("\n")

eff, _, max_gain, _ = DATA.get_field_stats(pat_ind=0, field="Gabs")

print(f"Efficiency [dB]: {10*np.log10(eff)}")
print(f'Efficiency 2 [dB]: {10*np.log10(DATA.get_tot_field_power(0, "Gabs"))}')
print(f"Peak gain [dBi]: {10*np.log10(max_gain)}")

brightness_model = {
    "type": "cone",
    "temperature": 290,
    "diameter": 70,
    "direction": [0, 0],
}
T = DATA.compute_antenna_temperature([0], brightness_model=brightness_model)

t_brightness = generate_brightness_temperature(model=brightness_model, shape=[37, 72])
t_brightness = close_sphere(t_brightness)

print(T)
print(T[0] * eff)


FIELD_TO_ANALYZE = "Gabs"
PATTERNS_TO_MANIPULATE = [0]

PLOT_FIELD, PLOT_PHASE, LEGENDS = DATA.prepare_for_plotting(
    PATTERNS_TO_MANIPULATE, FIELD_TO_ANALYZE
)

PLOT_FIELD.append(t_brightness)
LEGENDS.append("Brightness")
# Creates the plotting object
PLOTS = antenna_analysis.PatPlot(PLOT_FIELD, LEGENDS)

# Plot same thing but as polar
PLOTS.set_plot_params({"output_path_pad": "2d_polar", "plot_type": "polar"})
# PLOTS.make_standard_plot(data_select="XZ", linear_plot_offset=180)
# PLOTS.make_standard_plot(data_select="YZ", xtick_spacing=45)
# PLOTS.make_standard_plot(data_select=["phi", [30, 90, 180]])

PLOTS.make_countour_plot(0, tick_spacing=[45, 60])
PLOTS.make_countour_plot(1, tick_spacing=[45, 60])

# Complete the script
print("\n")
print("Done")
