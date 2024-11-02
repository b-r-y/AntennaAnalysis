import sys

sys.path.insert(0, "C:\\Users\\boyan\\git\\AntennaAnalysis")

import numpy as np

import antenna_analysis
from antenna_analysis.data_load_store import get_analytical_parameters
from antenna_analysis.imports.pattern_creator import pattern_creator
from antenna_analysis.imports.tgpp_antenna import TGPPAntenna

# from antenna_analysis.imports.tgpp_antenna import TGPPAntenna

# Initialize an object from the class
DATA = antenna_analysis.Patterns()

input_parameters = {
    "ant_type": "3gpp antenna",
    "resolution": [5, 5],
    # "frequency": 8.25e9,
    # "efficiency": 70,
    # "substrate_height": 1,
    # "epsilon_r": 3.66,
    # 'patch_radius': 86
}

input_parameters = get_analytical_parameters(input_parameters)
tgpp_ant: TGPPAntenna = pattern_creator.get("3GPP Antenna", **input_parameters)

# print(patch.describe())

print(tgpp_ant.wavelength0)

# Generate analytical circular patch
DATA.generate_analytical_pattern(input_parameters)

eff, minf, maxf, avgf = DATA.get_field_stats(0)

print(f"Maximum gain: {10*np.log10(maxf)}")
print(f"Efficiency: {10*np.log10(eff)}")

FIELD_TO_ANALYZE = "Gabs"
PATTERNS_TO_MANIPULATE = [0]

PLOT_FIELD, PLOT_PHASE, LEGENDS = DATA.prepare_for_plotting(
    PATTERNS_TO_MANIPULATE, FIELD_TO_ANALYZE
)
# Creates the plotting object
PLOTS = antenna_analysis.PatPlot(PLOT_FIELD, LEGENDS)

PLOTS.make_3d_plot(0)

# Plot some linear stuff
PLOTS.set_plot_params(
    {
        "dynamic_range": ["auto", 40, "log"],
        "plot_range": [[0, 180], [0, 360]],
        "plot_type": "polar",
        "output_path_pad": "2d_linear",
    }
)
PLOTS.make_standard_plot(data_select="XZ")

# Complete the script
print("\n")
print("Done")
