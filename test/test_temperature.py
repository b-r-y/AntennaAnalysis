import sys

sys.path.insert(0, "C:\\Users\\boyan\\git\\AntennaAnalysis")

import antenna_analysis as aa

brightness_model = {
    "type": "cone",
    "temperature": 290,
    "diameter": 70,
    "direction": [0, 0],
}
t_brightness = aa.includes.field_operations.generate_brightness_temperature(
    model=brightness_model
)

print(t_brightness)
