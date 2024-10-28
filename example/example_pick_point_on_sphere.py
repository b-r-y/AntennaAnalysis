"""Pick point on sphere example.

In this example it will shown how to load a pattern and fetch a point on
the sphere of that pattern for various supported fields.
"""
from antenna_analysis.patterns import Patterns

DATA = Patterns()

DATA.generate_analytical_pattern()  # generates half lambda dipole

field, phase = DATA.get_point_on_sphere(pat_ind=0, direction=[90, 180])

print(f'Field value: {field}')
print(f'Phase value: {phase}')
