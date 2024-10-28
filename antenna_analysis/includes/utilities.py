"""Created on Nov 9, 2018.

@author: bry
"""
import sys
from typing import Any, TypedDict
import numpy as np
from numpy.typing import NDArray, ArrayLike
from numpy import complex64, inf


class PatternType(TypedDict, total=False):
    """Pattern structure.

    Arguments:
        source (str): a string typically containing the folder the imported \
            filed were stored
        file (str): a string with the file name of the pattern(s) imported
        frequency (float): the pattern frequency. This is extracted from the \
            file name in case of CST files
        port (str): antenna port identification - currently unused
        th_step (int): Theta angle stepping for the pattern resolution. \
            This is always integer and constant for the entire pattern.
        ph_step (int): Phi angle stepping for the pattern resolution. \
            This is always integer and constant for the entire pattern.
        ff (array): ff stand for far frield. The `ff` stores the electric \
            field in E_theta and E_Phi components. The dimensions of `ff` \
            are 2 x n_theta x n_phi where n_theta/n_phi are the total number \
            of samples in Theta and Phi, spherical coordinate system as\
            defined in [IEEESTDTest]_. E_theta is stores as index 0 and E_phi \
            as index 1.
    """

    source: str
    file: str
    frequency: float
    port: int
    th_step: float
    ph_step: float
    rot_offset: list[float]
    ff: NDArray[np.complex64]


def get_p_struct() -> PatternType:
    """Initialize a new empty ``patterns`` list member.

    This ordered dictionary describes how patterns are stored internally and
    what information is kept at import.

    Arguments:
        source (str): a string typically containing the folder the \
            imported filed were stored
        file (str): a string with the file name of the pattern(s) imported
        frequency (float): the pattern frequency. This is extracted \
            from the file name in case of CST files
        port (str): antenna port identification - currently unused
        th_step (int): Theta angle stepping for the pattern resolution. \
            This is always integer and constant for the entire pattern.
        ph_step (int): Phi angle stepping for the pattern resolution. \
            This is always integer and constant for the entire pattern.
        ff (NDArray[np.complex64]): ff stand for far frield. \
            The `ff` stores the electric \
            field in E_theta and E_Phi components. \
            The dimensions of `ff` are 2 x n_theta x n_phi, where \
            n_theta/n_phi are the total number of samples in Theta and Phi, \
            spherical coordinate system as defined in [IEEESTDTest]_. \
            E_theta is stores as index 0 and E_phi as index 1.
    """
    p_struct: PatternType = PatternType()
    p_struct['source'] = ""
    p_struct['file'] = ""
    p_struct['frequency'] = 0
    p_struct['port'] = 0
    p_struct['th_step'] = 0
    p_struct['ph_step'] = 0
    p_struct['rot_offset'] = [0, 0, 0]
    # Default Far Field stored as Real/Imaginary Electric Fields in Theta/Phi
    # Linear polarization
    p_struct['ff'] = np.array([0], dtype=complex64)
    return p_struct


def lin2log(field: str, field_data: NDArray[np.float64]) \
        -> NDArray[np.float64]:
    """Convert the ``field_data`` of ``field`` from linear to logarithmic scale.

    Arguments:
        field (str): the field being converted. Options are the same as in \
            ``fetch_field``
        field_data (NDArray[np.complex64]): a matrix with the field \
            values along Theta and Phi

    Returns:
        new_field_data (ArrayLike): ``field_data`` converted to logarithmic \
            scale depending on the ``field``

    """
    new_field_data = np.zeros(field_data.shape)
    if any(field in s for s in ['Gabs', 'Gth', 'Gph', 'GL', 'GR', 'Eabs']):
        new_field_data = 10*np.log10(field_data)
    elif any(field in s for s in ['AR', 'Eth', 'Eph', 'El', 'ER']):
        new_field_data = 20*np.log10(np.abs(field_data))
    new_field_data[new_field_data < -100.0] = -100.0
    return new_field_data


def close_sphere(matrix):
    """Close the sphere for plotting.

    Used for interpolation and other computations
    by copying the first column and appending it as last.
    """
    matrix = np.append(matrix, np.transpose([matrix[:, 0]]), axis=1)
    return matrix


def open_sphere(matrix):
    """Remove the 360th degree column in matrix representation of a field.

    Used for example when statistics are to be calculated for data
    prepared for plotting.
    """
    matrix = matrix[:, 0:-1]
    return matrix


def apply_analysis_range(pat, analysis_range=None):
    """Remove the data points outside the ``analysis_range`` limits.

    It expects an open field structure.

    Arguments:
        pat (array): a numpy array with the field. \
            The field must be of open type, meaning that the\
            360th phi degree column is not present.
        analysis_range ([[int,int],[int,int]]): angular range to analyze, \
            where the first two\
            values are theta start and stop angles and the second two \
            values are phi start and stop angles.
    """
    if analysis_range is None:
        analysis_range = [[0, 180], [0, 360]]
    else:
        check_angles_input(analysis_range, 'theta/phi')
    th_step, ph_step = get_angle_steps(pat, sphere_type='open')

    num_points = np.prod(pat.shape)//2
    precision = 10000

    rng = np.random.default_rng()
    pick_th = rng.integers(precision, size=num_points, endpoint=True)/precision
    pick_th = ((np.arccos(2*pick_th - 1))*180/np.pi)//th_step
    remove_th_points = (pick_th < analysis_range[0][0]) | \
        (pick_th > analysis_range[0][1])

    pick_ph = rng.integers(precision, size=num_points, endpoint=True)/precision
    pick_ph = ((2*np.pi*pick_ph)*180/np.pi)//ph_step
    remove_ph_points = (pick_ph < analysis_range[1][0]) | \
        (pick_ph > analysis_range[1][1])

    pick_th = np.delete(pick_th, remove_th_points | remove_ph_points)
    pick_ph = np.delete(pick_ph, remove_th_points | remove_ph_points)
    pick_ph = np.where(pick_ph == 360 // ph_step, 0, pick_ph)

    selected_points = []
    for ind, picked_ph_val in enumerate(pick_ph):
        selected_points.append(pat[int(pick_th[ind]), int(picked_ph_val)])

    return selected_points


def get_angle_steps(pattern, sphere_type=None):
    """Return the theta and phi steps based on the pattern shape."""
    no_samples = pattern.shape
    th_step, ph_step = get_angle_steps_from_shape(no_samples, sphere_type)

    return int(th_step), int(ph_step)


def get_angle_steps_from_shape(shape, sphere_type=None):
    """Return the theta and phi steps based on ``shape`` input list."""
    if sphere_type is None:
        sys.exit('Unknown sphere type - this MUST be explicitly defined!')
    th_step = (180/(shape[0]-1))
    ph_step = 0
    if sphere_type == 'closed':
        ph_step = (360/(shape[1]-1))
    elif sphere_type == 'open':
        ph_step = (360/(shape[1]))
    return int(th_step), int(ph_step)


def get_angle_vectors(pattern, sphere_type=None, return_mesh=True):
    """Calculate the theta and phi angle vectors and meshes.

    Used for integration, interpolation and plotting.
    """
    if sphere_type is None:
        sys.exit('Unknown sphere type - this MUST be explicitly defined!')
    resolution = get_angle_steps(pattern, sphere_type)

    # Theta angle is always 0 to 180 including both poles
    vth = np.array(range(0, 180 + resolution[0], resolution[0])) \
        * np.pi / 180
    vph = 0

    # The Phi angle is most of the time open so that phi is between
    # 0 and 360-phi_step.
    # The full 0 to 360 degree phi is used for plotting and integration.
    if sphere_type == 'closed':
        vph = np.array(range(0, 360 + resolution[1], resolution[1]))\
                                         * np.pi / 180
    elif sphere_type == 'open':
        vph = np.array(range(0, 360, resolution[1]))\
                                         * np.pi / 180

    if return_mesh:
        vth, vph = np.meshgrid(vth, vph)

    return vth, vph


def update_params(ref_values: dict[str, Any], input_values: dict[str, Any]) \
        -> dict[str, Any]:
    """Update a default dictionary with inputs from new dictionary.

    Used to keep track and pass parameters.
    """
    # Assign inputs to replace defaults
    if input_values:
        for key, value in input_values.items():
            if key in ref_values:
                if key == 'plot_range':
                    check_angles_input(value, 'theta/phi')
                ref_values[key] = value
            else:
                sys.exit(f'Unknown input values requested! --> {key}')
    return ref_values


def fix_singularities(matrix: NDArray[np.float64]) -> NDArray[np.float64]:
    """Remove singularities from submitted matrix.

    When handling zeros in patterns with partially measured sphere or
    analytical idealistic patterns it is necessary to correct some
    singularities resulting from division by zero or zero divided by zero.
    Those result in NaN, -INF or +INF. In this function these are corrected
    in LINEAR domain so that NaN becomes zero, -INF becomes -100000
    (equivalent to -100 dB) and +INF becomes 100000 (equivalent to +100 dB)
    """
    matrix[matrix == -inf] = -100000
    matrix[matrix == inf] = 100000
    matrix[np.isnan(matrix) == 1] = 0
    return matrix


def check_angles_input(angles, angle_type):
    """Check angles input.

    Plot and analysis ranges for theta and phi are defined between 0 and
    180 and 0 and 360 correspondingly.
    This function checks that the inputs are within these ranges.
    """
    if angle_type == 'theta/phi':
        if not all([0 <= ang <= 180 for ang in angles[0]]):
            sys.exit(f'Theta angle is defined between 0 and 180 degrees! '
                     f'Input was {angles[0]}')
        if not all([0 <= ang <= 360 for ang in angles[1]]):
            sys.exit(f'Phi angle is defined between 0 and 360 degrees! '
                     f'Input was {angles[1]}')
    elif angle_type == 'theta/phi direction':
        if not all([0 <= angles[0] <= 180]):
            sys.exit(f'Theta angle is defined between 0 and 180 degrees! '
                     f'Input was {angles[0]}')
        if not all([0 <= angles[1] <= 180]):
            sys.exit(f'Phi angle is defined between 0 and 360 degrees! '
                     f'Input was {angles[1]}')
