'''
Created on Nov 9, 2018

@author: bry
'''
import sys
import numpy as np
from numpy import inf


def lin2log(field, field_data):
    '''
    Converts the ``field_data`` of ``field`` from linear to logarithmic scale.

    Arguments:
        field (str): the field being converted. Options are the same as in ``fetch_field``
        field_data (num,num): a matrix with the field values along Theta and Phi

    Returns:
        float: ``field_data`` converted to logarithmic scale depending on the ``field``

    '''
    if any(field in s for s in ['Gabs', 'Gth', 'Gph', 'GL', 'GR', 'Eabs']):
        with np.errstate(divide='ignore'):
            field_data = 10*np.log10(field_data)
    elif any(field in s for s in ['AR', 'Eth', 'Eph', 'El', 'ER']):
        with np.errstate(divide='ignore'):
            field_data = 20*np.log10(np.abs(field_data))
    field_data[field_data == -inf] = -100
    return field_data


def close_sphere(matrix):
    '''
    This method closes the sphere for plotting, interpolation and other computations by copying\
    the first column and appending it as last.
    '''
    matrix = np.append(matrix, np.transpose([matrix[:, 0]]), axis=1)
    return matrix

def open_sphere(matrix):
    '''
    Generic method to remove the 360th degree column in matrix representation of a field. This is
    used for example when statistics are to be calculated for data prepared for plotting.
    '''
    matrix = matrix[:,0:-1]
    return matrix

def apply_analysis_range(pat, analysis_range = None):
    '''
    Applies a filter removing the data points outside the ``analysis_range`` limits. It expects an
    open filed structure.
    Arguments:
        pat (array): a numpy array with the field. The field must be of open type, meaning that the\
            360th phi degree column is not present.
        analysis_range ([[int,int],[int,int]]): angular range to analyze where the first two\
            values are theta start and stop angles and the second two values are phi start and stop\
            angles.
    '''
    if analysis_range is None:
        analysis_range = [[0, 180], [0, 360]]
    else:
        check_angles_input(analysis_range, 'theta/phi')
    th_step, ph_step = get_angle_steps(pat, sphere_type = 'open')

    if analysis_range[0][0]==analysis_range[0][1]:
        analysis_range[0][1] = analysis_range[0][1] + th_step
    if analysis_range[1][0]==analysis_range[1][1]:
        analysis_range[1][1] = analysis_range[1][1] + ph_step

    pat = pat[analysis_range[0][0]//th_step:\
                               (analysis_range[0][1]+th_step)//th_step , \
                               analysis_range[1][0]//ph_step:(analysis_range[1][1])//ph_step]

    return pat

def get_angle_steps(pattern, sphere_type = None):
    '''
    Returns the theta and phi steps based on the pattern shape.
    '''
    no_samples = pattern.shape
    th_step, ph_step = get_angle_steps_from_shape(no_samples, sphere_type)

    return int(th_step), int(ph_step)

def get_angle_steps_from_shape(shape, sphere_type=None):
    '''
    Returns the theta and phi steps based on ``shape`` input list.
    '''
    if sphere_type is None:
        sys.exit('Unknown sphere type - this MUST be xplicitly defined!')
    th_step = (180/(shape[0]-1))
    if sphere_type == 'closed':
        ph_step = (360/(shape[1]-1))
    elif sphere_type == 'open':
        ph_step = (360/(shape[1]))
    return int(th_step), int(ph_step)

def get_angle_vectors(pattern, sphere_type = None, return_mesh = True):
    '''
    Calculates the theta and phi angle vectors and meshes for integration, interpolation and
    plotting.
    '''
    if sphere_type is None:
        sys.exit('Unknown sphere type - this MUST be xplicitly defined!')
    resolution = get_angle_steps(pattern, sphere_type)

    # Theta angle is always 0 to 180 including both poles
    vth = np.array(range(0, 180 + resolution[0], resolution[1]))\
                                         * np.pi / 180

    # The Phi angle is most of the time open so that phi is between 0 and 360-phi_step
    # the full 0 to 360 degree phi is used for plotting and integration.
    if sphere_type == 'closed':
        vph = np.array(range(0, 360 + resolution[1], resolution[1]))\
                                         * np.pi / 180
    elif sphere_type == 'open':
        vph = np.array(range(0, 360 , resolution[1]))\
                                         * np.pi / 180

    if return_mesh:
        vth, vph = np.meshgrid(vth, vph)

    return vth, vph

def update_params(ref_values, input_values):
    '''
    Updates a default dictionary with inputs from new dictionary. Used to keep track and pass
    parameters.
    '''
    # Assign inputs to replace defaults
    if input_values is not None:
        if isinstance(input_values,dict):
            for key, value in input_values.items():
                if key in ref_values:
                    if key == 'plot_range':
                        check_angles_input(value, 'theta/phi')
                    ref_values[key] = value
                else:
                    sys.exit('Unknown input values requested! --> {}'.format(key))
        else:
            sys.exit('This requires a dict input!!')

    return ref_values

def fix_singularities(matrix):
    '''
    When handling zeros in patterns with partialy measured sphere or analytical idealistic patterns
    it is necesary to correct some singularities resulting from division by zero or zero divided by
    zero. Those result in NaN, -INF or +INF. In this function these are corrected in LINEAR domain
    so that NaN becomes zero, -INF becomes -100000 (equivalent to -100 dB) and +INF becomes 100000
    (equivalent to +100 dB)
    '''
    matrix[matrix == -inf] = -100000
    matrix[matrix == inf] = 100000
    matrix[np.isnan(matrix) == 1] = 0
    return matrix

def check_angles_input(angles,angle_type):
    '''
    Plot and analysis ranges for theta and phi are defined between 0 and 180 and 0 and 360
    correspondingly. This function checks that the inputs are within these ranges.
    '''
    if angle_type == 'theta/phi':
        if not all([0 <= ang <= 180 for ang in angles[0]]):
            sys.exit('Theta angle is defined between 0 and 180 degrees! Input was {}'.\
                     format(str(angles[0])))
        if not all([0 <= ang <= 360 for ang in angles[1]]):
            sys.exit('Phi angle is defined between 0 and 360 degrees! Input was {}'.\
                     format(str(angles[1])))
    elif angle_type == 'theta/phi direction':
        if not all([0 <= angles[0] <= 180]):
            sys.exit('Theta angle is defined between 0 and 180 degrees! Input was {}'.\
                     format(str(angles[0])))
        if not all([0 <= angles[1] <= 180]):
            sys.exit('Phi angle is defined between 0 and 360 degrees! Input was {}'.\
                     format(str(angles[1])))
